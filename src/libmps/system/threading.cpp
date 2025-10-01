/*
* This file is part of MPSolve 3.2.2
*
* Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
* License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
*
* Authors:
*   Dario Andrea Bini <bini@dm.unipi.it>
*   Giuseppe Fiorentino <fiorent@dm.unipi.it>
*   Leonardo Robol <leonardo.robol@unipi.it>
*/


#include <float.h>
#include <mps/mps.h>
#include <stdio.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef __WINDOWS
#include <windows.h>
#endif

/**
* @brief This is the default shared thread pool that contexts of MPSolve use.
* It's shared here to improve the performance of more contexts using it.
*/
static mps_thread_pool* system_thread_pool = NULL;

/**
* @brief This lock is used when trying to detect if the thread pool
* has already been allocated or not.
*/
static mps_initialized_mutex(system_thread_pool_lock, PTHREAD_MUTEX_INITIALIZER);

/**
* @brief Get number of logic cores on the local machine, or
* 0 if that information is not available with the method
* known to this implementations.
*/
MPS_PRIVATE int
    mps_thread_get_code_number(mps_context* ctx)
{
    int cores = 0;
    char* cores_env = NULL;

#ifdef __WINDOWS
    SYSTEM_INFO windows_sys_info;
#endif

    if ((cores_env = getenv("MPS_JOBS")) != NULL)
    {
        /* Give reasonable bounds to the possible values of MPS_JOBS */
        cores = MAX(1, MIN(MPS_MAX_CORES, atoi(cores_env)));

        return cores;
    }

    /* Test for POSIX platforms */
#ifdef HAVE_SYSCONF
    cores = sysconf(_SC_NPROCESSORS_ONLN);
#endif

#ifdef __WINDOWS
    GetSystemInfo(&windows_sys_info);
    cores = windows_sys_info.dwNumberOfProcessors;
#endif

    if (cores != 0)
    {
        MPS_DEBUG_WITH_INFO(ctx, "Found %d cores on this system", cores);
    }

    /* In case no runtime method of finding the available cores
    * worked out, select a fixed value. */
    if (cores <= 0)
    {
        cores = 8;
        if (ctx->debug_level & MPS_DEBUG_INFO)
        {
            MPS_DEBUG(ctx, "No runtime information about available cores found");
            MPS_DEBUG(ctx, "Selecting a fixed number of %d threads", cores);
            MPS_DEBUG(ctx, "Use the MPS_JOBS environment variable to override this value");
        }
    }

    return cores;
}


/**
* @brief Create a new mps_thread_job_queue that can
* handle at most max_iter iterations for n_roots roots.
*/
mps_thread_job_queue*
    mps_thread_job_queue_new(mps_context* ctx)
{
    /* Space allocation and related jobs */

    mps_new_obj(mps_thread_job_queue, q, sizeof(mps_thread_job_queue));

    mps_mutex_init(q->job_queue_mutex);

    /* Set initial data */
    q->iter = 0;
    q->n_roots = ctx->n;
    q->max_iter = ctx->max_it;
    q->cluster_item = ctx->clusterization->first_cluster_item;
    q->root = q->cluster_item->cluster->first_root;
    return q;
}

/*
* @brief Free a mps_thread_job_queue previously allocated
* with mps_thread_job_queue_new ().
*/
void
    mps_thread_job_queue_free(mps_thread_job_queue* q)
{
    mps_mutex_destroy(q->job_queue_mutex);
    mps_del_obj(q);
}

/**
* @brief Obtain iter and i for the next available job.
*/
mps_thread_job
    mps_thread_job_queue_next(mps_context* ctx, mps_thread_job_queue* q)
{
    mps_thread_job j;

    mps_mutex_lock(q->job_queue_mutex);

    j.i = 0;
    j.cluster_item = NULL;

    if (q->iter == MPS_THREAD_JOB_EXCEP)
    {
        j.iter = MPS_THREAD_JOB_EXCEP;
    }
    else
    {
        /* Assigning the root */
        j.i = q->root->k;
        j.cluster_item = q->cluster_item;
        j.iter = q->iter;

        /* Get the next element of the cluster, incrementing the queue */
        q->root = q->root->next_root;

        /* Check if the previous one was the last element in the
        * cluster, and if that's the case pass to the next one. */
        if (q->root == NULL)
        {
            q->cluster_item = q->cluster_item->next_cluster_item;

            /* If we got to the end of the clusterization restart from
            * the first cluster and dump the iteration counter. */
            if (q->cluster_item == NULL)
            {
                q->cluster_item = ctx->clusterization->first_cluster_item;
                q->iter++;
            }

            q->root = q->cluster_item->cluster->first_root;


            /* Check if maximum number of iteration was reached and
            * if that was the case set j->iter to MPS_THREAD_JOB_EXCEP.  */
            if (j.iter == q->max_iter)
            {
                j.iter = MPS_THREAD_JOB_EXCEP;
                q->iter = MPS_THREAD_JOB_EXCEP;
            }
        }
    }

    mps_mutex_unlock(q->job_queue_mutex);
    return j;
}

// a worker thread's main loop - it loops until it is no longer alive 
// each loop either picks up a new task and does it to completion, or wait for a new task to appear in the pool's queue

MPS_PRIVATE void*
    mps_thread_mainloop(void* caller_thread)
{
    mps_thread* thread = (mps_thread*)caller_thread;
    mps_thread_pool* pool = thread->pool;
    thread->sys_thread_id = GetCurrentThreadId();

    while (thread->alive)
    {
        /* Try to pop a work item from the queue, if available. */
        mps_mutex_lock(pool->work_completed_mutex);
        mps_mutex_lock(pool->queue_changed_mutex);

        if (pool->queue->first_item != NULL)
        {
            // there is a job in the queue that has not been assigned
            mps_thread_pool_queue_item* item = pool->queue->first_item;

            if (!thread->busy)  // thread might already have been busy,but if this is the 1st job it is picking up
            {
                pool->busy_counter++;
                thread->busy = true;
            }

            /* Pop item from the queue and release the lock on it */
            pool->queue->first_item = item->next_item;
            if (pool->queue->first_item == NULL)  // queue is suddenly empty
            {
                pool->queue->last_item = item;
                // mps_thread_pool_wait wants to know that queue was just emptied
                mps_cond_signal(pool->work_completed_cond, pool->wait_work_completed_mutex);
            }
            mps_mutex_unlock(pool->queue_changed_mutex);
            mps_mutex_unlock(pool->work_completed_mutex);

            item->work(item->work_args);   // do the work then return

            mps_del_obj(item);
        }
        else  // this pool's job queue is empty
        {
            /* Check if other threads are sleeping */
            if (thread->busy)  // if I just finished a job, I am not now busy
            {
                pool->busy_counter--;
                thread->busy = false;
            }

            // signal to a thread in the pool mps_thread_pool_wait that this thread is available but unassigned
            // the pool itself may not be busy
            mps_cond_signal(pool->work_completed_cond, pool->wait_work_completed_mutex);
            mps_mutex_unlock(pool->work_completed_mutex);

            if (!thread->alive)  // have I just been told to die?
            {
                mps_mutex_unlock(pool->queue_changed_mutex);
                break;
            }

            // this thread is still alive, but had no work.  Sleep until I have been told to die or there is a job in the pool queue
            // wait until this thread is told to stop, or there is a new task in the queue
            mps_cond_wait(pool->queue_changed_cond, pool->wait_queue_changed_mutex, pool->queue_changed_mutex, !thread->alive || pool->queue->first_item != NULL);

            mps_mutex_unlock(pool->queue_changed_mutex);
            // the return could be "spurious" so repeat the loop
        }
    }

    mps_thread_processing_exit(thread->sys_thread);
    return NULL;
}

/**
 * @brief Start the thread mainloop.
 */
MPS_PRIVATE void
mps_thread_start_mainloop (mps_context * ctx, mps_thread * thread)
{
    mps_thread_create(thread->sys_thread, &mps_thread_mainloop, thread)
}

/**
* @brief Limit the maximum number of threads that can be used in the thread pool.
*/
void mps_thread_pool_set_concurrency_limit(mps_context* ctx, mps_thread_pool* pool,
    unsigned int concurrency_limit)
{
    if (pool == NULL)
        pool = ctx->pool;

    if (concurrency_limit == 0)
        concurrency_limit = mps_thread_get_code_number(ctx);

    if (concurrency_limit < pool->concurrency_limit)
    {
        mps_thread* old_first = pool->first_thread;
        mps_thread* thread;
        unsigned int i = 0;

        for (thread = pool->first_thread; i < (pool->concurrency_limit - concurrency_limit); thread = thread->next_thread, i++)
            ;

        pool->first_thread = thread;  // new head pointer to the threads that will be allowed to remain
        pool->n = concurrency_limit;

        i = 0;
        //  kill the threads that were disallowed
        for (thread = old_first; i < (pool->concurrency_limit - concurrency_limit); i++)
        {
            mps_thread* next_thread = thread->next_thread;
            mps_thread_free(ctx, thread);
            thread = next_thread;
        }
    }
    else
    {
        unsigned int i = 0;
        for (i = 0; i < concurrency_limit - pool->concurrency_limit; i++)
            mps_thread_pool_insert_new_thread(ctx, ctx->pool);
    }

    pool->concurrency_limit = concurrency_limit;
}

// a thread is picking up and starting on the "work" task
void
    mps_thread_pool_assign(mps_context* ctx, mps_thread_pool* pool, mps_thread_work work,
        void* args)
{
    if (!pool)
        pool = ctx->pool;

    if (pool->n == 1 && !pool->strict_async)
    {
        (*work)(args);
        return;
    }
    /* Insert the job in the queue */
    mps_mutex_lock(pool->queue_changed_mutex);

    mps_new_obj(mps_thread_pool_queue_item, item, sizeof(mps_thread_pool_queue_item));

    item->work = work;
    item->work_args = args;

    if (pool->queue->first_item == NULL)
    {
        pool->queue->first_item = pool->queue->last_item = item;
        item->next_item = NULL;
    }
    else
    {
        pool->queue->last_item->next_item = item;
        pool->queue->last_item = item;
        item->next_item = NULL;
    }

    // mps_thread_mainloop wants to know that queue has been altered
    mps_cond_signal(pool->queue_changed_cond, pool->wait_queue_changed_mutex);
    mps_mutex_unlock(pool->queue_changed_mutex);
}

/**
* @brief Wait for an entire thread pool to complete all its jobs
* OR to be told to die.
* There are multiple pools where each pool can have multiple jobs, and it may not be the ctx->pool being referenced
*/
void
    mps_thread_pool_wait(mps_context* ctx, mps_thread_pool* pool)
{
    mps_mutex_lock(pool->work_completed_mutex);

    while (true)
    {
        if (pool->busy_counter == 0 && pool->queue->first_item == NULL)  // no in-process and no new unassigned jobs
        {
            // die
            mps_mutex_unlock(pool->work_completed_mutex);
            return;
        }
        else
        {
            // this thread is busy or the queue is empty.
            // wait until the pool is empty and there is no new work in the queue
            mps_cond_wait(pool->work_completed_cond, pool->wait_work_completed_mutex, pool->work_completed_mutex, pool->busy_counter == 0 && pool->queue->first_item == NULL);
            // the return could be "spurious" so repeat the loop
        }
    }
}

/**
* @brief Allocate a new <code>mps_thread</code> and start its mainloop.
*/
mps_thread*
    mps_thread_new(mps_context* ctx, mps_thread_pool* pool)
{
    if (!pool)
        pool = ctx->pool;
    mps_new_obj(mps_thread, thread, sizeof(mps_thread));

    /* Set the initial values in the thread */
    thread->data = NULL;
    mps_mutex_init(thread->busy_mutex);
    mps_cond_init(thread->start_condition);
    thread->work = NULL;
    thread->args = NULL;
    thread->alive = true;
    thread->pool = pool;
    thread->busy = false;

    /* Create a thread object and start the thread mainloop */
    mps_start_new_thread(thread->sys_thread, mps_thread_mainloop, thread);

    return thread;
}

/**
* @brief Free a thread asking it to stop.
* Either the pool was disallowed due to concurreny limitations,
* or the pool itself is dying (it has completed mps_thread_pool_wait which let all its threads die)
*/

void
    mps_thread_free(mps_context* ctx, mps_thread* thread)
{
    /* Wait for the thread to finish its work, if it is doing something */
    /* mps_mutex_lock (&thread->busy_mutex); */
    /* mps_mutex_unlock (&thread->busy_mutex); */
    mps_mutex_lock(thread->pool->queue_changed_mutex);

    thread->alive = false;

    // the target thread is no longer 'alive' but is waiting for a new task.  
    // (It should not actually be working or else we wouldn't be telling that thread to free its resources).
    // Make it stop waiting and die.
    // "broadcast' may awaken some different threads that are just waiting for a job but they will simply re-wait themselves

    // mps_thread_mainloop wants to know that the thread is no longer alive
    mps_cond_broadcast(thread->pool->queue_changed_cond, thread->pool->wait_queue_changed_mutex);
    mps_mutex_unlock(thread->pool->queue_changed_mutex);

    mps_thread_join(thread->sys_thread);  // the alive flag is now false so the thread will soon check it and exit

    mps_mutex_destroy(thread->busy_mutex);
    mps_cond_destroy(thread->start_condition);

    // then after the thread has exited
    if (thread->sys_thread)
    {
        mps_del_obj(thread->sys_thread);
    }
    mps_del_obj(thread);
}

/**
* @brief Create a new thread and add it to the specified thread pool.
*/
void
    mps_thread_pool_insert_new_thread(mps_context* ctx, mps_thread_pool* pool)
{
    if (!pool)
        pool = ctx->pool;

    mps_thread* thread = mps_thread_new(ctx, pool);

    thread->next_thread = pool->first_thread;
    pool->first_thread = thread;
    pool->n++;
}

/**
* @brief Obtain a pointer to the default shared thread pool on
* this system.
*/
mps_thread_pool*
    mps_thread_pool_get_system_pool(mps_context* ctx)
{
    mps_mutex_lock(system_thread_pool_lock);

    if (system_thread_pool == NULL)
        system_thread_pool = mps_thread_pool_new(ctx, 0);

    mps_mutex_unlock(system_thread_pool_lock);

    return system_thread_pool;
}

/**
* @brief Set the value of the internal field strict_async of the pool.
*
* This is used to control the optimizations performed in the case where
* the number of threads is 1. In case strict_async is set to true (that is
* currently the default value) every call to mps_thread_pool_assign() will not return
* immediately but instead complete the job.
*/
void
    mps_thread_pool_set_strict_async(mps_thread_pool* pool, mps_boolean strict_async)
{
    pool->strict_async = strict_async;
}

/**
* @brief Allocate a new thread pool and return a pointer to it,
* with a number of threads suitable for this system.
*/
mps_thread_pool*
    mps_thread_pool_new(mps_context* ctx, int n_threads)
{
    int i, threads;

    mps_new_obj(mps_thread_pool, pool, sizeof(mps_thread_pool));

    if (n_threads != 0)
    {
        threads = n_threads;
    }
    else
    {
        threads = mps_thread_get_code_number(ctx);
    }

    pool->n = 0;
    pool->first_thread = NULL;

    mps_new_defined_obj(mps_thread_pool_queue, pool->queue, sizeof(mps_thread_pool_queue));
    pool->queue->first_item = pool->queue->last_item = NULL;

    mps_mutex_init(pool->queue_changed_mutex);
    mps_cond_init(pool->queue_changed_cond);

    mps_mutex_init(pool->work_completed_mutex);
    mps_cond_init(pool->work_completed_cond);

    pool->busy_counter = 0;
    pool->strict_async = false;

    for (i = 0; i < threads; i++)
        mps_thread_pool_insert_new_thread(ctx, pool);

    pool->concurrency_limit = threads;

    mps_thread_pool_wait(ctx, pool);

    return pool;
}

/**
* @brief Free a thread pool and all its threads, waiting for them
* to terminate.
*/
void
    mps_thread_pool_free(mps_context* ctx, mps_thread_pool* pool)
{
    if (!pool)
        pool = ctx->pool;

    mps_thread* thread = pool->first_thread;
    mps_thread* next_thread;

    while (thread)
    {
        next_thread = thread->next_thread;
        mps_thread_free(ctx, thread);
        thread = next_thread;
    }

    mps_mutex_destroy(pool->queue_changed_mutex);
    mps_cond_destroy(pool->queue_changed_cond);

    mps_mutex_destroy(pool->work_completed_mutex);
    mps_cond_destroy(pool->work_completed_cond);

    mps_del_obj(pool->queue);

    mps_del_obj(pool);
}

int mps_thread_get_id(mps_context* ctx, mps_thread_pool* pool)
{
    mps_thread_id_t self = mps_thread_self();
    int i = 0;

    mps_thread* thread = pool->first_thread;

    while (thread)
    {
        if (mps_thread_equal(*thread->sys_thread, self))
        {
            return i;
        }
        i++;
        thread = thread->next_thread;
    }
    return -1;
}
