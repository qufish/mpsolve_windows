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

#include <mps/mps.h>

/**
* @brief Worker for the fpolzer routine.
*/
void*
    mps_thread_fpolzer_worker(void* data_ptr)
{
    mps_thread_worker_data* data = (mps_thread_worker_data*)data_ptr; // an upstream method will free this
    mps_context* ctx = data->ctx;
    mps_polynomial* p = ctx->active_poly;
    int i, iter;
    cplx_t corr, abcorr, froot;
    double rad1_double, modcorr;
    mps_thread_job job;

    while (!(*data->excep) && (*data->nzeros) < data->required_zeros)
    {
        job = mps_thread_job_queue_next(ctx, data->queue);
        i = job.i;
        iter = job.iter;

        /* Check if we got over the maximum number of iterations */
        if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
            (*data->excep) = true;
            return 0;
        }

        /* Lock this roots to make sure that we are the only one working on it */
        mps_mutex_lock(data->roots_mutex[i]);

        if (ctx->approx_root[i]->again)
        {
            /* Check if, while we were waiting, excep condition has been reached */
            if (*data->excep || !ctx->approx_root[i]->again || (*data->nzeros > data->required_zeros))
            {
                mps_mutex_unlock(data->roots_mutex[i]);
                return 0;
            }

            (*data->it)++;
            rad1_double = ctx->approx_root[i]->frad;

            /* Make a local copy of the root */
            mps_mutex_lock(data->aberth_mutex[i]);
            cplx_set(froot, ctx->approx_root[i]->fvalue);
            mps_mutex_unlock(data->aberth_mutex[i]);

            mps_polynomial_fnewton(ctx, p, ctx->approx_root[i], corr);

            if (cplx_check_fpe(corr))
            {
                /* If we get a floating point exception we need to switch to DPE
                * arithmetic. */
                ctx->approx_root[i]->frad = rad1_double;
                ctx->skip_float = true;
                ctx->approx_root[i]->again = false;
            }

            if (iter == 0 && !ctx->approx_root[i]->again && ctx->approx_root[i]->frad > rad1_double && rad1_double != 0)
                ctx->approx_root[i]->frad = rad1_double;
            /***************************************
                The above condition is needed to cope with the case
                where at the first iteration the starting point
                is already in the root neighbourhood and the actually
                computed radius is too big since the value of the first
                derivative is too small.
                In this case the previous radius bound, obtained by
                means of Rouche' is more reliable and strict
            **************************************/

            if (ctx->approx_root[i]->again
                /* the correction is performed only if iter!=1 or rad(i)!=rad1_double */
                || iter != 0 || ctx->approx_root[i]->frad != rad1_double)
            {
                mps_faberth(ctx, ctx->approx_root[i], abcorr);

                cplx_mul_eq(abcorr, corr);
                cplx_sub(abcorr, cplx_one, abcorr);

                if (cplx_eq_zero(abcorr))
                {
                    MPS_DEBUG(ctx, "Aberth correction is zero");
                    cplx_set_d(abcorr, DBL_EPSILON, 0);
                }

                cplx_div(abcorr, corr, abcorr);
                cplx_sub_eq(froot, abcorr);
                modcorr = cplx_mod(abcorr);
                ctx->approx_root[i]->frad += modcorr;

                mps_mutex_lock(data->aberth_mutex[i]);
                cplx_set(ctx->approx_root[i]->fvalue, froot);
                mps_mutex_unlock(data->aberth_mutex[i]);
            }

            /* check for new approximated roots */
            if (!ctx->approx_root[i]->again)
            {
                (*data->nzeros)++;
                if (*data->nzeros >= data->required_zeros)
                {
                    mps_mutex_unlock(data->roots_mutex[i]);
                    return 0;
                }
            }
        }

        mps_mutex_unlock(data->roots_mutex[i]);
    }

    return NULL;
}

/**
* @brief Drop-in replacement for the stock fpolzer routine.
* This version adds multithread support.
*/
MPS_PRIVATE void
    mps_thread_fpolzer(mps_context* ctx, int* it, mps_boolean* excep, int required_zeros)
{
    int i, nzeros = 0, n_threads = ctx->n_threads;

    mps_new_array_obj(mps_mutex_t, aberth_mutex, sizeof(mps_mutex_t), ctx->n);
    mps_new_array_obj(mps_mutex_t, roots_mutex, sizeof(mps_mutex_t), ctx->n);

    for (i = 0; i < ctx->n; i++)
    {
        mps_mutex_init(roots_mutex[i]);
        mps_mutex_init(aberth_mutex[i]);
    }

    /* Create a new job queue */
    mps_thread_job_queue* queue = mps_thread_job_queue_new(ctx);

    *it = 0;
    *excep = false;

    /* count the number of approximations in the root neighbourhood */
    for (i = 0; i < ctx->n; i++)
        if (!ctx->approx_root[i]->again)
            nzeros++;
    if (nzeros == ctx->n)
    {
        mps_del_array_obj(roots_mutex);
        mps_del_array_obj(aberth_mutex);
        mps_thread_job_queue_free(queue);
        return;
    }

    mps_new_array_obj(mps_thread_worker_data, data_array, sizeof(mps_thread_worker_data), n_threads);

    for (i = 0; i < n_threads; i++)
    {
        data_array[i].it = it;
        data_array[i].nzeros = &nzeros;
        data_array[i].ctx = ctx;
        data_array[i].excep = excep;
        data_array[i].thread = i;
        data_array[i].n_threads = n_threads;
        data_array[i].aberth_mutex = aberth_mutex;
        data_array[i].roots_mutex = roots_mutex;
        data_array[i].queue = queue;
        data_array[i].required_zeros = required_zeros;
        /* mps_thread_create (&threads[i], NULL, &mps_thread_fpolzer_worker, */
        /* data_array[i]); */
        mps_thread_pool_assign(ctx, ctx->pool, mps_thread_fpolzer_worker, &data_array[i]);
    }

    mps_thread_pool_wait(ctx, ctx->pool);

    for (i = 0; i < ctx->n; i++)
    {
        mps_mutex_destroy(roots_mutex[i]);
        mps_mutex_destroy(aberth_mutex[i]);
    }
    mps_del_array_obj(data_array);
    mps_del_array_obj(roots_mutex);
    mps_del_array_obj(aberth_mutex);
    mps_thread_job_queue_free(queue);
}

/**
* @brief Multithread worker for mps_thread_dpolzer ()
*/
static void*
    mps_thread_dpolzer_worker(void* data_ptr)
{
    int iter, i;
    rdpe_t rad1_rdpe_t, rtmp;
    cdpe_t corr, abcorr;

    /* Parse input data */
    mps_thread_worker_data* data = (mps_thread_worker_data*)data_ptr; // an upstream method will free this
    mps_context* ctx = data->ctx;
    mps_polynomial* p = ctx->active_poly;
    mps_thread_job job;

    while (!(*data->excep) && (*data->nzeros < data->required_zeros))
    {
        job = mps_thread_job_queue_next(ctx, data->queue);
        i = job.i;
        iter = job.iter;

        /* Check if we got over the maximum number of iterations */
        if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
            (*data->excep) = true;
            return 0;
        }

        /* Make sure that we are the only one iterating on this root */
        if (ctx->pool->n > 1)
            mps_mutex_lock(data->roots_mutex[i]);

        if (ctx->approx_root[i]->again)
        {
            /* Check if, while we were waiting, excep condition has been reached */
            if (*data->excep || !ctx->approx_root[i]->again || (*data->nzeros > data->required_zeros))
            {
	        if (ctx->pool->n > 1)
                    mps_mutex_unlock(data->roots_mutex[i]);
                return 0;
            }

            (*data->it)++;
            rdpe_set(rad1_rdpe_t, ctx->approx_root[i]->drad);

            mps_polynomial_dnewton(ctx, p, ctx->approx_root[i], corr);
            if (iter == 0 && !ctx->approx_root[i]->again && rdpe_gt(ctx->approx_root[i]->drad, rad1_rdpe_t)
                && rdpe_ne(rad1_rdpe_t, rdpe_zero))
                rdpe_set(ctx->approx_root[i]->drad, rad1_rdpe_t);
            /************************************************
                The above condition is needed to manage with the case where
                at the first iteration the starting point is already in the
                root neighbourhood and the actually computed radius is too
                big since the value of the first derivative is too small.
                In this case the previous radius bound, obtained by means of
                Rouche' is more reliable and strict
            **********************************************/

            if (ctx->approx_root[i]->again
                /* the correction is performed only if iter!=1 or rad(i)!=rad1_rdpe_t */
                || iter != 0
                || rdpe_ne(ctx->approx_root[i]->drad, rad1_rdpe_t))
            {
                mps_daberth(ctx, ctx->approx_root[i], abcorr);
                cdpe_mul_eq(abcorr, corr);
                cdpe_sub(abcorr, cdpe_one, abcorr);
                if (cdpe_eq_zero(abcorr))
                {
                    MPS_DEBUG(ctx, "Aberth correction is zero.");
                    ctx->lastphase = dpe_phase;
                    cdpe_set_d(abcorr, DBL_EPSILON, 0);
                }

                cdpe_div(abcorr, corr, abcorr);
                cdpe_sub_eq(ctx->approx_root[i]->dvalue, abcorr);
                cdpe_mod(rtmp, abcorr);
                rdpe_add_eq(ctx->approx_root[i]->drad, rtmp);
            }

            /* check for new approximated roots */
            if (!ctx->approx_root[i]->again)
            {
                (*data->nzeros)++;
                if ((*data->nzeros) >= data->required_zeros)
                {
		    if (ctx->pool->n > 1)
                        mps_mutex_unlock(data->roots_mutex[i]);
                    return 0;
                }
            }
        }

        if (ctx->pool->n > 1)
            mps_mutex_unlock(data->roots_mutex[i]);
    }

    return NULL;
}

/**
* @brief Multithread version of mps_dpolzer ().
*/
MPS_PRIVATE void
    mps_thread_dpolzer(mps_context* ctx, int* it, mps_boolean* excep, int required_zeros)
{
    int i, nzeros = 0;

    /* initialize the iteration counter */
    *it = 0;
    *excep = false;

    /* count the number of approximations in the root neighbourhood */
    for (i = 0; i < ctx->n; i++)
        if (!ctx->approx_root[i]->again)
            nzeros++;
    if (nzeros == ctx->n)
        return;

    /* Prepare queue */
    mps_thread_job_queue* queue = mps_thread_job_queue_new(ctx);

    /* Allocate space for thread data */
    mps_new_array_obj(mps_thread_worker_data, data_array, sizeof(mps_thread_worker_data), ctx->n_threads);

    /* Allocate mutexes and init them */
    mps_new_array_obj(mps_mutex_t, aberth_mutex, sizeof(mps_mutex_t), ctx->n);
    mps_new_array_obj(mps_mutex_t, roots_mutex, sizeof(mps_mutex_t), ctx->n);
    for (i = 0; i < ctx->n; i++)
    {
        if (ctx->pool->n > 1)
	  {
            mps_mutex_init(aberth_mutex[i]);
            mps_mutex_init(roots_mutex[i]);
          }
    }

    /* Start spawning thread */
    for (i = 0; i < ctx->n_threads; i++)
    {
        data_array[i].aberth_mutex = aberth_mutex;
        data_array[i].excep = excep;
        data_array[i].it = it;
        data_array[i].n_threads = ctx->n_threads;
        data_array[i].nzeros = &nzeros;
        data_array[i].queue = queue;
        data_array[i].roots_mutex = roots_mutex;
        data_array[i].ctx = ctx;
        data_array[i].thread = i;
        data_array[i].required_zeros = required_zeros;
        mps_thread_pool_assign(ctx, ctx->pool, mps_thread_dpolzer_worker, &data_array[i]);
    }

    /* Wait for the thread to complete */
    mps_thread_pool_wait(ctx, ctx->pool);

    for (i = 0; i < ctx->n; i++)
    {
        if (ctx->pool->n > 1)
        {
            mps_mutex_destroy(roots_mutex[i]);
            mps_mutex_destroy(aberth_mutex[i]);
        }
    }
    mps_del_array_obj(roots_mutex);
    mps_del_array_obj(aberth_mutex);
    mps_del_array_obj(data_array);  // the sub-workers did not delete their individual pointers
    mps_thread_job_queue_free(queue);
}

/**
* @brief Worker for the mpolzer routine.
*/
static void*
    mps_thread_mpolzer_worker(void* data_ptr)
{
    mps_thread_worker_data* data = (mps_thread_worker_data*)data_ptr; // an upstream method will free this
    mps_context* ctx = data->ctx;
    mps_polynomial* p = ctx->active_poly;
    mps_thread_job job;
    int iter, l;
    mpc_t corr, abcorr, mroot, diff;
    rdpe_t eps, rad1_rdpe_t, rtmp;
    cdpe_t ctmp;

    mpc_init2(abcorr, ctx->mpwp);
    mpc_init2(corr, ctx->mpwp);
    mpc_init2(mroot, ctx->mpwp);
    mpc_init2(diff, ctx->mpwp);

    rdpe_mul_d(eps, ctx->mp_epsilon, (double)4 * ctx->n);

    /* Continue to iterate while exception condition has not
    * been reached and there more roots to approximate   */
    while ((*data->nzeros) < data->required_zeros)
    {
        /* Get next job for this thread */
        job = mps_thread_job_queue_next(ctx, data->queue);

        /* Set variables to be used in the rest of the code */
        iter = job.iter;

        /* Check if we exceeded the maximum number of iterations */
        if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
            (*data->excep) = true;
            goto endfun;
        }

        l = job.i;

        /* Lock roots_mutex to assure that we are the only thread
        * working on this root. Parallel computation on the same
        * root is not useful, since we would be performing the
        * same computations.                                  */
        if (ctx->pool->n > 1)
            mps_mutex_lock(data->roots_mutex[l]);

        /* MPS_DEBUG (ctx, "Iterating on root %d, iter %d", l, job.iter); */

        if (ctx->approx_root[l]->again)
        {
            /* Check if, while we were waiting, excep condition has been reached,
            * or all the zeros has been approximated.                         */
            if (*data->excep || (*data->nzeros) >= data->required_zeros)
            {
	        if (ctx->pool->n > 1)
                    mps_mutex_unlock(data->roots_mutex[l]);
                goto endfun;
            }

            /* Increment total iteration counter */
            (*data->it)++;

            /* Copy locally the root to work on */
	    if (ctx->pool->n > 1)
              mps_mutex_lock(data->aberth_mutex[l]);
            mpc_set(mroot, ctx->approx_root[l]->mvalue);
	    if (ctx->pool->n > 1)
                mps_mutex_unlock(data->aberth_mutex[l]);

            /* sparse/dense polynomial */
            rdpe_set(rad1_rdpe_t, ctx->approx_root[l]->drad);

            mps_polynomial_mnewton(ctx, p, ctx->approx_root[l], corr, ctx->mpwp);

            if (iter == 0 && !ctx->approx_root[l]->again && rdpe_gt(ctx->approx_root[l]->drad, rad1_rdpe_t)
                && rdpe_ne(rad1_rdpe_t, rdpe_zero))
                rdpe_set(ctx->approx_root[l]->drad, rad1_rdpe_t);

            /************************************************
                The above condition is needed to cope with the case
                where at the first iteration the starting point is
                already in the root neighbourhood and the actually
                computed radius is too big since the value of the
                first derivative is too small.
                In this case the previous radius bound, obtained by
                means of Rouche' is more reliable and strict
            ***********************************************/

            if (ctx->approx_root[l]->again
                /* the correction is performed only if iter!=1 or rad[l]!=rad1_rdpe_t */
                || iter != 0
                || rdpe_ne(ctx->approx_root[l]->drad, rad1_rdpe_t))
            {
                /* Global lock to aberth step to reach a real Gauss-Seidel iteration */
	        if (ctx->pool->n > 1)
                    mps_mutex_lock(*data->global_aberth_mutex);

                /* Compute Aberth correction with locks so we can lock the
                * roots while reading them.                          */
                mps_maberth_s_wl(ctx, l, job.cluster_item->cluster, abcorr,
                    data->aberth_mutex);

                /* Apply aberth correction that has been computed */
                mpc_mul_eq(abcorr, corr);
                mpc_neg_eq(abcorr);
                mpc_add_eq_ui(abcorr, 1, 0);
                mpc_div(abcorr, corr, abcorr);
                mpc_sub_eq(mroot, abcorr);
                mpc_get_cdpe(ctmp, abcorr);
                cdpe_mod(rtmp, ctmp);
                rdpe_add_eq(ctx->approx_root[l]->drad, rtmp);

                /* Lock aberth_mutex and copy the computed root back
                * to its place                                   */
	        if (ctx->pool->n > 1)
                    mps_mutex_lock(data->aberth_mutex[l]);
                mpc_set(ctx->approx_root[l]->mvalue, mroot);
	        if (ctx->pool->n > 1)
                    mps_mutex_unlock(data->aberth_mutex[l]);

                /* Go with others aberth iterations */
	        if (ctx->pool->n > 1)
                    mps_mutex_unlock(*data->global_aberth_mutex);
            }

            /* check for new approximated roots */
            if (!ctx->approx_root[l]->again)
            {
                (*data->nzeros)++;
                if ((*data->nzeros) >= data->required_zeros)
                {
	            if (ctx->pool->n > 1)
                        mps_mutex_unlock(data->roots_mutex[l]);
                    goto endfun;
                }
            }
        }

	if (ctx->pool->n > 1)
            mps_mutex_unlock(data->roots_mutex[l]);

        /* MPS_DEBUG_MPC (ctx, 15, ctx->approx_root[l]->mvalue, "ctx->mroot[%d]", l); */
        /* MPS_DEBUG_RDPE (ctx, ctx->approx_root[l]->drad, "ctx->drad[%d]", l); */

        if ((*data->nzeros) == ctx->n)
        {
            goto endfun;
        }
    }

endfun:                        /* free local MP variables */
    mpc_clear(corr);
    mpc_clear(abcorr);
    mpc_clear(mroot);
    mpc_clear(diff);

    return NULL;
}

/**
* @brief Drop-in threaded replacement for the stock mpolzer.
*/
MPS_PRIVATE void
    mps_thread_mpolzer(mps_context* ctx, int* it, mps_boolean* excep, int required_zeros)
{
    int i, nzeros = 0, n_threads = ctx->n_threads;

    *it = 0;
    *excep = false;

    /* Check if we have already approxmiated roots */
    for (i = 0; i < ctx->n; i++)
        if (!ctx->approx_root[i]->again)
            nzeros++;
    if (nzeros == ctx->n)
    {
        return;
    }

    /* Lower the number of threads if there are a lot of approximated roots */
    if (ctx->n_threads > (ctx->n - nzeros))
        n_threads = ctx->n - nzeros;
    else
        n_threads = ctx->n_threads;

    MPS_DEBUG_WITH_INFO(ctx, "Spawning %d worker", n_threads);

    /* Allocate and the init mutexes needed by the routine */
    mps_new_array_obj(mps_mutex_t, roots_mutex, sizeof(mps_mutex_t), ctx->n);
    mps_new_array_obj(mps_mutex_t, aberth_mutex, sizeof(mps_mutex_t), ctx->n);
    mps_initialized_mutex(global_aberth_mutex, PTHREAD_MUTEX_INITIALIZER);

    for (i = 0; i < ctx->n; i++)
    {
        mps_mutex_init(aberth_mutex[i]);
        mps_mutex_init(roots_mutex[i]);
    }

    /* Create a new work queue */
    mps_thread_job_queue* queue = mps_thread_job_queue_new(ctx);

    mps_new_array_obj(mps_thread_worker_data, data_array, sizeof(mps_thread_worker_data), n_threads);

    /* Set data to be passed to every thread and actually spawn the threads. */
    for (i = 0; i < n_threads; i++)
    {
        data_array[i].it = it;
        data_array[i].nzeros = &nzeros;
        data_array[i].ctx = ctx;
        data_array[i].excep = excep;
        data_array[i].thread = i;
        data_array[i].n_threads = n_threads;
        data_array[i].aberth_mutex = aberth_mutex;
        data_array[i].global_aberth_mutex = &global_aberth_mutex;
        data_array[i].queue = queue;
        data_array[i].roots_mutex = roots_mutex;
        data_array[i].required_zeros = required_zeros;
        mps_thread_pool_assign(ctx, ctx->pool, mps_thread_mpolzer_worker, &data_array[i]);
    }

    /* Wait for the threads to complete */
    mps_thread_pool_wait(ctx, ctx->pool);

    /* Free data and exit */
    mps_del_array_obj(data_array);  // the sub-workers did not delete their individual pointers
    for (i = 0; i < ctx->n; i++)
    {
        mps_mutex_destroy(roots_mutex[i]);
        mps_mutex_destroy(aberth_mutex[i]);
    }
    mps_del_array_obj(roots_mutex);
    mps_del_array_obj(aberth_mutex);
    mps_thread_job_queue_free(queue);
}
