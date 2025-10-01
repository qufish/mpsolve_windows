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
#include <math.h>
#include <string.h>

/* Routine called to perform the floating point newton iterations. */
static void*
    __mps_secular_ga_fiterate_worker(void* data_ptr)
{
    mps_thread_worker_data* data = (mps_thread_worker_data*)data_ptr; // an upstream method will free this
    mps_context* ctx = data->ctx;
    int i;
    cplx_t corr, abcorr;
    double modcorr;
    mps_thread_job job;

    while (true && !ctx->exit_required)
    {
        job = mps_thread_job_queue_next(ctx, data->queue);
        i = job.i;

        if (job.iter == MPS_THREAD_JOB_EXCEP || *data->nzeros >= ctx->n)
            goto cleanup;

        mps_mutex_lock(data->roots_mutex[i]);

        if (job.iter == MPS_THREAD_JOB_EXCEP || *data->nzeros >= ctx->n)
        {
            mps_mutex_unlock(data->roots_mutex[i]);
            goto cleanup;
        }

        if (ctx->approx_root[i]->again && !ctx->approx_root[i]->approximated)
        {
            /* Increment the number of performed iterations */
#if defined(__GCC__)
            __sync_add_and_fetch(data->it, 1);
#else
            mps_mutex_lock(*data->gs_mutex);
            (*data->it)++;
            mps_mutex_unlock(*data->gs_mutex);
#endif
            cdpe_set_x(ctx->approx_root[i]->dvalue, ctx->approx_root[i]->fvalue);

            mps_secular_fnewton(ctx, MPS_POLYNOMIAL(ctx->secular_equation), ctx->approx_root[i], corr);

            if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_NOT_FLOAT)
            {
                *data->excep = true;
                mps_mutex_unlock(data->roots_mutex[i]);
                break;
            }

            /* Apply Aberth correction */
            mps_faberth_wl(ctx, i, abcorr, data->aberth_mutex);

            if (isnan(cplx_Re(abcorr)) || isnan(cplx_Im(abcorr)))
            {
                ctx->approx_root[i]->again = false;
                mps_mutex_unlock(data->roots_mutex[i]);
                continue;
            }

            cplx_mul_eq(abcorr, corr);
            cplx_sub(abcorr, cplx_one, abcorr);
            cplx_div(abcorr, corr, abcorr);

            if (cplx_check_fpe(abcorr))
            {
                ctx->approx_root[i]->again = false;
                mps_mutex_unlock(data->roots_mutex[i]);
                continue;
            }

            if (!ctx->approx_root[i]->again || ctx->approx_root[i]->approximated)
            {
                if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
                    MPS_DEBUG(ctx, "Root %d again was set to false on iteration %d by thread %d", i, *data->it, data->thread);

#if defined(__GCC__)
                __sync_add_and_fetch(data->nzeros, 1);
#else
                mps_mutex_lock(*data->gs_mutex);
                (*data->nzeros)++;
                mps_mutex_unlock(*data->gs_mutex);
#endif
            }
            else
            {
                mps_mutex_lock(data->aberth_mutex[i]);
                cplx_sub_eq(ctx->approx_root[i]->fvalue, abcorr);
                mps_mutex_unlock(data->aberth_mutex[i]);

                /* Correct the radius */
                modcorr = cplx_mod(abcorr);
                ctx->approx_root[i]->frad += modcorr;
            }
        }

        mps_mutex_unlock(data->roots_mutex[i]);
    }

cleanup:

    return NULL;
}

/**
* @brief Routine that performs a block of iteration
* in floating point on the secular equation.
*
* @param ctx the pointer to the mps_context struct.
* @param maxit Maximum number of iteration to perform.
* @param just_regenerated true if this is the first iteration after a coefficient
* regeneration. If just_regenerated is true and the iteration packet is completed
* in less than 2 * (n - computed_roots) iterations that best_approx is set to true
* in ctx->secular_equation so a raise in the precision will be triggered.
* @return The number of approximated roots after the iteration.
*/
MPS_PRIVATE int
    mps_secular_ga_fiterate(mps_context* ctx, int maxit, mps_boolean just_regenerated)
{
    int computed_roots = 0;
    int approximated_roots = 0;
    int root_neighborhood_roots = 0;
    int i;
    int nit = 0;
    int it_threshold = 0;
    mps_boolean excep = false;

#ifndef DISABLE_DEBUG
    clock_t* my_clock = mps_start_timer();
#endif

    ctx->operation = MPS_OPERATION_ABERTH_FP_ITERATIONS;

    mps_new_array_obj(mps_mutex_t, aberth_mutex, sizeof(mps_mutex_t), ctx->n);
    mps_new_array_obj(mps_mutex_t, roots_mutex, sizeof(mps_mutex_t), ctx->n);
    mps_initialized_mutex(gs_mutex, PTHREAD_MUTEX_INITIALIZER);

    for (i = 0; i < ctx->n; i++)
    {
        mps_mutex_init(roots_mutex[i]);
        mps_mutex_init(aberth_mutex[i]);
    }

    mps_new_array_obj(mps_thread_worker_data, data_array, sizeof(mps_thread_worker_data), ctx->n_threads);

    MPS_DEBUG_THIS_CALL(ctx);

    /* Mark the approximated roots as ready for output */
    for (i = 0; i < ctx->n; i++)
    {
        /* Set again to false if the root is already approximated. If a root is approximated but
        * it has less digits than the current precision don't stop the iterations on that component. */
        if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_ISOLATED ||
            ctx->approx_root[i]->status == MPS_ROOT_STATUS_APPROXIMATED)
        {
            if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
            {
                MPS_DEBUG_WITH_INFO(ctx, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
            }
            ctx->approx_root[i]->again = false;

            if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_APPROXIMATED)
                ctx->approx_root[i]->approximated = true;
        }

        if (!ctx->approx_root[i]->again || ctx->approx_root[i]->approximated)
            computed_roots++;
    }

    MPS_DEBUG_WITH_INFO(ctx, "%d roots %s already approximated at the start of the packet",
        computed_roots,
        (computed_roots == 1) ? "is" : "are");

    it_threshold = ctx->n - computed_roots;

    mps_thread_job_queue* queue = mps_thread_job_queue_new(ctx);

    for (i = 0; i < ctx->n_threads; i++)
    {
        data_array[i].it = &nit;
        data_array[i].nzeros = &computed_roots;
        data_array[i].ctx = ctx;
        data_array[i].thread = i;
        data_array[i].n_threads = ctx->n_threads;
        data_array[i].aberth_mutex = aberth_mutex;
        data_array[i].roots_mutex = roots_mutex;
        data_array[i].queue = queue;
        data_array[i].gs_mutex = &gs_mutex;
        data_array[i].excep = &excep;

        mps_thread_pool_assign(ctx, ctx->pool, __mps_secular_ga_fiterate_worker,
            &data_array[i]);
    }

    mps_thread_pool_wait(ctx, ctx->pool);

    /* Check if the roots are improvable in floating point */
    MPS_DEBUG_WITH_INFO(ctx, "Performed %d iterations with floating point arithmetic",
        nit);

    if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        mps_dump(ctx);

    /* Check if we need to get higher precision for the roots */
    ctx->best_approx = true;
    for (i = 0; i < ctx->n; i++)
    {
        if (!ctx->approx_root[i]->approximated)
            ctx->best_approx = false;
        if (ctx->approx_root[i]->approximated)
            approximated_roots++;
        if (!ctx->approx_root[i]->again)
            root_neighborhood_roots++;
    }

    if (just_regenerated && (nit <= it_threshold))
        ctx->best_approx = true;

    MPS_DEBUG_WITH_INFO(ctx, "%d roots are approximated with the current precision", approximated_roots);
    MPS_DEBUG_WITH_INFO(ctx, "%d roots are in the root neighborhood", root_neighborhood_roots);
    MPS_DEBUG_WITH_INFO(ctx, "%d roots have reached a stop condition", computed_roots);

    if (excep)
    {
        MPS_DEBUG_WITH_INFO(ctx, "Switching to DPE arithmetic since there are roots not representable in standard floating point");
        for (i = 0; i < ctx->n; i++)
        {
            cdpe_set_x(ctx->approx_root[i]->dvalue, ctx->approx_root[i]->fvalue);
            rdpe_set_d(ctx->approx_root[i]->drad, ctx->approx_root[i]->frad);
            ctx->approx_root[i]->status = MPS_ROOT_STATUS_CLUSTERED;
        }
        ctx->lastphase = dpe_phase;
    }

    /* Compute the inclusion radii with Gerschgorin so we can compute
    * clusterizations for the roots. */
    /* mps_fradii (ctx, fradii); */
    /* mps_fcluster (ctx, fradii, 2.0 * ctx->n);  */
    /* mps_fmodify (ctx, false);  */

    /* These lines are used to debug the again vector, but are not useful
    * at the moment being */
    if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
        __MPS_DEBUG(ctx, "Again vector = ");
        for (i = 0; i < ctx->n; i++)
        {
            fprintf(ctx->logstr, "%d ", ctx->approx_root[i]->again);
        }
        fprintf(ctx->logstr, "\n");
    }

    if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        mps_dump(ctx);

    /* Count time taken  */
#ifndef DISABLE_DEBUG
    ctx->fp_iteration_time += mps_stop_timer(my_clock);
#endif

    for (i = 0; i < ctx->n; i++)
    {
        mps_mutex_destroy(roots_mutex[i]);
        mps_mutex_destroy(aberth_mutex[i]);
    }
    
    mps_mutex_destroy(gs_mutex);
    
    mps_thread_job_queue_free(queue);
    mps_del_array_obj(data_array);  // the sub-workers did not delete their individual pointers

    mps_del_array_obj(roots_mutex);
    mps_del_array_obj(aberth_mutex);

    /* Return the number of approximated roots */
    return computed_roots;
}

/* Routine called to perform the floating point newton iterations with DPE */
static void*
    __mps_secular_ga_diterate_worker(void* data_ptr)
{
    mps_thread_worker_data* data = (mps_thread_worker_data*)data_ptr; // an upstream method will free this
    mps_context* ctx = data->ctx;
    int i;
    cdpe_t corr, abcorr, droot;
    rdpe_t modcorr;
    mps_thread_job job;

    while (true && !ctx->exit_required)
    {
        job = mps_thread_job_queue_next(ctx, data->queue);
        i = job.i;

        if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
            return NULL;
        }

        mps_mutex_lock(data->roots_mutex[i]);

        if (ctx->approx_root[i]->again && !ctx->approx_root[i]->approximated)
        {
            /* Lock this roots to make sure that we are the only one working on it */
            cdpe_set(droot, ctx->approx_root[i]->dvalue);

            (*data->it)++;

            mps_secular_dnewton(ctx, MPS_POLYNOMIAL(ctx->secular_equation), ctx->approx_root[i], corr);

            /* Apply Aberth correction */
            mps_daberth_wl(ctx, i, abcorr, data->aberth_mutex);
            cdpe_mul_eq(abcorr, corr);
            cdpe_sub(abcorr, cdpe_one, abcorr);
            cdpe_div(abcorr, corr, abcorr);

            cdpe_sub_eq(droot, abcorr);

            /* Correct the radius */
            if (ctx->approx_root[i]->again)
            {
                cdpe_mod(modcorr, abcorr);
                rdpe_add_eq(ctx->approx_root[i]->drad, modcorr);
            }

            if (!ctx->approx_root[i]->again || ctx->approx_root[i]->approximated)
            {
                if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
                    MPS_DEBUG(ctx, "Root %d again was set to false on iteration %d by thread %d", i, *data->it, data->thread);
                (*data->nzeros)++;
            }
            else
                cdpe_set(ctx->approx_root[i]->dvalue, droot);
        }

        mps_mutex_unlock(data->roots_mutex[i]);
    }

    return NULL;
}

/**
* @brief Routine that performs a block of iteration
* in floating point on the secular equation using
* CDPE
*
* @param ctx the pointer to the mps_context struct.
* @param maxit Maximum number of iteration to perform.
* @return The number of approximated roots after the iteration.
* @param just_regenerated true if this is the first iteration after a coefficient
* regeneration. If just_regenerated is true and the iteration packet is completed
* in less than 2 * (n - computed_roots) iterations that best_approx is set to true
* in ctx->secular_equation so a raise in the precision will be triggered.
*/
MPS_PRIVATE int
    mps_secular_ga_diterate(mps_context* ctx, int maxit, mps_boolean just_regenerated)
{
    int computed_roots = 0;
    int root_neighborhood_roots = 0;
    int approximated_roots = 0;
    int i;
    int nit = 0;
    int it_threshold = 0;

    ctx->operation = MPS_OPERATION_ABERTH_DPE_ITERATIONS;

#ifndef DISABLE_DEBUG
    clock_t* my_clock = mps_start_timer();
#endif

    mps_new_array_obj(mps_mutex_t, aberth_mutex, sizeof(mps_mutex_t), ctx->n);
    mps_new_array_obj(mps_mutex_t, roots_mutex, sizeof(mps_mutex_t), ctx->n);

    for (i = 0; i < ctx->n; i++)
    {
        mps_mutex_init(roots_mutex[i]);
        mps_mutex_init(aberth_mutex[i]);
    }

    mps_new_array_obj(mps_thread_worker_data, data_array, sizeof(mps_thread_worker_data), ctx->n_threads);

    MPS_DEBUG_THIS_CALL(ctx);

    ctx->best_approx = false;

    /* Mark the approximated roots as ready for output */
    for (i = 0; i < ctx->n; i++)
    {
        /* Set again to false if the root is already approximated. If a root is approximated but
        * it has less digits than the current precision don't stop the iterations on that component. */
        if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_ISOLATED ||
            ctx->approx_root[i]->status == MPS_ROOT_STATUS_APPROXIMATED)
        {
            if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
            {
                MPS_DEBUG_WITH_INFO(ctx, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
            }
            ctx->approx_root[i]->again = false;
            ctx->approx_root[i]->approximated = true;
        }

        if (!ctx->approx_root[i]->again || ctx->approx_root[i]->approximated)
            computed_roots++;
    }

    it_threshold = (ctx->n - computed_roots);

    MPS_DEBUG_WITH_INFO(ctx, "%d roots %s already approximated at the start of the packet",
        computed_roots,
        (computed_roots == 1) ? "is" : "are");

    mps_thread_job_queue* queue = mps_thread_job_queue_new(ctx);

    for (i = 0; i < ctx->n_threads; i++)
    {
        data_array[i].it = &nit;
        data_array[i].nzeros = &computed_roots;
        data_array[i].ctx = ctx;
        data_array[i].thread = i;
        data_array[i].n_threads = ctx->n_threads;
        data_array[i].aberth_mutex = aberth_mutex;
        data_array[i].roots_mutex = roots_mutex;
        data_array[i].queue = queue;

        mps_thread_pool_assign(ctx, ctx->pool, __mps_secular_ga_diterate_worker,
            &data_array[i]);
    }

    mps_thread_pool_wait(ctx, ctx->pool);

    /* Check if the roots are improvable in floating point */
    MPS_DEBUG_WITH_INFO(ctx, "Performed %d iterations with CDPE arithmetic",
        nit);

    if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        mps_dump(ctx);

    /* Check if we need to get higher precision for the roots */
    ctx->best_approx = true;
    for (i = 0; i < ctx->n; i++)
    {
        if (!ctx->approx_root[i]->approximated)
            ctx->best_approx = false;
        if (ctx->approx_root[i]->approximated)
            approximated_roots++;
        if (!ctx->approx_root[i]->again)
            root_neighborhood_roots++;
    }

    if (just_regenerated && (nit <= it_threshold))
        ctx->best_approx = true;

    MPS_DEBUG_WITH_INFO(ctx, "%d roots are approximated with the current precision", approximated_roots);
    MPS_DEBUG_WITH_INFO(ctx, "%d roots are in the root neighborhood", root_neighborhood_roots);
    MPS_DEBUG_WITH_INFO(ctx, "%d roots have reached a stop condition", computed_roots);

    /* Clock the routine */
#ifndef DISABLE_DEBUG
    ctx->dpe_iteration_time += mps_stop_timer(my_clock);
#endif

    for (i = 0; i < ctx->n; i++)
    {
        mps_mutex_destroy(roots_mutex[i]);
        mps_mutex_destroy(aberth_mutex[i]);
    }

    mps_thread_job_queue_free(queue);
    mps_del_array_obj(aberth_mutex);
    mps_del_array_obj(roots_mutex);
 
    mps_del_array_obj(data_array);  // the sub-workers did not delete their individual pointers

    /* Return the number of approximated roots */
    return computed_roots;
}

/* Routine called to perform the floating point newton iterations with MP */
static void*
    __mps_secular_ga_miterate_worker(void* data_ptr)
{
    mps_thread_worker_data* data = (mps_thread_worker_data*)data_ptr; // an upstream method will free this
    mps_context* ctx = data->ctx;
    /* mps_secular_equation * sec = ctx->secular_equation; */
    int i;
    mpc_t corr, abcorr;
    mpc_t mroot;
    rdpe_t modcorr;
    mps_thread_job job;

    mps_cluster* cluster = NULL;

    mpc_init2(corr, ctx->mpwp);
    mpc_init2(abcorr, ctx->mpwp);
    mpc_init2(mroot, ctx->mpwp);

    /* Get a copy of the MP coefficients that is local to this thread */
    while (true && !ctx->exit_required)
    {
        job = mps_thread_job_queue_next(ctx, data->queue);
        i = job.i;

        if (job.iter == MPS_THREAD_JOB_EXCEP || *data->nzeros >= ctx->n)
            goto cleanup;

        mps_mutex_lock(data->roots_mutex[i]);

        if (job.iter == MPS_THREAD_JOB_EXCEP || *data->nzeros >= ctx->n)
        {
            mps_mutex_unlock(data->roots_mutex[i]);
            goto cleanup;
        }

        /* printf ("Thread %d iterating on root %d\n", data->thread, i); */

        cluster = job.cluster_item->cluster;

        if (ctx->approx_root[i]->again && !ctx->approx_root[i]->approximated)
        {
            /* Lock this roots to make sure that we are the only one working on it */
            mps_mutex_lock(data->aberth_mutex[i]);
            mpc_set(mroot, ctx->approx_root[i]->mvalue);
            mps_mutex_unlock(data->aberth_mutex[i]);

            /* Check if, while we were waiting, excep condition has been reached,
            * or all the zeros has been approximated.                         */
            if ((*data->nzeros) >= ctx->n)
            {
                mps_mutex_unlock(data->roots_mutex[i]);
                goto cleanup;
            }

            /* mps_mutex_lock (* data->gs_mutex); */
            (*data->it)++;
            /* mps_mutex_unlock (* data->gs_mutex); */

            mps_secular_mnewton(ctx, MPS_POLYNOMIAL(ctx->secular_equation), ctx->approx_root[i],
                corr, mpc_get_prec(ctx->approx_root[i]->mvalue));

            /* Apply Aberth correction */
            mps_maberth_s_wl(ctx, i, cluster, abcorr, data->aberth_mutex);
            mpc_mul_eq(abcorr, corr);
            mpc_ui_sub(abcorr, 1U, 0U, abcorr);

            if (!mpc_eq_zero(abcorr))
            {
                mpc_div(abcorr, corr, abcorr);

                mps_mutex_lock(data->aberth_mutex[i]);
                mpc_sub_eq(mroot, abcorr);
                mps_mutex_unlock(data->aberth_mutex[i]);
            }
            else
                ctx->approx_root[i]->again = true;


            if (!ctx->approx_root[i]->again || ctx->approx_root[i]->approximated)
            {
                if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
                    MPS_DEBUG(ctx, "Root %d again was set to false on iteration %d by thread %d", i, *data->it, data->thread);

                (*data->nzeros)++;
            }
            else
            {
                mps_mutex_lock(data->aberth_mutex[i]);
                mpc_set(ctx->approx_root[i]->mvalue, mroot);
                mps_mutex_unlock(data->aberth_mutex[i]);

                /* Correct the radius */
                mpc_rmod(modcorr, abcorr);
                rdpe_add_eq(ctx->approx_root[i]->drad, modcorr);

                mpc_rmod(modcorr, mroot);
                rdpe_mul_eq(modcorr, ctx->mp_epsilon);
                rdpe_add_eq(ctx->approx_root[i]->drad, modcorr);
            }
        }

        mps_mutex_unlock(data->roots_mutex[i]);
    }

cleanup:
    mpc_clear(mroot);
    mpc_clear(abcorr);
    mpc_clear(corr);

    return NULL;
}

/**
* @brief Routine that performs a block of iteration
* in floating point on the secular equation using
* CDPE
*
* @param ctx the pointer to the mps_context struct.
* @param maxit Maximum number of iteration to perform.
* @param just_regenerated true if this is the first iteration after a coefficient
* regeneration. If just_regenerated is true and the iteration packet is completed
* in less than 2 * (n - computed_roots) iterations that best_approx is set to true
* in ctx->secular_equation so a raise in the precision will be triggered.
* @return The number of approximated roots after the iteration.
*/
MPS_PRIVATE int
    mps_secular_ga_miterate(mps_context* ctx, int maxit, mps_boolean just_regenerated)
{
    int computed_roots = 0;
    int approximated_roots = 0;
    int root_neighborhood_roots = 0;
    int i;
    int nit = 0;
    int it_threshold = 0;

    ctx->operation = MPS_OPERATION_ABERTH_MP_ITERATIONS;

#ifndef DISABLE_DEBUG
    clock_t* my_clock = mps_start_timer();
#endif

    mps_new_array_obj(mps_mutex_t, aberth_mutex, sizeof(mps_mutex_t), ctx->n);
    mps_new_array_obj(mps_mutex_t, roots_mutex, sizeof(mps_mutex_t), ctx->n);

    static mps_initialized_mutex(gs_mutex, PTHREAD_MUTEX_INITIALIZER);

    for (i = 0; i < ctx->n; i++)
    {
        mps_mutex_init(roots_mutex[i]);
        mps_mutex_init(aberth_mutex[i]);
    }

    mps_new_array_obj(mps_thread_worker_data, data_array, sizeof(mps_thread_worker_data), ctx->n_threads);

    MPS_DEBUG_THIS_CALL(ctx);

    ctx->best_approx = false;

    it_threshold = ctx->n - computed_roots;

    /* Mark the approximated roots as ready for output */
    for (i = 0; i < ctx->n; i++)
    {
        /* Set again to false if the root is already approximated. If a root is approximated but
        * it has less digits than the current precision don't stop the iterations on that component. */
        if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_ISOLATED ||
            ctx->approx_root[i]->status == MPS_ROOT_STATUS_APPROXIMATED)
        {
            if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
            {
                MPS_DEBUG_WITH_INFO(ctx, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
            }
            ctx->approx_root[i]->again = false;
            ctx->approx_root[i]->approximated = true;
        }

        if (!ctx->approx_root[i]->again || ctx->approx_root[i]->approximated)
            computed_roots++;
    }

    mps_thread_job_queue* queue = mps_thread_job_queue_new(ctx);

    for (i = 0; i < ctx->n_threads; i++)
    {
        data_array[i].it = &nit;
        data_array[i].nzeros = &computed_roots;
        data_array[i].ctx = ctx;
        data_array[i].thread = i;
        data_array[i].n_threads = ctx->n_threads;
        data_array[i].aberth_mutex = aberth_mutex;
        data_array[i].roots_mutex = roots_mutex;
        data_array[i].queue = queue;
        data_array[i].gs_mutex = &gs_mutex;

        mps_thread_pool_assign(ctx, ctx->pool, __mps_secular_ga_miterate_worker,
            &data_array[i]);
    }

    mps_thread_pool_wait(ctx, ctx->pool);

    /* Check if the roots are improvable in floating point */
    MPS_DEBUG_WITH_INFO(ctx, "Performed %d iterations with MP arithmetic",
        nit);

    if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        mps_dump(ctx);

    /* Check if we need to get higher precision for the roots */
    ctx->best_approx = true;
    for (i = 0; i < ctx->n; i++)
    {
        if (!ctx->approx_root[i]->approximated)
            ctx->best_approx = false;
        if (ctx->approx_root[i]->approximated)
            approximated_roots++;
        if (!ctx->approx_root[i]->again)
            root_neighborhood_roots++;
    }

    if (just_regenerated && (nit <= it_threshold))
        ctx->best_approx = true;

    MPS_DEBUG_WITH_INFO(ctx, "%d roots are approximated with the current precision", approximated_roots);
    MPS_DEBUG_WITH_INFO(ctx, "%d roots are in the root neighborhood", root_neighborhood_roots);
    MPS_DEBUG_WITH_INFO(ctx, "%d roots have reached a stop condition", computed_roots);

    /* These lines are used to debug the again vector, but are not useful
    * at the moment being */
    if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
        __MPS_DEBUG(ctx, "Again vector = ");
        for (i = 0; i < ctx->n; i++)
        {
            fprintf(ctx->logstr, "%d ", ctx->approx_root[i]->again);
        }
        fprintf(ctx->logstr, "\n");
    }

    /* Clock the routine */
#ifndef DISABLE_DEBUG
    ctx->mp_iteration_time += mps_stop_timer(my_clock);
#endif

    for (i = 0; i < ctx->n; i++)
    {
        mps_mutex_destroy(roots_mutex[i]);
        mps_mutex_destroy(aberth_mutex[i]);
    }

    mps_thread_job_queue_free(queue);
    mps_del_array_obj(aberth_mutex);
    mps_del_array_obj(roots_mutex);

    mps_del_array_obj(data_array);  // the sub-workers did not delete their individual pointers

    /* Return the number of approximated roots */
    return computed_roots;
}
