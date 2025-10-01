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
#include <mps/link.h>
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef MPS_CATCH_FPE
#include <fenv.h>
int feenableexcept(int excepts);
#endif

#if HAVE_CONFIG_H
# include <config.h>
#endif

/**
* @brief Perform some preliminary checks and setup before starting the real
* MPSolve loop.
*/
static void
    mps_preliminary_setup(mps_context* ctx)
{
    /* Make sure that non thread safe polynomial implementations are handled
    * in a safe way. */
    if (!ctx->active_poly->thread_safe)
    {
        mps_thread_pool_set_concurrency_limit(ctx, NULL, 1);
    }
}

/**
* @brief Call the real polynomial (or secular equation, or whatever) solver
* and do the computation.
*
* The algorithm used must be selected before this call with <code>mps_select_algorithm</code>
* and the data (the coefficients, or whatever the algorithm may require) should be provided
* after that.
*
* Roots can then be obtained with the functions <code>mps_context_get_roots_*</code>
*
*/
void
    mps_mpsolve(mps_context* ctx)
{
#ifdef MPS_CATCH_FPE
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

    if (mps_context_has_errors(ctx))
        return;

    mps_preliminary_setup(ctx);

    (*ctx->mpsolve_ptr)(ctx);
}

static void*
    mps_caller(mps_context* ctx)
{
    if (!mps_context_has_errors(ctx))
    {
        mps_preliminary_setup(ctx);
        ctx->mpsolve_ptr(ctx);
    }

    /* Call user defined callback if available */
    if (ctx->callback == NULL)
        return NULL;
    else
        return (*ctx->callback)(ctx, ctx->user_data);
}

void
    mps_mpsolve_async(mps_context* ctx, mps_callback callback, void* user_data)
{
#ifdef MPS_CATCH_FPE
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

    /* Set up callbacks */
    ctx->callback = callback;
    ctx->user_data = user_data;

    mps_thread_pool* private_pool = mps_thread_pool_new(ctx, 1);
    mps_thread_pool_set_strict_async(private_pool, true);
    ctx->self_thread_pool = private_pool;
    mps_thread_pool_assign(ctx, private_pool, (mps_thread_work)mps_caller, ctx);
}
