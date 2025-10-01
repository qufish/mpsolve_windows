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
#include <float.h>

/**
* @brief Globally set the current precision of mp variables
*
* @param ctx The <code>mps_context</code> of the computation.
* @param prec The precision that is desired for the next MP computations.
*/
MPS_PRIVATE void
    mps_mp_set_prec(mps_context* ctx, long int prec)
{
    long int min_prec = mps_context_get_minimum_precision(ctx);

    ctx->mpwp = (prec / min_prec + 1) * min_prec;
    rdpe_set_2dl(ctx->mp_epsilon, 1.0, -ctx->mpwp);

    if (ctx->debug_level & MPS_DEBUG_MEMORY)
    {
        MPS_DEBUG_RDPE(ctx, ctx->mp_epsilon, "Increased precision to %ld bits. Machine epsilon set to eps", ctx->mpwp);
    }
}

/**
* @brief Allocate all the data needed by MPSolve. Must be called after setting
* the degree of the polynomial (or, more generally, the number of root of the
* equation) in <code>ctx->deg</code>.
*
* @param ctx The <code>mps_context</code> of the computation.
*/
MPS_PRIVATE void
    mps_allocate_data(mps_context* ctx)
{
    MPS_DEBUG_THIS_CALL(ctx);
    int i;

    if (ctx->initialized)
        return;

    ctx->approx_root = mps_newv(mps_approximation*, ctx->n);
    for (i = 0; i < ctx->n; i++)
        ctx->approx_root[i] = mps_approximation_new(ctx);

    /* Reset the cluster structure, so we can start without assumption
    * on the location of the roots. */
    mps_cluster_reset(ctx);

    ctx->order = int_valloc(ctx->deg);

    ctx->fppc1 = cplx_valloc(ctx->deg + 1);

    ctx->mfpc1 = mpc_valloc(ctx->deg + 1);
    for (i = 0; i <= ctx->deg; i++)
        mpc_init2(ctx->mfpc1[i], 0);

    ctx->mfppc1 = mpc_valloc(ctx->deg + 1);
    for (i = 0; i <= ctx->deg; i++)
        mpc_init2(ctx->mfppc1[i], 0);

    /* temporary vectors */
    ctx->spar1 = mps_boolean_valloc(ctx->deg + 2);
    ctx->again_old = mps_boolean_valloc(ctx->deg);

    ctx->fap1 = double_valloc(ctx->deg + 1);
    ctx->fap2 = double_valloc(ctx->deg + 1);

    ctx->dap1 = rdpe_valloc(ctx->deg + 1);
    ctx->dpc1 = cdpe_valloc(ctx->deg + 1);
    ctx->dpc2 = cdpe_valloc(ctx->deg + 1);

    /* Setting some default here, that were not settable because we didn't know
    * the degree of the polynomial */
    for (i = 0; i < ctx->n; i++)
        ctx->approx_root[i]->wp = DBL_DIG * (int)LOG2_10;

    /* Init the mutex that need it */
    mps_mutex_init(ctx->precision_mutex);

    /* Other mt variables in status */
    mps_mutex_init(ctx->data_prec_max.value_mutex);

    ctx->initialized = true;
}

/**
* @brief Raise precision performing a real computation of the data.
*
* @param ctx The <code>mps_context</code> of the computation.
* @param prec The desired precision.
* @return The precision set (that may be different from the one requested
* since GMP works only with precision divisible by some given integer.
*/
MPS_PRIVATE long int
    mps_raise_data(mps_context* ctx, long int prec)
{
    int k;
    mps_polynomial* p = ctx->active_poly;

    /* raise the precision of  mroot */
    for (k = 0; k < ctx->n; k++)
        mpc_set_prec(ctx->approx_root[k]->mvalue, prec);

    /* raise the precision of auxiliary variables */
    for (k = 0; k < ctx->n + 1; k++)
    {
        mpc_set_prec(ctx->mfpc1[k], prec);
        mpc_set_prec(ctx->mfppc1[k], prec);
    }

    mps_polynomial_raise_data(ctx, p, prec);

    return mpc_get_prec(ctx->approx_root[0]->mvalue);
}

/**
* @brief The same of <code>mps_raise_data()</code> but using
* raw routines of GMP, that will not change allocations.
*
* @param ctx The <code>mps_context</code> of the computation.
* @param prec The desired precision.
*/
MPS_PRIVATE void
    mps_raise_data_raw(mps_context* ctx, long int prec)
{
    int k;

    if (!MPS_IS_MONOMIAL_POLY(ctx->active_poly))
        return;

    mps_monomial_poly* p = MPS_MONOMIAL_POLY(ctx->active_poly);

    /* raise the precision of  mroot */
    for (k = 0; k < ctx->n; k++)
        mpc_set_prec_raw(ctx->approx_root[k]->mvalue, prec);

    /* raise the precision of  mfpc */
    if (MPS_IS_MONOMIAL_POLY(ctx->active_poly))
        for (k = 0; k < ctx->n + 1; k++)
            mpc_set_prec_raw(p->mfpc[k], prec);

    /* Raise the precision of sparse vectors */
    if (MPS_DENSITY_IS_SPARSE(ctx->active_poly->density))
        for (k = 0; k < ctx->n; k++)
            if (p->spar[k + 1])
                mpc_set_prec_raw(p->mfppc[k], prec);

    /* raise the precision of auxiliary variables */
    for (k = 0; k < ctx->n + 1; k++)
    {
        mpc_set_prec_raw(ctx->mfpc1[k], prec);
        mpc_set_prec_raw(ctx->mfppc1[k], prec);
    }
}

/**
* @brief Compute the mp_complex values of the coefficients of p(x)
* with the  current precision of mpwds words, given the
* rational or integer coefficients.
*
* @param ctx The <code>mps_context</code> of the computation.
* @param prec The precision that should be set and to which the data should
* be adjusted.
*/
MPS_PRIVATE void
    mps_prepare_data(mps_context* ctx, long int prec)
{
    MPS_DEBUG_THIS_CALL(ctx);

    mps_mutex_lock( ctx->precision_mutex);

    if (ctx->debug_level & MPS_DEBUG_MEMORY)
        MPS_DEBUG(ctx, "Increasing working precision to %ld bits", prec);

    mps_mutex_lock(ctx->data_prec_max.value_mutex);

    if (prec > ctx->data_prec_max.value)
    {
        ctx->data_prec_max.value = mps_raise_data(ctx, prec);
    }
    else
    {
        mps_polynomial_raise_data(ctx, ctx->active_poly, prec);
    }

    mps_mutex_unlock(ctx->data_prec_max.value_mutex);

    mps_mutex_unlock(ctx->precision_mutex);
}

/**
* @brief Resets the data to the highest used precision
*
* @param ctx The <code>mps_context</code> of the computation.
*/
MPS_PRIVATE void
    mps_restore_data(mps_context* ctx)
{
    mps_mutex_lock(ctx->data_prec_max.value_mutex);
    if (ctx->debug_level & MPS_DEBUG_MEMORY)
        MPS_DEBUG(ctx, "Restore data to %ld bits", ctx->data_prec_max.value);

    if (ctx->data_prec_max.value)
    {
        mps_mutex_unlock(ctx->data_prec_max.value_mutex);
        mps_raise_data_raw(ctx, ctx->data_prec_max.value);
    }
    else
        mps_mutex_unlock(ctx->data_prec_max.value_mutex);
}

/**
* @brief Free all the data allocated with <code>mps_allocate_data()</code>
*
* @param ctx The <code>mps_context</code> of the computation.
*/
MPS_PRIVATE void
    mps_free_data(mps_context* ctx)
{
    int i;

    if (ctx->debug_level & MPS_DEBUG_MEMORY)
    {
        MPS_DEBUG(ctx, "Deallocating data");
    }

    if (ctx->bmpc)
    {
        mpc_vclear(ctx->bmpc, ctx->n * ctx->pool->n);
        mps_free(ctx->bmpc);
        ctx->bmpc = NULL;
    }

    /* Release our reference to any active polynomial without
    * freeing it, since that is responsability of the user.
    *
    * Note: This is really important since this context may be
    * recycled at a later time, and having a pointer to a possibly
    * not-anymore-valid polynomial will cause a lot of issues. */
    if (ctx->active_poly)
        ctx->active_poly = NULL;

    mps_clusterization_free(ctx, ctx->clusterization);
    ctx->clusterization = NULL;

    mps_free(ctx->order);

    for (i = 0; i < ctx->n; i++)
        mps_approximation_free(ctx, ctx->approx_root[i]);
    mps_free(ctx->approx_root);

    for (i = 0; i <= ctx->deg; i++)
        mpc_clear(ctx->mfpc1[i]);
    mpc_vfree(ctx->mfpc1);

    cplx_vfree(ctx->fppc1);
    for (i = 0; i <= ctx->deg; i++)
    {
        mpc_clear(ctx->mfppc1[i]);
    }

    mps_mutex_destroy(ctx->data_prec_max.value_mutex);
    mps_mutex_destroy(ctx->precision_mutex);

    mps_free(ctx->mfppc1);

    /* free temporary vectors */
    mps_free(ctx->spar1);
    mps_free(ctx->again_old);

    mps_free(ctx->fap1);
    mps_free(ctx->fap2);

    rdpe_vfree(ctx->dap1);
    cdpe_vfree(ctx->dpc1);
    cdpe_vfree(ctx->dpc2);
}
