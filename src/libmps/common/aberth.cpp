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
* @brief Compute Aberth correction for j-th root, without
* selective correction.
*/
MPS_PRIVATE void
    mps_faberth(mps_context* ctx, mps_approximation* root, cplx_t abcorr)
{
    cplx_t z;
    int i;

    cplx_set(abcorr, cplx_zero);
    for (i = 0; i < ctx->n; i++)
    {
        if (ctx->approx_root[i] == root)
            continue;

        cplx_sub(z, root->fvalue, ctx->approx_root[i]->fvalue);
        cplx_inv_eq(z);
        cplx_add_eq(abcorr, z);
    }
}

MPS_PRIVATE void
    mps_faberth_wl(mps_context* ctx, int j, cplx_t abcorr, mps_mutex_t* aberth_mutexes)
{
    int i;
    cplx_t z, froot;

    mps_mutex_lock(aberth_mutexes[j]);
    cplx_set(froot, ctx->approx_root[j]->fvalue);
    mps_mutex_unlock(aberth_mutexes[j]);

    cplx_set(abcorr, cplx_zero);
    for (i = 0; i < ctx->n; i++)
    {
        if (i == j)
            continue;

        mps_mutex_lock(aberth_mutexes[i]);
        cplx_sub(z, froot, ctx->approx_root[i]->fvalue);
        mps_mutex_unlock(aberth_mutexes[i]);

        cplx_inv_eq(z);
        cplx_add_eq(abcorr, z);
    }
}


/**
* @brief Compute Aberth correction for j-th root, without
* selective correction.
*/
MPS_PRIVATE void
    mps_daberth(mps_context* ctx, mps_approximation* root, cdpe_t abcorr)
{
    int i;
    cdpe_t z;

    cdpe_set(abcorr, cdpe_zero);
    for (i = 0; i < ctx->n; i++)
    {
        if (ctx->approx_root[i] == root)
            continue;
        cdpe_sub(z, root->dvalue, ctx->approx_root[i]->dvalue);
        cdpe_inv_eq(z);
        cdpe_add_eq(abcorr, z);
    }
}

/**
* @brief Compute Aberth correction for j-th root, without
* selective correction.
*/
MPS_PRIVATE void
    mps_maberth(mps_context* ctx, mps_approximation* root, mpc_t abcorr)
{
    int i;
    cdpe_t z, temp;
    mpc_t diff;

    mpc_init2(diff, ctx->mpwp);

    cdpe_set(temp, cdpe_zero);
    for (i = 0; i < ctx->n; i++)
    {
        if (ctx->approx_root[i] == root)
            continue;
        mpc_sub(diff, root->mvalue, ctx->approx_root[i]->mvalue);
        mpc_get_cdpe(z, diff);
        cdpe_inv_eq(z);
        cdpe_add_eq(temp, z);
    }
    mpc_set_cdpe(abcorr, temp);

    mpc_clear(diff);
}

/**
* @brief Compute Aberth correction for the j-th root,
* but only with other roots of the <code>jc</code>-th
* cluster.
*/
MPS_PRIVATE void
    mps_faberth_s(mps_context* ctx, mps_approximation* ab_root, mps_cluster* cluster, cplx_t abcorr)
{
    cplx_t z;
    mps_root* root;

    cplx_set(abcorr, cplx_zero);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        mps_approximation* appr = ctx->approx_root[root->k];
        if (appr == ab_root)
            continue;
        cplx_sub(z, ab_root->fvalue, appr->fvalue);
        cplx_inv_eq(z);
        cplx_add_eq(abcorr, z);
    }
}

/**
* @brief Compute Aberth correction for the j-th root,
* but only with other roots of the <code>jc</code>-th
* cluster.
*/
MPS_PRIVATE void
    mps_daberth_s(mps_context* ctx, mps_approximation* ab_root, mps_cluster* cluster, cdpe_t abcorr)
{
    mps_root* root;
    cdpe_t z;

    cdpe_set(abcorr, cdpe_zero);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        mps_approximation* appr = ctx->approx_root[root->k];
        if (appr == ab_root)
            continue;
        cdpe_sub(z, ab_root->dvalue, appr->dvalue);
        cdpe_inv_eq(z);
        cdpe_add_eq(abcorr, z);
    }
}

/**
* @brief Compute Aberth correction for the j-th root,
* but only with other roots of the <code>jc</code>-th
* cluster.
*/
MPS_PRIVATE void
    mps_maberth_s(mps_context* ctx, mps_approximation* ab_root, mps_cluster* cluster, mpc_t abcorr)
{
    mps_root* root;
    cdpe_t z, temp;
    mpc_t diff;

    mpc_init2(diff, ctx->mpwp);

    cdpe_set(temp, cdpe_zero);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        mps_approximation* appr = ctx->approx_root[root->k];
        if (appr == ab_root)
            continue;
        mpc_sub(diff, ab_root->mvalue, appr->mvalue);
        mpc_get_cdpe(z, diff);
        cdpe_inv_eq(z);
        cdpe_add_eq(temp, z);
    }
    mpc_set_cdpe(abcorr, temp);

    mpc_clear(diff);
}

MPS_PRIVATE void
    mps_daberth_wl(mps_context* ctx, int j, cdpe_t abcorr, mps_mutex_t* aberth_mutexes)
{
    int i;
    cdpe_t z, droot;

    mps_mutex_lock(aberth_mutexes[j]);
    cdpe_set(droot, ctx->approx_root[j]->dvalue);
    mps_mutex_unlock(aberth_mutexes[j]);

    cdpe_set(abcorr, cdpe_zero);
    for (i = 0; i < ctx->n; i++)
    {
        if (i == j)
            continue;

        mps_mutex_lock(aberth_mutexes[i]);
        cdpe_sub(z, droot, ctx->approx_root[i]->dvalue);
        mps_mutex_unlock(aberth_mutexes[i]);

        cdpe_inv_eq(z);
        cdpe_add_eq(abcorr, z);
    }
}

MPS_PRIVATE void
    mps_maberth_s_wl(mps_context* ctx, int j, mps_cluster* cluster, mpc_t abcorr,
        mps_mutex_t* aberth_mutexes)
{
    mps_root* root;
    cdpe_t z, temp;
    mpc_t diff, mroot;

    mpc_init2(mroot, ctx->mpwp);
    mpc_init2(diff, ctx->mpwp);

    mps_mutex_lock(aberth_mutexes[j]);
    mpc_set(mroot, ctx->approx_root[j]->mvalue);
    mps_mutex_unlock(aberth_mutexes[j]);

    cdpe_set(temp, cdpe_zero);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        int k;
        k = root->k;
        if (k == j)
            continue;

        mps_mutex_lock(aberth_mutexes[k]);
        mpc_sub(diff, mroot, ctx->approx_root[k]->mvalue);
        mps_mutex_unlock(aberth_mutexes[k]);
        mpc_get_cdpe(z, diff);

        if (cdpe_eq_zero(z))
            continue;

        cdpe_inv_eq(z);
        cdpe_add_eq(temp, z);
    }
    mpc_set_cdpe(abcorr, temp);

    mpc_clear(mroot);
    mpc_clear(diff);
}
