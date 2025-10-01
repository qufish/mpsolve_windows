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

/**
* @brief Check if the target set has been reached or not, and update
* the field ctx->root_inclusion[i] for every root.
*/
MPS_PRIVATE void
    mps_fupdate_inclusions(mps_context* ctx)
{
    mps_cluster_item* cluster_item;
    mps_cluster* cluster;
    mps_root* root;
    int i, nf = 2 * ctx->n;

    MPS_DEBUG_THIS_CALL(ctx);

    /* Scan the inclusion depending on the selected search set. */
    for (cluster_item = ctx->clusterization->first_cluster_item; cluster_item != NULL;
        cluster_item = cluster_item->next_cluster_item)
    {
        cluster = cluster_item->cluster;

        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            i = root->k;

            /* First check if the root has already recongnized as part of
            * a set (or out of it) and if that's true skip to the next one. */
            if (ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
                switch (ctx->output_config->search_set)
                {
                case MPS_SEARCH_SET_COMPLEX_PLANE:
                    ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                    break;

                case MPS_SEARCH_SET_UNITARY_DISC:
                    if (!mps_ftouchunit(ctx, nf, i))
                        ctx->approx_root[i]->inclusion = (cplx_mod(ctx->approx_root[i]->fvalue) < 1) ? MPS_ROOT_INCLUSION_IN :
                        MPS_ROOT_INCLUSION_OUT;
                    break;

                case MPS_SEARCH_SET_UNITARY_DISC_COMPL:
                    if (!mps_ftouchunit(ctx, nf, i))
                        ctx->approx_root[i]->inclusion = (cplx_mod(ctx->approx_root[i]->fvalue) > 1) ?
                        MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    break;

                case MPS_SEARCH_SET_NEGATIVE_REAL_PART:
                    if (!mps_ftouchimag(ctx, nf, i))
                        ctx->approx_root[i]->inclusion = (cplx_Re(ctx->approx_root[i]->fvalue) < 0) ?
                        MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    break;

                case MPS_SEARCH_SET_POSITIVE_REAL_PART:
                    if (!mps_ftouchimag(ctx, nf, i))
                        ctx->approx_root[i]->inclusion = (cplx_Re(ctx->approx_root[i]->fvalue) > 0) ?
                        MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    break;

                case MPS_SEARCH_SET_NEGATIVE_IMAG_PART:
                    if (!mps_ftouchreal(ctx, nf, i))
                        ctx->approx_root[i]->inclusion = (cplx_Im(ctx->approx_root[i]->fvalue) < 0) ?
                        MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    break;

                case MPS_SEARCH_SET_POSITIVE_IMAG_PART:
                    if (!mps_ftouchreal(ctx, nf, i))
                        ctx->approx_root[i]->inclusion = (cplx_Im(ctx->approx_root[i]->fvalue) > 0) ?
                        MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    break;

                case MPS_SEARCH_SET_REAL:
                    /* Quite a particular case since we are requiring detection
                    * of the reality of the roots. */
                    if (cluster->n == 1)
                    {
                        if (mps_ftouchreal(ctx, 1, i))
                        {
                            if (MPS_STRUCTURE_IS_REAL(ctx->active_poly->structure) ||
                                (::log(ctx->approx_root[i]->frad) < ctx->sep - ctx->n * ctx->lmax_coeff))
                            {
                                ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                                ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_REAL;
                            }
                        }
                        else
                        {
                            ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                            ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                        }
                    }

                    break;

                case MPS_SEARCH_SET_IMAG:
                    /* The same as the real case, a part from the fact that
                    * we don't support a pure imaginary coefficients polynomial. */
                    if (cluster->n == 1)
                    {
                        if (mps_ftouchimag(ctx, 1, i))
                        {
                            if (::log(ctx->approx_root[i]->frad) < ctx->sep - ctx->n * ctx->lmax_coeff)
                            {
                                ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                                ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_IMAG;
                            }
                            else
                            {
                                ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                                ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                            }
                        }
                    }
                    break;

                case MPS_SEARCH_SET_CUSTOM:
                    break;
                }
        }
    }

    /* Recheck all the clusters and if a cluster with an uncertain root is found reset
    * all the roots in it as uncertaing. */
    for (cluster_item = ctx->clusterization->first_cluster_item; cluster_item != NULL;
        cluster_item = cluster_item->next_cluster_item)
    {
        cluster = cluster_item->cluster;

        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            i = root->k;
            if (ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            {
                for (root = cluster->first_root; root != NULL; root = root->next_root)
                    ctx->approx_root[root->k]->inclusion = MPS_ROOT_INCLUSION_UNKNOWN;
                break;
            }
        }
    }
}

/**
* @brief Check if the target set has been reached or not, and update
* the field ctx->approx_root[i]->inclusion for every root.
*/
MPS_PRIVATE void
    mps_dupdate_inclusions(mps_context* ctx)
{
    mps_cluster_item* cluster_item;
    mps_cluster* cluster;
    mps_root* root;
    int i, nf = 2 * ctx->n;
    rdpe_t mod;

    MPS_DEBUG_THIS_CALL(ctx);

    /* Scan the inclusion depending on the selected search set. */
    for (cluster_item = ctx->clusterization->first_cluster_item; cluster_item != NULL;
        cluster_item = cluster_item->next_cluster_item)
    {
        cluster = cluster_item->cluster;

        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            i = root->k;

            /* First check if the root has already recongnized as part of
            * a set (or out of it) and if that's true skip to the next one. */
            if (ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
                switch (ctx->output_config->search_set)
                {
                case MPS_SEARCH_SET_COMPLEX_PLANE:
                    ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                    break;

                case MPS_SEARCH_SET_UNITARY_DISC:
                    if (!mps_dtouchunit(ctx, nf, i))
                    {
                        cdpe_mod(mod, ctx->approx_root[i]->dvalue);
                        ctx->approx_root[i]->inclusion = (rdpe_le(mod, rdpe_one)) ?
                            MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    }
                    break;

                case MPS_SEARCH_SET_UNITARY_DISC_COMPL:
                    if (!mps_dtouchunit(ctx, nf, i))
                    {
                        cdpe_mod(mod, ctx->approx_root[i]->dvalue);
                        ctx->approx_root[i]->inclusion = (rdpe_ge(mod, rdpe_one)) ?
                            MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    }
                    break;

                case MPS_SEARCH_SET_NEGATIVE_REAL_PART:
                    if (!mps_dtouchimag(ctx, nf, i))
                    {
                        rdpe_set(mod, cdpe_Re(ctx->approx_root[i]->dvalue));
                        ctx->approx_root[i]->inclusion = (rdpe_le(mod, rdpe_zero)) ?
                            MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    }
                    break;

                case MPS_SEARCH_SET_POSITIVE_REAL_PART:
                    if (!mps_dtouchimag(ctx, nf, i))
                    {
                        rdpe_set(mod, cdpe_Re(ctx->approx_root[i]->dvalue));
                        ctx->approx_root[i]->inclusion = (rdpe_ge(mod, rdpe_zero)) ?
                            MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    }
                    break;

                case MPS_SEARCH_SET_NEGATIVE_IMAG_PART:
                    if (!mps_dtouchreal(ctx, nf, i))
                    {
                        rdpe_set(mod, cdpe_Im(ctx->approx_root[i]->dvalue));
                        ctx->approx_root[i]->inclusion = (rdpe_le(mod, rdpe_zero)) ?
                            MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    }
                    break;

                case MPS_SEARCH_SET_POSITIVE_IMAG_PART:
                {
                    rdpe_set(mod, cdpe_Im(ctx->approx_root[i]->dvalue));
                    if (!mps_dtouchreal(ctx, nf, i))
                        ctx->approx_root[i]->inclusion = (rdpe_ge(mod, rdpe_zero)) ?
                        MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                }
                break;

                case MPS_SEARCH_SET_REAL:
                    /* Quite a particular case since we are requiring detection
                    * of the reality of the roots. */
                    if (cluster->n == 1)
                    {
                        if (mps_dtouchreal(ctx, 1, i))
                        {
                            if (MPS_STRUCTURE_IS_REAL(ctx->active_poly->structure) ||
                                (rdpe_log(ctx->approx_root[i]->drad) < ctx->sep - ctx->n * ctx->lmax_coeff))
                            {
                                ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                                ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_REAL;
                            }
                        }
                        else
                        {
                            ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                            ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                        }
                    }

                    break;

                case MPS_SEARCH_SET_IMAG:
                    /* The same as the real case, a part from the fact that
                    * we don't support a pure imaginary coefficients polynomial. */
                    if (cluster->n == 1)
                    {
                        if (mps_dtouchimag(ctx, 1, i))
                        {
                            if (rdpe_log(ctx->approx_root[i]->drad) < ctx->sep - ctx->n * ctx->lmax_coeff)
                            {
                                ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                                ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_IMAG;
                            }
                            else
                            {
                                ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                                ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                            }
                        }
                    }
                    break;

                case MPS_SEARCH_SET_CUSTOM:
                    break;
                }
        }
    }

    /* Recheck all the clusters and if a cluster with an uncertain root is found reset
    * all the roots in it as uncertaing. */
    for (cluster_item = ctx->clusterization->first_cluster_item; cluster_item != NULL;
        cluster_item = cluster_item->next_cluster_item)
    {
        cluster = cluster_item->cluster;

        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            i = root->k;
            if (ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            {
                for (root = cluster->first_root; root != NULL; root = root->next_root)
                    ctx->approx_root[root->k]->inclusion = MPS_ROOT_INCLUSION_UNKNOWN;
                break;
            }
        }
    }
}

/**
* @brief Check if the target set has been reached or not, and update
* the field ctx->approx_root[i]->inclusion for every root.
*/
MPS_PRIVATE void
    mps_mupdate_inclusions(mps_context* ctx)
{
    mps_cluster_item* cluster_item;
    mps_cluster* cluster;
    mps_root* root;
    int i, nf = 2 * ctx->n;
    cdpe_t cmod;
    rdpe_t mod;

    MPS_DEBUG_THIS_CALL(ctx);

    /* Scan the inclusion depending on the selected search set. */
    for (cluster_item = ctx->clusterization->first_cluster_item; cluster_item != NULL;
        cluster_item = cluster_item->next_cluster_item)
    {
        cluster = cluster_item->cluster;

        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            i = root->k;

            /* Get a CDPE representation of ctx->approx_root[i]->mvalue */
            mpc_get_cdpe(cmod, ctx->approx_root[i]->mvalue);

            /* First check if the root has already recongnized as part of
            * a set (or out of it) and if that's true skip to the next one. */
            if (ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
                switch (ctx->output_config->search_set)
                {
                case MPS_SEARCH_SET_COMPLEX_PLANE:
                    ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                    break;

                case MPS_SEARCH_SET_UNITARY_DISC:
                    if (!mps_mtouchunit(ctx, nf, i))
                    {
                        cdpe_mod(mod, cmod);
                        ctx->approx_root[i]->inclusion = (rdpe_le(mod, rdpe_one)) ?
                            MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    }
                    break;

                case MPS_SEARCH_SET_UNITARY_DISC_COMPL:
                    if (!mps_mtouchunit(ctx, nf, i))
                    {
                        cdpe_mod(mod, cmod);
                        ctx->approx_root[i]->inclusion = (rdpe_ge(mod, rdpe_one)) ?
                            MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    }
                    break;

                case MPS_SEARCH_SET_NEGATIVE_REAL_PART:
                    if (!mps_mtouchimag(ctx, nf, i))
                    {
                        rdpe_set(mod, cdpe_Re(cmod));
                        ctx->approx_root[i]->inclusion = (rdpe_le(mod, rdpe_zero)) ?
                            MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    }
                    break;

                case MPS_SEARCH_SET_POSITIVE_REAL_PART:
                    if (!mps_mtouchimag(ctx, nf, i))
                    {
                        rdpe_set(mod, cdpe_Re(cmod));
                        ctx->approx_root[i]->inclusion = (rdpe_ge(mod, rdpe_zero)) ?
                            MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    }
                    break;

                case MPS_SEARCH_SET_NEGATIVE_IMAG_PART:
                    if (!mps_mtouchreal(ctx, nf, i))
                    {
                        rdpe_set(mod, cdpe_Im(cmod));
                        ctx->approx_root[i]->inclusion = (rdpe_le(mod, rdpe_zero)) ?
                            MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                    }
                    break;

                case MPS_SEARCH_SET_POSITIVE_IMAG_PART:
                {
                    rdpe_set(mod, cdpe_Im(cmod));
                    if (!mps_mtouchreal(ctx, nf, i))
                        ctx->approx_root[i]->inclusion = (rdpe_ge(mod, rdpe_zero)) ?
                        MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                }
                break;

                case MPS_SEARCH_SET_REAL:
                    /* Quite a particular case since we are requiring detection
                    * of the reality of the roots. */
                    if (cluster->n == 1)
                    {
                        if (mps_mtouchreal(ctx, 1, i))
                        {
                            if (MPS_STRUCTURE_IS_REAL(ctx->active_poly->structure) ||
                                (rdpe_log(ctx->approx_root[i]->drad) < ctx->sep - ctx->n * ctx->lmax_coeff))
                            {
                                ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                                ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_REAL;
                            }
                        }
                        else
                        {
                            ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                            ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                        }
                    }

                    break;

                case MPS_SEARCH_SET_IMAG:
                    /* The same as the real case, a part from the fact that
                    * we don't support a pure imaginary coefficients polynomial. */
                    if (cluster->n == 1)
                    {
                        if (mps_mtouchimag(ctx, 1, i))
                        {
                            if (rdpe_log(ctx->approx_root[i]->drad) < ctx->sep - ctx->n * ctx->lmax_coeff)
                            {
                                ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                                ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_IMAG;
                            }
                            else
                            {
                                ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                                ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                            }
                        }
                    }
                    break;

                case MPS_SEARCH_SET_CUSTOM:
                    break;
                }
        }
    }

    /* Recheck all the clusters and if a cluster with an uncertain root is found reset
    * all the roots in it as uncertaing. */
    for (cluster_item = ctx->clusterization->first_cluster_item; cluster_item != NULL;
        cluster_item = cluster_item->next_cluster_item)
    {
        cluster = cluster_item->cluster;

        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            i = root->k;
            if (ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            {
                for (root = cluster->first_root; root != NULL; root = root->next_root)
                    ctx->approx_root[root->k]->inclusion = MPS_ROOT_INCLUSION_UNKNOWN;
                break;
            }
        }
    }
}
