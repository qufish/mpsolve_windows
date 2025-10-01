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

#include <mps/gmptools.h>
#include <mps/mps.h>


/**
* @brief Check the validity of the inclusions disks.
*
* @param ctx A pointer to the current mps_context.
*/
MPS_PRIVATE mps_boolean
    mps_inclusion(mps_context* ctx)
{
    int i, j, k, oldnclust;
    rdpe_t rad, difr;
    cdpe_t difc;
    mpc_t tmp;
    rdpe_t ap, az, temp, ep, apeps;
    cdpe_t temp1;
    mpc_t p;
    mps_monomial_poly* poly = MPS_MONOMIAL_POLY(ctx->active_poly);

    /* add inclusion code here */
    if (!ctx->chkrad || ctx->lastphase != mp_phase)
    {
        if (ctx->DOLOG)
            fprintf(ctx->logstr, "Skipping inclusion disks check.\n");
        return true;
    }

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "Checking inclusion disks...\n");

    if (ctx->DOLOG)
    {
        fprintf(ctx->logstr, "Old radii\n");
        for (i = 0; i < ctx->n; i++)
        {
            fprintf(ctx->logstr, "r(%d)=", i);
            rdpe_outln_str(ctx->logstr, ctx->approx_root[i]->drad);
        }
    }

    /* save old radii */
    for (i = 0; i < ctx->n; i++)
        rdpe_set(ctx->dap1[i], ctx->approx_root[i]->drad);

    mpc_init2(p, ctx->mpwp);
    rdpe_mul_d(ep, ctx->mp_epsilon, (double)(ctx->n * 4));

    mpc_init2(tmp, ctx->mpwp);

    for (i = 0; i < ctx->n; i++)
    {
        /* compute denominator */
        rdpe_set(rad, rdpe_one);
        for (j = 0; j < ctx->n; j++)
        {
            if (i == j)
                continue;
            mpc_sub(tmp, ctx->approx_root[j]->mvalue, ctx->approx_root[i]->mvalue);
            mpc_get_cdpe(difc, tmp);
            cdpe_smod(difr, difc);
            rdpe_mul_eq(rad, difr);
        }
        rdpe_sqrt_eq(rad);
        rdpe_mul_eq(rad, poly->dap[ctx->n]);

        /* compute numerator */
        if (MPS_DENSITY_IS_SPARSE(ctx->active_poly->density))
        {                       /* case of sparse polynomial */
        /* compute p(mroot[i]) */
            mps_polynomial_meval(ctx, MPS_POLYNOMIAL(poly), ctx->approx_root[i]->mvalue, p, ap);
            rdpe_div_eq(ap, ctx->mp_epsilon);
        }
        else
        {                       /*  dense polynomial */
        /* commpute p(mroot[i]) and p'(mroot[i]) */
            mpc_set(p, poly->mfpc[ctx->n]);
            for (k = ctx->n - 1; k > 0; k--)
            {
                mpc_mul(p, p, ctx->approx_root[i]->mvalue);
                mpc_add(p, p, poly->mfpc[k]);
            }
            mpc_mul(p, p, ctx->approx_root[i]->mvalue);
            mpc_add(p, p, poly->mfpc[0]);

            /* compute bound to the error */
            rdpe_set(ap, poly->dap[ctx->n]);
            mpc_get_cdpe(temp1, ctx->approx_root[i]->mvalue);
            cdpe_mod(az, temp1);
            for (k = ctx->n - 1; k >= 0; k--)
            {
                rdpe_mul(temp, ap, az);
                rdpe_add(ap, temp, poly->dap[k]);
            }
        }

        /* common part */
        mpc_get_cdpe(difc, p);
        cdpe_mod(difr, difc);
        rdpe_mul(apeps, ap, ep);
        rdpe_add_eq(apeps, difr);
        rdpe_mul_eq_d(apeps, (double)ctx->n);

        /* compute ratio */
        rdpe_div(ctx->approx_root[i]->drad, apeps, rad);

        if (ctx->DOLOG)
        {
            fprintf(ctx->logstr, "New r(%d)=", i);
            rdpe_outln_str(ctx->logstr, ctx->approx_root[i]->drad);
        }
    }

    oldnclust = ctx->clusterization->n;

    rdpe_t* newton_radii = rdpe_valloc(ctx->n);
    for (i = 0; i < ctx->n; i++)
        rdpe_set(newton_radii[i], ctx->approx_root[i]->drad);

    mps_mcluster(ctx, newton_radii, 2 * ctx->n);
    mps_free(newton_radii);

    if (ctx->clusterization->n >= oldnclust)
    {
        /* choose the smallest radius */
        for (i = 0; i < ctx->n; i++)
            if (rdpe_lt(ctx->dap1[i], ctx->approx_root[i]->drad))
                rdpe_set(ctx->approx_root[i]->drad, ctx->dap1[i]);
        /* update(); */
    }
    else
        mps_warn(ctx, "Some roots might be not approximated");

    mpc_clear(tmp);
    mpc_clear(p);

    return true;
}
