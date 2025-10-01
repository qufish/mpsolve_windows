/*
* This file is part of MPSolve 3.2.2
*
* Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
* License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
*
* Authors:
*   Leonardo Robol <robol@mail.dm.unipi.it>
*/

#include <mps/mps.h>

MPS_PRIVATE mps_boolean
    mps_standard_regeneration_driver_update_fsecular_equation(mps_context* ctx,
        mps_polynomial* p,
        mps_approximation** approximations,
        mps_secular_equation* sec)
{
    mps_boolean successful_regeneration = true;
    cplx_t* old_a = NULL, * old_b = NULL;
    cdpe_t* old_db = NULL;
    mpc_t* old_mb = NULL;
    int i;

    ctx->mpwp = MPS_SECULAR_EQUIVALENT_FP_PRECISION;

    /* Allocate old_a and old_b */
    old_mb = mpc_valloc(ctx->n);
    for (i = 0; i < ctx->n; i++)
        mpc_init2(old_mb[i], approximations[i]->wp);

    old_a = cplx_valloc(ctx->n);
    old_b = cplx_valloc(ctx->n);
    old_db = cdpe_valloc(ctx->n);

    /* Copy the old coefficients, and set the new
    * b_i with the current approximationss approximations. */
    for (i = 0; i < ctx->n; i++)
    {
        cplx_set(old_a[i], sec->afpc[i]);
        cplx_set(old_b[i], sec->bfpc[i]);
        cdpe_set_x(old_db[i], old_b[i]);
        mpc_set_cplx(old_mb[i], old_b[i]);
        mpc_set_prec(sec->bmpc[i], ctx->mpwp);
        mpc_set(sec->bmpc[i], approximations[i]->mvalue);
    }

    /* Regeneration */
    if (!(successful_regeneration = mps_secular_ga_regenerate_coefficients_mp(ctx, old_db, old_mb)))
    {
        for (i = 0; i < ctx->n; i++)
        {
            cplx_set(sec->afpc[i], old_a[i]);
            cplx_set(sec->bfpc[i], old_b[i]);
        }
    }
    else
    {
        mps_secular_ga_update_coefficients(ctx);
        for (i = 0; i < ctx->n; i++)
        {
            /* We may risk that NaN or inf have been introduced because of huge
            * coefficients computed, so let's check it and in the case of failure
            * switch to DPE. */
            if (cplx_check_fpe(sec->afpc[i]) || cplx_check_fpe(sec->bfpc[i]) ||
                (cplx_mod(sec->afpc[i]) > 1.0e300) ||
                (cplx_mod(sec->bfpc[i]) > 1.0e300))
            {
                successful_regeneration = false;
                if (ctx->debug_level & MPS_DEBUG_REGENERATION)
                    MPS_DEBUG(ctx, "Found floating point exception in regenerated coefficients, reusing old ones.");

                for (i = 0; i < ctx->n; i++)
                {
                    cplx_set(sec->afpc[i], old_a[i]);
                    cplx_set(sec->bfpc[i], old_b[i]);
                }
                break;
            }

            if (ctx->debug_level & MPS_DEBUG_REGENERATION)
            {
                MPS_DEBUG_CPLX(ctx, sec->afpc[i], "sec->afpc[%d]", i);
                MPS_DEBUG_CPLX(ctx, sec->bfpc[i], "sec->bfpc[%d]", i);
            }

            mpc_set_cplx(approximations[i]->mvalue, approximations[i]->fvalue);
        }
    }

    cplx_vfree(old_a);
    cplx_vfree(old_b);
    cdpe_vfree(old_db);

    mpc_vclear(old_mb, ctx->n);
    mpc_vfree(old_mb);

    return successful_regeneration;
}

MPS_PRIVATE mps_boolean
    mps_standard_regeneration_driver_update_dsecular_equation(mps_context* ctx,
        mps_polynomial* p,
        mps_approximation** approximations,
        mps_secular_equation* sec)
{
    mps_boolean successful_regeneration = true;
    int i;
    cdpe_t* old_da = NULL;
    cdpe_t* old_db = NULL;
    mpc_t* old_mb = NULL;

    old_mb = mpc_valloc(ctx->n);
    for (i = 0; i < ctx->n; i++)
        mpc_init2(old_mb[i], approximations[i]->wp);

    ctx->mpwp = MPS_SECULAR_EQUIVALENT_FP_PRECISION;

    /* Allocate old_a and old_b */
    old_da = cdpe_valloc(ctx->n);
    old_db = cdpe_valloc(ctx->n);

    /* Copy the old coefficients, and set the new
    * b_i with the current approximationss approximations. */
    for (i = 0; i < ctx->n; i++)
    {
        cdpe_set(old_da[i], sec->adpc[i]);
        cdpe_set(old_db[i], sec->bdpc[i]);
        mpc_get_cdpe(sec->bdpc[i], approximations[i]->mvalue);
        mpc_set_cdpe(old_mb[i], old_db[i]);
        mpc_set_prec(sec->bmpc[i], ctx->mpwp);
        mpc_set(sec->bmpc[i], approximations[i]->mvalue);
    }

    /* Regeneration */
    if (!(successful_regeneration = mps_secular_ga_regenerate_coefficients_mp(ctx, old_db, old_mb)))
    {
        MPS_DEBUG(ctx, "Regeneration failed");
        for (i = 0; i < ctx->n; i++)
        {
            cdpe_set(sec->adpc[i], old_da[i]);
            cdpe_set(sec->bdpc[i], old_db[i]);
            mpc_set_cdpe(old_mb[i], old_db[i]);
            mpc_set_cdpe(sec->ampc[i], old_da[i]);
            mpc_set_cdpe(sec->bmpc[i], old_db[i]);
        }

        mps_secular_ga_update_coefficients(ctx);
        successful_regeneration = false;
        goto cleanup;
    }
    else
    {
        mps_secular_ga_update_coefficients(ctx);
        for (i = 0; i < ctx->n; i++)
            mpc_set_cdpe(approximations[i]->mvalue, approximations[i]->dvalue);
        mps_secular_set_radii(ctx);
    }

    if (ctx->debug_level & MPS_DEBUG_REGENERATION)
    {
        for (i = 0; i < ctx->n; i++)
        {
            MPS_DEBUG_CDPE(ctx, sec->bdpc[i], "sec->bdpc[%d]", i);
            MPS_DEBUG_CDPE(ctx, sec->adpc[i], "sec->adpc[%d]", i);
        }
    }

    /* Free data */
cleanup:

    cdpe_vfree(old_da);
    cdpe_vfree(old_db);

    mpc_vclear(old_mb, MPS_POLYNOMIAL(sec)->degree);
    mpc_vfree(old_mb);

    return successful_regeneration;
}

MPS_PRIVATE mps_boolean
    mps_standard_regeneration_driver_update_msecular_equation(mps_context* ctx,
        mps_polynomial* p,
        mps_approximation** approximations,
        mps_secular_equation* sec)
{
    mps_boolean successful_regeneration = true;
    int i;
    mpc_t* old_ma = NULL, * old_mb = NULL;
    cdpe_t* old_db = NULL;

    /* Allocate old_a and old_b */
    old_ma = mpc_valloc(ctx->n);
    old_mb = mpc_valloc(ctx->n);
    old_db = cdpe_valloc(ctx->n);

    mpc_vinit2(old_ma, ctx->n, ctx->mpwp);
    mpc_vinit2(old_mb, ctx->n, ctx->mpwp);

    /* Copy the old coefficients, and set the new
    * b_i with the current roots approximations. */
    for (i = 0; i < ctx->n; i++)
    {
        mpc_set(old_ma[i], sec->ampc[i]);
        mpc_set(old_mb[i], sec->bmpc[i]);
        mpc_set_prec(sec->bmpc[i], mpc_get_prec(ctx->approx_root[i]->mvalue));
        mpc_set(sec->bmpc[i], ctx->approx_root[i]->mvalue);
        mpc_get_cdpe(old_db[i], old_mb[i]);
    }

    /* Regeneration */
    if ((successful_regeneration = mps_secular_ga_regenerate_coefficients_mp(ctx, old_db, old_mb)))
    {
        mps_secular_ga_update_coefficients(ctx);
        /* Finally set radius according to new computed a_i coefficients,
        * if they are convenient   */
        mps_secular_set_radii(ctx);
    }
    else
        MPS_DEBUG(ctx, "Regeneration failed");

    if (ctx->debug_level & MPS_DEBUG_REGENERATION)
    {
        MPS_DEBUG(ctx, "Dumping regenerated coefficients");
        for (i = 0; i < ctx->n; i++)
        {
            MPS_DEBUG_MPC(ctx, ctx->mpwp, sec->ampc[i], "ampc[%d]", i);
            MPS_DEBUG_MPC(ctx, ctx->mpwp, sec->bmpc[i], "bmpc[%d]", i);
        }
    }

    mpc_vclear(old_ma, ctx->n);
    mpc_vfree(old_ma);
    rdpe_vfree(old_db);

    return successful_regeneration;
}

static mps_regeneration_driver _mps_standard_regeneration_driver_instance = {
mps_standard_regeneration_driver_update_fsecular_equation,
mps_standard_regeneration_driver_update_dsecular_equation,
mps_standard_regeneration_driver_update_msecular_equation,
NULL
};

mps_regeneration_driver*
    mps_regeneration_driver_new_standard(mps_context* ctx)
{
    return &_mps_standard_regeneration_driver_instance;
}

void
    mps_regeneration_driver_free(mps_context* ctx, mps_regeneration_driver* rd)
{
    if (rd->freeres)
        rd->freeres(ctx, rd);

    mps_free(rd);
}
