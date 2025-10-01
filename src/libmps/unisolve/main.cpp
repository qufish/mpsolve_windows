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


/**
* @file
* @brief File with the implementation of the driver routines
* for MPSolve.
*/

#include <math.h>
#include <float.h>
#include <mps/mps.h>

/**
* @brief Main routine of the program that implements the algorithm
* in the standard polynomial version.
*
* The program is divided into many parts
* - Check the correctness of data, scale coefficients if
*   needed, and select cases: the variable <code>which_case</code> is
*   <code>'f'</code> or <code>'d'</code>  according to float or dpe case.
* - Call msolve or dsolve according to the value of which_case.
* - Allocate MP variables mfpc, mroot, drad (if needed).
* - Start MPsolve loop
*   - prepare data according to the current precision
*      and to the data_type (density/sparsity/user)
*   - Call msolve with the current precision
* - check for termination
*/
MPS_PRIVATE void
    mps_standard_mpsolve(mps_context* ctx)
{
    int i, nzc;
    char which_case;
    mps_boolean d_after_f, computed;

#ifndef DISABLE_DEBUG
    clock_t* my_timer = mps_start_timer();
#endif

    mps_allocate_data(ctx);

    if (ctx->DOLOG)
        ctx->debug_level |= MPS_DEBUG_TRACE;

    /* == 1 ==  Setup variables, i.e. copy coefficients
        into dpr, dpc and similar. */
    mps_setup(ctx);

    ctx->lastphase = no_phase;
    computed = false;
    ctx->over_max = false;

    /* == 2 ==  Resume from pre-computed roots */
    if (ctx->resume)
    {
        mps_error(ctx, "Resume not supported yet");
#ifndef DISABLE_DEBUG
        mps_stop_timer(my_timer);
#endif
        return;
    }

    /* == 3 ==  Check data and get starting phase */
    if (ctx->skip_float)
        which_case = 'd';
    else
        which_case = 'f';

    /* This variable is true if we need a dpe phase after the
    * float phase */
    d_after_f = false;

    /* Check if a dpe phase is needed and deflate polynomial */
    mps_check_data(ctx, &which_case);

    /* Check for errors in check data */
    if (mps_context_has_errors(ctx))
    {
#ifndef DISABLE_DEBUG
        mps_stop_timer(my_timer);
#endif
        return;
    }

    rdpe_set_2dl(ctx->eps_out, 1.0, -ctx->output_config->prec);

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "Which_case = %c, skip_float= %d\n", which_case,
            ctx->skip_float);

    /* == 4 ==  Float phase */
    if (which_case == 'f')
    {
        if (ctx->DOLOG)
            fprintf(ctx->logstr, "Float phase ...\n");
        mps_fsolve(ctx, &d_after_f);
        ctx->lastphase = float_phase;

        if (ctx->DOLOG)
            mps_dump(ctx);

        computed = mps_check_stop(ctx);
        if (computed && ctx->output_config->goal != MPS_OUTPUT_GOAL_APPROXIMATE)
            goto exit_sub;
        /* stop for COUNT and ISOLATE goals */
    }

    /* == 5 ==  DPE phase */
    if (which_case == 'd' || d_after_f)
    {                           /* DPE phase */
        if (ctx->DOLOG)
            fprintf(ctx->logstr, "DPE phase ...\n");
        /* If we are arriving from a float phase copy the floating points
        * roots approximations in the DPE root approximations. */
        if (d_after_f)
            for (i = 0; i < ctx->n; i++)
            {
                rdpe_set_d(ctx->approx_root[i]->drad, ctx->approx_root[i]->frad);
                cdpe_set_x(ctx->approx_root[i]->dvalue, ctx->approx_root[i]->fvalue);
            }
        ctx->lastphase = dpe_phase;
        mps_dsolve(ctx, d_after_f);

        if (ctx->DOLOG)
            mps_dump(ctx);

        computed = mps_check_stop(ctx);
        if (computed && ctx->output_config->goal != MPS_OUTPUT_GOAL_APPROXIMATE)
            goto exit_sub;
    }

    /* == 6 ==   Allocate MP variables mfpc, mroot, drad, mfppc, mfppc1
    * (the real input case is not implemented yet ) */
    MPS_DEBUG(ctx, "Starting MP phase");

    ctx->lastphase = mp_phase;

    /* ==== 6.1 initialize mp variables */
    mps_mp_set_prec(ctx, 2 * DBL_MANT_DIG);

    /* Prepare data according to the current working precision */
    mps_prepare_data(ctx, ctx->mpwp);

    /* ==== 6.2 set initial values for mp variables */
    for (i = 0; i < ctx->n; i++)
    {
        if (which_case == 'd' || d_after_f)
            mpc_set_cdpe(ctx->approx_root[i]->mvalue, ctx->approx_root[i]->dvalue);
        else
        {
            mpc_set_cplx(ctx->approx_root[i]->mvalue, ctx->approx_root[i]->fvalue);
            rdpe_set_d(ctx->approx_root[i]->drad, ctx->approx_root[i]->frad);
        }
    }
    if (computed && ctx->output_config->goal == MPS_OUTPUT_GOAL_APPROXIMATE)
    {
        MPS_DEBUG(ctx, "Exiting since the approximation are computed and the goal is MPS_OUTPUT_GOAL_APPROXIMATE");
        goto exit_sub;
    }

    MPS_DEBUG(ctx, "ctx->mpwp = %ld, ctx->mpwp_max = %ld", ctx->mpwp, ctx->mpwp_max);
    MPS_DEBUG(ctx, "ctx->input_config->prec = %ld", ctx->active_poly->prec);

    /* == 7 ==  Start MPSolve loop */
    ctx->mpwp = mps_context_get_minimum_precision(ctx);

    /* Poor man GMP - machine precision detection. We need that min_prec is contained
    * in the interval [ DBL_MANT_DIG , 2 * DBL_MANT_DIG ]. This is probably true on most
    * architectures with the instruction above, but we want to be sure. */
    while (ctx->mpwp < DBL_MANT_DIG)
        ctx->mpwp <<= 1;
    while (ctx->mpwp > 2 * DBL_MANT_DIG)
        ctx->mpwp >>= 1;

    while (!computed && ctx->mpwp < ctx->mpwp_max)
    {
        ctx->mpwp *= 2;

        if (ctx->mpwp > ctx->mpwp_max)
        {
            ctx->mpwp = ctx->mpwp_max;
            ctx->over_max = true;
        }

        if (ctx->DOLOG)
            fprintf(ctx->logstr, "MAIN: mp_loop: mpwp=%ld\n", ctx->mpwp);

        /* == 7.1 ==   prepare data according to the current precision */
        mps_mp_set_prec(ctx, ctx->mpwp);
        mps_prepare_data(ctx, ctx->mpwp);

        /* == 7.2 ==   Call msolve with the current precision */
        if (ctx->DOLOG)
            fprintf(ctx->logstr, "MAIN: now call msolve nclust=%ld\n", ctx->clusterization->n);
        mps_msolve(ctx);
        ctx->lastphase = mp_phase;

        /* if (ctx->DOLOG) dump(logstr); */

        if (ctx->DOLOG)
        {                       /* count isolated zeros */
            nzc = 0;
            for (i = 0; i < ctx->n; i++)
            {
                if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_ISOLATED ||
                    ctx->approx_root[i]->status == MPS_ROOT_STATUS_APPROXIMATED)
                    nzc++;
            }
            fprintf(ctx->logstr, "MAIN: isolated %d roots\n", nzc);
            fprintf(ctx->logstr, "MAIN: after msolve check stop\n");
        }

        /* == 7.3 ==  Check the stop condition */
        computed = mps_check_stop(ctx);
        mps_mmodify(ctx, true);

        /* == 7.4 ==  reset the status vector */
        for (i = 0; i < ctx->n; i++)
            if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                ctx->approx_root[i]->status = MPS_ROOT_STATUS_CLUSTERED;
    }

    /* == 8 ==  Check for termination */
    if (!computed)
    {
        if (ctx->over_max)
        {
            ctx->over_max = true;
            /* mps_error (ctx, "Reached the maximum working precision"); */
            MPS_DEBUG(ctx, "Reached the maximum working precision");
            goto exit_sub;
        }
        else
        {
            /* mps_warn (ctx, "Reached the input precision"); */
            MPS_DEBUG(ctx, "Reached the input precision");
            goto exit_sub;
        }
    }

exit_sub:

    /* == 9 ==  Check inclusion disks */
    if (computed && ctx->clusterization->n < ctx->n)
        if (!mps_inclusion(ctx))
        {
            mps_error(ctx, "Unable to compute inclusion disks");
            return;
        }

    /* == 10 ==  Refine roots */
    if (computed && !ctx->over_max && ctx->output_config->goal == MPS_OUTPUT_GOAL_APPROXIMATE)
    {
        ctx->lastphase = mp_phase;
        mps_improve(ctx);
    }

    /* == 11 == Check inclusions */
    /* This step is disabled since it cause problems with the lar* kind of polynomials.
    * To be re-enabled a careful check of the necessary precision to avoid NULL DERIVATIVE
    * warnings should be implemented.
    * if (ctx->active_poly->prec > 0)
    *   mps_validate_inclusions (ctx);
    */

    /* == 12 ==  Restore to highest used precision */
    if (ctx->lastphase == mp_phase)
        mps_restore_data(ctx);

#ifndef DISABLE_DEBUG
    {
        unsigned long time = mps_stop_timer(my_timer);
        MPS_DEBUG(ctx, "Total time using MPSolve: %lu ms", time);
    }
#endif

    /* Finally copy the roots ready for output */
    mps_copy_roots(ctx);
}

/**
* @brief Setup vectors and variables
*/
MPS_PRIVATE void
    mps_setup(mps_context* ctx)
{
    int i;
    mps_polynomial* p = ctx->active_poly;
    mpf_t mptemp;
    mpc_t mptempc;

    if (ctx->DOLOG)
    {
        /* fprintf (ctx->logstr, "Goal      = %5s\n", ctx->goal); */
        /* fprintf (ctx->logstr, "Data type = %3s\n", ctx->data_type); */
        fprintf(ctx->logstr, "Degree    = %d\n", ctx->n);
        fprintf(ctx->logstr, "Input prec.  = %ld digits\n", (long)(ctx->active_poly->prec
            * LOG10_2));
        fprintf(ctx->logstr, "Output prec. = %ld digits\n", (long)(ctx->output_config->prec
            * LOG10_2));
    }

    /* setup temporary vectors */
    if (MPS_IS_MONOMIAL_POLY(p) && MPS_DENSITY_IS_SPARSE(ctx->active_poly->density))
    {
        mps_monomial_poly* mp = MPS_MONOMIAL_POLY(p);

        for (i = 0; i <= p->degree; i++)
        {
            mp->fap[i] = 0.0;
            mp->fpr[i] = 0.0;
            rdpe_set(mp->dap[i], rdpe_zero);
            cplx_set(mp->fpc[i], cplx_zero);
            rdpe_set(mp->dpr[i], rdpe_zero);
            cdpe_set(mp->dpc[i], cdpe_zero);
        }
    }

    /* Indexes of the first (and only) cluster start from
    * 0 and reach n */
    mps_cluster_reset(ctx);

    /* set input and output epsilon */
    rdpe_set_2dl(ctx->eps_in, 1.0, 1 - ctx->active_poly->prec);
    rdpe_set_2dl(ctx->eps_out, 1.0, 1 - ctx->output_config->prec);

    /* precision of each root */
    for (i = 0; i < ctx->n; i++)
        ctx->approx_root[i]->wp = 53;

    /* output order info */
    for (i = 0; i < ctx->n; i++)
        ctx->order[i] = i;

    /* reset root counts */
    ctx->count[0] = ctx->count[1] = ctx->count[2] = 0;

    /* compute DPE approximations */
    if (MPS_IS_MONOMIAL_POLY(p))
    {
        mps_monomial_poly* mp = MPS_MONOMIAL_POLY(p);

        /* init temporary mp variables */
        mpf_init2(mptemp, DBL_MANT_DIG);
        mpc_init2(mptempc, DBL_MANT_DIG);

        /* main loop */
        ctx->skip_float = false;
        for (i = 0; i <= ctx->n; i++)
        {
            if (MPS_DENSITY_IS_SPARSE(ctx->active_poly->density) && !mp->spar[i])
                continue;

            if (MPS_STRUCTURE_IS_REAL(ctx->active_poly->structure))
            {
                if (MPS_STRUCTURE_IS_RATIONAL(ctx->active_poly->structure) ||
                    MPS_STRUCTURE_IS_INTEGER(ctx->active_poly->structure))
                {
                    mpf_set_q(mptemp, mp->initial_mqp_r[i]);
                    mpf_get_rdpe(mp->dpr[i], mptemp);
                    /*#G GMP 2.0.2 bug begin */
                    if (rdpe_sgn(mp->dpr[i]) != mpq_sgn(mp->initial_mqp_r[i]))
                        rdpe_neg_eq(mp->dpr[i]);
                    /*#G GMP bug end */
                }

                if (MPS_STRUCTURE_IS_FP(ctx->active_poly->structure))
                    mpf_get_rdpe(mp->dpr[i], mpc_Re(mp->mfpc[i]));

                cdpe_set_e(mp->dpc[i], mp->dpr[i], rdpe_zero);

                /* compute dap[i] and check for float phase */
                rdpe_abs(mp->dap[i], mp->dpr[i]);
                rdpe_abs(mp->dap[i], mp->dpr[i]);
                if (rdpe_gt(mp->dap[i], rdpe_maxd)
                    || rdpe_lt(mp->dap[i], rdpe_mind))
                    ctx->skip_float = true;
            }
            else if (MPS_STRUCTURE_IS_COMPLEX(ctx->active_poly->structure))
            {
                if (MPS_STRUCTURE_IS_RATIONAL(ctx->active_poly->structure) ||
                    MPS_STRUCTURE_IS_INTEGER(ctx->active_poly->structure))
                {
                    mpc_set_q(mptempc, mp->initial_mqp_r[i], mp->initial_mqp_i[i]);
                    mpc_get_cdpe(mp->dpc[i], mptempc);
                    /*#G GMP 2.0.2 bug begin */
                    if (rdpe_sgn(cdpe_Re(mp->dpc[i])) != mpq_sgn(mp->initial_mqp_r[i]))
                        rdpe_neg_eq(cdpe_Re(mp->dpc[i]));
                    if (rdpe_sgn(cdpe_Im(mp->dpc[i])) != mpq_sgn(mp->initial_mqp_i[i]))
                        rdpe_neg_eq(cdpe_Im(mp->dpc[i]));
                    /*#G GMP bug end */
                }
                else if (MPS_STRUCTURE_IS_FP(ctx->active_poly->structure))
                    mpc_get_cdpe(mp->dpc[i], mp->mfpc[i]);

                /* compute dap[i] */
                cdpe_mod(mp->dap[i], mp->dpc[i]);
                if (rdpe_gt(mp->dap[i], rdpe_maxd)
                    || rdpe_lt(mp->dap[i], rdpe_mind))
                    ctx->skip_float = true;
            }
        }

        /* free temporary mp variables */
        mpf_clear(mptemp);
        mpc_clear(mptempc);

        /* adjust input data type */
        if (MPS_STRUCTURE_IS_FP(ctx->active_poly->structure) && ctx->skip_float)
            ctx->active_poly->structure = MPS_STRUCTURE_IS_REAL(ctx->active_poly->structure) ?
            MPS_STRUCTURE_REAL_BIGFLOAT : MPS_STRUCTURE_COMPLEX_BIGFLOAT;

        /* prepare floating point vectors */
        if (!ctx->skip_float)
            for (i = 0; i <= MPS_POLYNOMIAL(p)->degree; i++)
            {
                if (MPS_DENSITY_IS_SPARSE(ctx->active_poly->density) || !mp->spar[i])
                    continue;
                if (MPS_STRUCTURE_IS_REAL(ctx->active_poly->structure))
                {
                    mp->fpr[i] = rdpe_get_d(mp->dpr[i]);
                    mp->fap[i] = fabs(mp->fpr[i]);
                    cplx_set_d(mp->fpc[i], mp->fpr[i], 0.0);
                }
                else
                {
                    cdpe_get_x(mp->fpc[i], mp->dpc[i]);
                    mp->fap[i] = cplx_mod(mp->fpc[i]);
                }
            }
    }
}

/**
* @brief Check consistency of data and makes some basic adjustments.
*
* This routine check, for example, if there are zero roots in the polynomial
* (i.e. no costant term) and deflates the polynomial if necessary (shifting
* the coefficients).
*
* It sets the value of the parameter <code>which_case</code> to <code>'f'</code>
* if a floating point phase is enough, or to <code>'d'</code> if
* a <code>dpe</code> phase is needed.
*
* @param ctx The <code>mps_context</code> associated with the current computation.
* @param which_case the address of the variable which_case;
*/
MPS_PRIVATE void
    mps_check_data(mps_context* ctx, char* which_case)
{
    rdpe_t min_coeff, max_coeff, tmp;
    mps_monomial_poly* p = NULL;
    int i;

    /* case of user-defined polynomial */
    if (!MPS_IS_MONOMIAL_POLY(ctx->active_poly))
    {
        if (ctx->output_config->multiplicity)
            mps_error(ctx,
                "Multiplicity detection not yet implemented for user polynomial");
        if (ctx->output_config->root_properties)
            mps_error(ctx,
                "Real/imaginary detection not yet implemented for user polynomial");
        *which_case = 'd';
        return;
    }
    else
        p = MPS_MONOMIAL_POLY(ctx->active_poly);

    /* Check consistency of input */
    if (rdpe_eq(p->dap[ctx->n], rdpe_zero))
    {
        mps_warn(ctx, "The leading coefficient is zero");
        do
            (ctx->n)--;
        while (rdpe_eq(p->dap[ctx->n], rdpe_zero));

        MPS_POLYNOMIAL(p)->degree = ctx->n;
    }

    /* Compute min_coeff */
    if (rdpe_lt(p->dap[0], p->dap[ctx->n]))
        rdpe_set(min_coeff, p->dap[0]);
    else
        rdpe_set(min_coeff, p->dap[ctx->n]);

    /* Compute max_coeff and its logarithm */
    rdpe_set(max_coeff, p->dap[0]);
    for (i = 1; i <= ctx->n; i++)
        if (rdpe_lt(max_coeff, p->dap[i]))
            rdpe_set(max_coeff, p->dap[i]);
    ctx->lmax_coeff = rdpe_log(max_coeff);

    /*  Multiplicity and sep */
    if (ctx->output_config->multiplicity)
    {
        if (MPS_STRUCTURE_IS_INTEGER(ctx->active_poly->structure))
        {
            mps_compute_sep(ctx);
        }
        else if (MPS_STRUCTURE_IS_RATIONAL(ctx->active_poly->structure))
        {
            mps_warn(ctx, "The multiplicity option has not been yet implemented");
            ctx->sep = 0.0;
        }
        else
        {
            mps_warn(ctx, "The input polynomial has neither integer nor rational");
            mps_warn(ctx, " coefficients: unable to compute multiplicities");
            ctx->sep = 0.0;
        }
    }

    /* Real/Imaginary detection */
    if (ctx->output_config->root_properties ||
        ctx->output_config->search_set == MPS_SEARCH_SET_REAL ||
        ctx->output_config->search_set == MPS_SEARCH_SET_IMAG)
    {
        if (MPS_STRUCTURE_IS_INTEGER(ctx->active_poly->structure))
        {
            mps_compute_sep(ctx);
        }
        else if (MPS_STRUCTURE_IS_RATIONAL(ctx->active_poly->structure))
        {
            mps_error(ctx,
                "The real/imaginary option has not been yet implemented for rational input");
            return;
        }
        else
        {
            mps_error(ctx, "The input polynomial has neither integer nor rational "
                "coefficients: unable to perform real/imaginary options");
            return;
        }
    }

    /* Select cases (dpe or floating point)
    * First normalize the polynomial (only the float version) */
    rdpe_div(tmp, max_coeff, min_coeff);
    rdpe_mul_eq_d(tmp, (double)(ctx->n + 1));
    rdpe_mul_eq(tmp, rdpe_mind);
    rdpe_div_eq(tmp, rdpe_maxd);

    if (rdpe_lt(tmp, rdpe_one))
    {
        mpc_t m_min_coeff;
        cdpe_t c_min_coeff;

        /* if  (n+1)*max_coeff/min_coeff < dhuge/dtiny -  float case */
        *which_case = 'f';
        rdpe_mul_eq(min_coeff, max_coeff);
        rdpe_mul(tmp, rdpe_mind, rdpe_maxd);
        rdpe_div(min_coeff, tmp, min_coeff);
        rdpe_sqrt_eq(min_coeff);

        rdpe_set(cdpe_Re(c_min_coeff), min_coeff);
        rdpe_set(cdpe_Im(c_min_coeff), rdpe_zero);

        mpc_init2(m_min_coeff, mpc_get_prec(p->mfpc[0]));
        mpc_set_cdpe(m_min_coeff, c_min_coeff);

        /* min_coeff = ::sqrt(dhuge*dtiny/(min_coeff*max_coeff))
        * NOTE: This is enabled for floating point polynomials only
        * for the moment, but it may work nicely also for other representations. */
        {
            for (i = 0; i <= ctx->n; i++)
            {
                /* Multiply the MP leading coefficient */
                mpc_mul_eq(p->mfpc[i], m_min_coeff);

                rdpe_mul(tmp, p->dap[i], min_coeff);
                rdpe_set(p->dap[i], tmp);
                p->fap[i] = rdpe_get_d(tmp);

                mpc_get_cdpe(p->dpc[i], p->mfpc[i]);
                cdpe_get_x(p->fpc[i], p->dpc[i]);
            }
        }

        mpc_clear(m_min_coeff);
    }
    else
        *which_case = 'd';
}

/**
* @brief Compute the minimum distance that can separate two roots of the input
* polynomial.
*/
MPS_PRIVATE void
    mps_compute_sep(mps_context* ctx)
{
    ctx->sep = ctx->n * ctx->lmax_coeff;
    ctx->sep = -(ctx->sep) - ctx->n * (1 + log((double)ctx->n)) / LOG2;
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "Sep = %f\n", ctx->sep);
}
