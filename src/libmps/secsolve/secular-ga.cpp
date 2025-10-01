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

MPS_PRIVATE mps_boolean
    mps_context_has_floating_point_exceptions(mps_context* ctx)
{
    int i;

    for (i = 0; i < ctx->n; i++)
        if (cplx_check_fpe(ctx->approx_root[i]->fvalue) ||
            isnan(ctx->approx_root[i]->frad) ||
            isinf(ctx->approx_root[i]->frad))
        {
            return true;
        }

    return false;
}

/**
* @brief Update all the coefficients of the secular equation, and their
* moduli, using the recomputed one stored in the multiprecision version.
*
* @param ctx The mps_context of the computation.
*/
MPS_PRIVATE void
    mps_secular_ga_update_coefficients(mps_context* ctx)
{
    int i;
    mps_secular_equation* sec = ctx->secular_equation;

    for (i = 0; i < ctx->n; ++i)
    {
        mpc_get_cplx(sec->afpc[i], sec->ampc[i]);
        mpc_get_cplx(sec->bfpc[i], sec->bmpc[i]);

        mpc_get_cdpe(sec->adpc[i], sec->ampc[i]);
        mpc_get_cdpe(sec->bdpc[i], sec->bmpc[i]);

        cdpe_mod(sec->aadpc[i], sec->adpc[i]);
        cdpe_mod(sec->abdpc[i], sec->bdpc[i]);

        sec->aafpc[i] = cplx_mod(sec->afpc[i]);
        sec->abfpc[i] = cplx_mod(sec->bfpc[i]);
    }
}

/**
* @brief Check if iterations can terminate, i.e. if newton
* isolation has been reached, if the target was approximate.
* If the target was approximation then
* <code>mps_secular_improve ()</code> should be used to reach
* the required precision.
*
* @param ctx The mps_context of the computation.
*/
MPS_PRIVATE mps_boolean
    mps_secular_ga_check_stop(mps_context* ctx)
{
    MPS_DEBUG_THIS_CALL(ctx);

    int i;

    if (ctx->exit_required)
        return true;

    for (i = 0; i < ctx->n; i++)
    {
        switch (ctx->lastphase)
        {
            /* Float case */
        case float_phase:
            if (!MPS_ROOT_STATUS_IS_COMPUTED(ctx->approx_root[i]->status))
            {
                MPS_DEBUG_WITH_INFO(ctx, "Root %d is not isolated, nor approximated, so we can't stop now.", i);
                return false;
            }
            break;

            /* Multiprecision and DPE case are the same, since the radii
            * are always RDPE. */
        case mp_phase:
            if (!MPS_ROOT_STATUS_IS_COMPUTED(ctx->approx_root[i]->status))
            {
                MPS_DEBUG_WITH_INFO(ctx, "Root %d is not isolated, nor approximated, so we can't stop now.", i);
                MPS_DEBUG_WITH_INFO(ctx, "Status of root %d: %s", i, MPS_ROOT_STATUS_TO_STRING(ctx->approx_root[i]->status));
                return false;
            }
            break;

        case dpe_phase:
            MPS_DEBUG(ctx, "Status of root %d: %s", i, MPS_ROOT_STATUS_TO_STRING(ctx->approx_root[i]->status));
            if (!MPS_ROOT_STATUS_IS_COMPUTED(ctx->approx_root[i]->status))
            {
                MPS_DEBUG_WITH_INFO(ctx, "Root %d is not isolated, nor approximated, so we can't stop now.", i);
                return false;
            }
            break;

        default:
            break;
        }
    }

    MPS_DEBUG_WITH_INFO(ctx, "Stop conditions were satisfied");
    return true;
}

#define EXIT_ON_ERRORS(ctx) if (mps_context_has_errors (ctx)) { goto cleanup; }

/**
* @brief MPSolve main function for the secular equation solving
* using Gemignani's approach.
*
* @param ctx The mps_context of the computation.
*/
MPS_PRIVATE void
    mps_secular_ga_mpsolve(mps_context* ctx)
{
    int roots_computed = 0;
    int packet;
    int i;
    mps_boolean skip_check_stop = false;
    mps_boolean just_regenerated = false;
    mps_secular_equation* sec = mps_secular_equation_from_status(ctx);

    /* Deflate polynomial before starting, in case it's possible */
    MPS_DEBUG_WITH_INFO(ctx, "Checking if the input polynomial needs to be deflated");
    if (MPS_IS_MONOMIAL_POLY(ctx->active_poly))
        mps_monomial_poly_deflate(ctx, ctx->active_poly);

    /* Check if the secular equation is allocated or if only the
    * polynomial is present. In the last case, allocate an empty
    * secular equation to hold the data during the computation. */
    if (!sec)
    {
        ctx->secular_equation = mps_secular_equation_new_raw(ctx, ctx->n);
        sec = mps_secular_equation_from_status(ctx);
    }

    mps_allocate_data(ctx);

    ctx->just_raised_precision = true;

#ifndef DISABLE_DEBUG
    /* Reset all time counters */
    ctx->regeneration_time = 0;
    ctx->fp_iteration_time = 0;
    ctx->dpe_iteration_time = 0;
    ctx->mp_iteration_time = 0;
    clock_t* total_clock = mps_start_timer();
#endif

    /* Set the output desired for the output */
    rdpe_set_2dl(ctx->eps_out, 1.0, -ctx->output_config->prec);

    /* Set degree and allocate polynomial-related variables
    * to allow initializitation to be performed. */
    ctx->deg = ctx->n = ctx->active_poly->degree;

    MPS_DEBUG(ctx, "Degree = %d", ctx->deg);

    /* Manually set FILE* pointer for streams.
    * More refined options will be added later. */
    packet = 0;

    /* Set the maximum possible radius even in DPE, since
    * we may be starting directly from the DPE phase */
    for (i = 0; i < ctx->n; i++)
    {
        ctx->approx_root[i]->frad = DBL_MAX;
        rdpe_set(ctx->approx_root[i]->drad, RDPE_BIG);
        /* rdpe_set_d (ctx->approx_root[i]->drad, DBL_MAX) */
    }

    /* Set initial cluster structure as no cluster structure. */
    mps_cluster_reset(ctx);

    /* Set phase */
    ctx->lastphase = ctx->input_config->starting_phase;

    /* Set the number of roots in and out from the target set. Since we do not
    * support this yet, we can say that all the roots are in the set. */
    ctx->count[0] = ctx->n;
    ctx->count[1] = 0;
    ctx->count[2] = 0;

    /* If the input was polynomial we need to determine the secular
    * coefficients */
    if (!MPS_IS_SECULAR_EQUATION(ctx->active_poly))
    {
        mps_polynomial* p = ctx->active_poly;

        for (i = 0; i < ctx->n; i++)
        {
            cplx_set(sec->bfpc[i], cplx_zero);
            cdpe_set(sec->bdpc[i], cdpe_zero);
        }

        /* Check data first */
        if (ctx->input_config->starting_phase == no_phase)
        {
            char which_case;
            mps_check_data(ctx, &which_case);

            if (mps_context_has_errors(ctx))
            {
#ifndef DISABLE_DEBUG
                mps_stop_timer(total_clock);
#endif
                return;
            }

            MPS_DEBUG_WITH_INFO(ctx, "Check data suggests starting phase should be %s", (which_case == 'f') ? "floating point" : "DPE phase");

            if (which_case == 'f')
                ctx->lastphase = float_phase;
            else
                ctx->lastphase = dpe_phase;
        }
        else
            ctx->lastphase = ctx->input_config->starting_phase;

    preliminary_aberth_packet:

        MPS_DEBUG_WITH_INFO(ctx, "Computing starting points and performing first Aberth packet");

        /* Perform a packet of Aberth iterations */
        switch (ctx->lastphase)
        {
        case float_phase:
            mps_polynomial_fstart(ctx, p, ctx->approx_root);
            EXIT_ON_ERRORS(ctx);

            if (p->fnewton)
                mps_faberth_packet(ctx, p, false);

            /* Check for floating point exceptions. In case we had some
            * at a certain point it is necessary to restart the algorithm
            * in DPE. */
            if (mps_context_has_floating_point_exceptions(ctx))
            {
                /* The cluster reset is needed because at this point
                * we cannot trust the information stored in the cluster
                * analysis, and it is not convenient to do a usual phase
                * transition anyway: we just restart from scratch. */
                mps_cluster_reset(ctx);

                MPS_DEBUG_WITH_INFO(ctx, "Aberth has failed due to floating point exceptions");
                MPS_DEBUG_WITH_INFO(ctx, "The iteration will be restarted using DPE");

                ctx->lastphase = dpe_phase;
                goto preliminary_aberth_packet;
            }

            break;

        case dpe_phase:
            mps_polynomial_dstart(ctx, p, ctx->approx_root);
            EXIT_ON_ERRORS(ctx);

            if (p->dnewton)
                mps_daberth_packet(ctx, p, false);
            break;

        default:
            mps_error(ctx, "Unrecognized starting phase");
#ifndef DISABLE_DEBUG
            mps_stop_timer(total_clock);
#endif
            return;
        }

        EXIT_ON_ERRORS(ctx);

        if (ctx->crude_approximation_mode)
            goto cleanup;

        mps_cluster_analysis(ctx, p);

        if (mps_secular_ga_check_stop(ctx))
            goto cleanup;

        /* In the case where we started in DPE but the initial approximation are
        * representable as standard floating point numbers, go back to float_phase. */
        if (ctx->lastphase == dpe_phase)
        {
            mps_boolean really_need_dpe = false;
            rdpe_t module;

            for (i = 0; i < ctx->n; i++)
            {
                cdpe_mod(module, ctx->approx_root[i]->dvalue);
                if (rdpe_gt(module, rdpe_maxd) || rdpe_lt(module, rdpe_mind))
                    really_need_dpe = true;
            }

            if (!really_need_dpe)
            {
                MPS_DEBUG_WITH_INFO(ctx, "Going back to float_phase because all the approximations "
                    "are representable as standard floating point numbers.");
                ctx->lastphase = float_phase;
                for (i = 0; i < ctx->n; i++)
                {
                    cdpe_get_x(ctx->approx_root[i]->fvalue, ctx->approx_root[i]->dvalue);
                    ctx->approx_root[i]->frad = DBL_MAX;
                }
            }
        }

        /* Check if we can manage to perform the recomputation of the
        * coefficients. If in floating point, switch do DPE if it fail.
        */
        if (!mps_secular_ga_regenerate_coefficients(ctx))
        {
            if (ctx->lastphase == float_phase)
            {
                MPS_DEBUG(ctx, "Switching to DPE phase since initial regeneration of the coefficients did not succeed.");
                ctx->lastphase = dpe_phase;
                mps_polynomial_dstart(ctx, p, ctx->approx_root);

                if (!mps_secular_ga_regenerate_coefficients(ctx))
                {
                    MPS_DEBUG_WITH_INFO(ctx, "Initial generation of the secular equation coefficients did not succeed");
                    mps_error(ctx, "Unable to perform initial regeneration of the secular equation. \n"
                        "This may be caused by a bug in MPSolve or an inadequate precision of the input coefficients.");
                    return;
                }
            }
            else
            {
                MPS_DEBUG_WITH_INFO(ctx, "Initial generation of the secular equation coefficients did not succeed");
                mps_error(ctx, "Unable to perform initial regeneration of the secular equation. \n"
                    "This may be caused by a bug in MPSolve or an inadequate precision of the input coefficients.");
                return;
            }
            just_regenerated = true;
        }
    }
    else
    {
        MPS_DEBUG_WITH_INFO(ctx, "Generated initial coefficients for the secular equation");

        for (i = 0; i < ctx->n; i++)
        {
            cplx_set(ctx->secular_equation->afpc[i], MPS_SECULAR_EQUATION(ctx->active_poly)->afpc[i]);
            cplx_set(ctx->secular_equation->bfpc[i], MPS_SECULAR_EQUATION(ctx->active_poly)->bfpc[i]);
            cdpe_set(ctx->secular_equation->adpc[i], MPS_SECULAR_EQUATION(ctx->active_poly)->adpc[i]);
            cdpe_set(ctx->secular_equation->bdpc[i], MPS_SECULAR_EQUATION(ctx->active_poly)->bdpc[i]);
            mpc_set(ctx->secular_equation->ampc[i], MPS_SECULAR_EQUATION(ctx->active_poly)->ampc[i]);
            mpc_set(ctx->secular_equation->bmpc[i], MPS_SECULAR_EQUATION(ctx->active_poly)->bmpc[i]);

            MPS_POLYNOMIAL(ctx->secular_equation)->degree = ctx->active_poly->degree;
            MPS_POLYNOMIAL(ctx->secular_equation)->structure = MPS_STRUCTURE_COMPLEX_FP;
        }

        ctx->lastphase = float_phase;
        for (i = 0; i < ctx->n; i++)
            ctx->approx_root[i]->wp = 53;
    }

    /* Select initial approximations using the custom secular
    * routine and based on the phase selected by the user. */
    if (!just_regenerated)
    {
        MPS_DEBUG_WITH_INFO(ctx, "Computing starting points");
        switch (ctx->lastphase)
        {
        case  float_phase:
            mps_secular_fstart(ctx, sec, ctx->approx_root);
            break;

        case dpe_phase:
            mps_secular_dstart(ctx, sec, ctx->approx_root);
            break;

        case mp_phase:
            mps_secular_mstart(ctx, sec, ctx->approx_root);
            break;

        default:
            break;
        }
    }

    EXIT_ON_ERRORS(ctx);

    for (i = 0; i < ctx->n; i++)
    {
        ctx->approx_root[i]->again = true;
        ctx->approx_root[i]->approximated = false;
    }

    /* Check if we need to exit */
    if (ctx->exit_required)
    {
        mps_error(ctx, "Exit forced by the caller");
        return;
    }

    /* Cycle until approximated */
    do
    {
        skip_check_stop = false;
        ctx->best_approx = false;

        /* Perform an iteration of floating point Aberth method */
        switch (ctx->lastphase)
        {
        case float_phase:
            MPS_DEBUG_WITH_INFO(ctx, "Starting floating point iterations");

            if (ctx->jacobi_iterations)
                roots_computed = mps_faberth_packet(ctx, MPS_POLYNOMIAL(sec), just_regenerated);
            else
                roots_computed = mps_secular_ga_fiterate(ctx, ctx->max_it, just_regenerated);

            /* If the computation fails we need to switch to DPE so do not
            * break here, but continue the cycle. */
            if (roots_computed != -1)
                break;

        case dpe_phase:
            MPS_DEBUG_WITH_INFO(ctx, "Starting DPE iterations");

            if (ctx->jacobi_iterations)
                roots_computed = mps_daberth_packet(ctx, MPS_POLYNOMIAL(sec), just_regenerated);
            else
                roots_computed = mps_secular_ga_diterate(ctx, ctx->max_it, just_regenerated);

            break;

        case mp_phase:
            MPS_DEBUG_WITH_INFO(ctx, "Starting MP iterations");

            if (ctx->jacobi_iterations)
                roots_computed = mps_maberth_packet(ctx, MPS_POLYNOMIAL(sec), just_regenerated);
            else
                roots_computed = mps_secular_ga_miterate(ctx, ctx->max_it, just_regenerated);

            break;

        default:
            break;
        }

        /* Increase the packet counter */
        packet++;

        /* Check if we need to exit */
        if (ctx->exit_required)
        {
            mps_error(ctx, "Exit forced by the caller");
            return;
        }

        /* Check that we haven't passed the maximum number of allowed iterations */
        if (packet > ctx->max_pack)
        {
            mps_error(ctx, "Maximum number of iteration passed. Aborting.");
#ifndef DISABLE_DEBUG
            mps_stop_timer(total_clock);
#endif
            return;
        }

        /* Check if all roots were approximated with the
        * given input precision                      */
        if (!just_regenerated)
        {
            if (mps_secular_ga_check_stop(ctx))
                break;
            else
                skip_check_stop = true;
        }

        /* If the iterations has ended in less than 2 * not_computed_roots iterations
        * and we have just regenerated the coefficients, we should increase precision. */
        if (ctx->best_approx)
        {
            skip_check_stop = false;

            /* If the user wants to avoid multiprecision we should go directly
            * to the exit stage. */
            if (ctx->avoid_multiprecision)
            {
                MPS_DEBUG_WITH_INFO(ctx, "Multiprecision has been manually disabled, jumping to the exit stage");
                goto cleanup;
            }

            /* Going to multiprecision if we're not there yet */
            if (ctx->lastphase != mp_phase)
                mps_secular_switch_phase(ctx, mp_phase);
            else
            {
                /* Raising precision otherwise */
                mps_secular_raise_precision(ctx, 2 * ctx->mpwp);
            }

            /* Check if we need to exit */
            if (ctx->exit_required)
            {
                mps_error(ctx, "Exit forced by the caller");
                return;
            }

            if (MPS_IS_MONOMIAL_POLY(ctx->active_poly))
            {
                MPS_DEBUG(ctx, "Performing restart phase");
                mps_secular_restart(ctx);
            }

            if (!mps_secular_ga_regenerate_coefficients(ctx))
            {
                MPS_DEBUG(ctx, "Regeneration failed");
            }
            else
                just_regenerated = true;

            /* Check if we need to exit */
            if (ctx->exit_required)
            {
                mps_error(ctx, "Exit forced by the caller");
                return;
            }

            /* just_regenerated = true; */
            ctx->best_approx = false;

            /* Set the packet counter to zero, we are restarting */
            packet = 0;
        }

        /* Check if we need to exit */
        if (ctx->exit_required)
        {
            mps_error(ctx, "Exit forced by the caller");
            return;
        }

        /* If we can't stop recompute coefficients in higher precision and
        * continue to iterate, unless the best approximation possible in
        * this precision has been reached. In that case increase the precision
        * of the computation. */
        ctx->just_raised_precision = false;

        /* Check if all the roots are approximated or, if we have done more than 4 packets
        * of iterations without finding all of them, if at least we are near to the result. */
        {
            if (mps_secular_ga_regenerate_coefficients(ctx))
            {
                skip_check_stop = false;
                just_regenerated = true;
            }
            else
            {
                MPS_DEBUG(ctx, "Raising precision because regeneration failed");

                skip_check_stop = false;

                /* Going to multiprecision if we're not there yet */
                if (ctx->lastphase != mp_phase)
                {
                    mps_secular_switch_phase(ctx, mp_phase);
                }
                else
                {
                    /* Raising precision otherwise */
                    mps_secular_raise_precision(ctx, 2 * ctx->mpwp);
                    mps_secular_ga_regenerate_coefficients(ctx);
                }

                /* just_regenerated = true; */
                ctx->best_approx = false;

                /* Set the packet counter to zero, we are restarting */
                packet = 0;
            }

            /* Check if we need to exit */
            if (ctx->exit_required)
            {
                mps_error(ctx, "Exit forced by the caller");
                return;
            }
        }
    } while (skip_check_stop || !mps_secular_ga_check_stop(ctx));

cleanup:

    if (mps_context_has_errors(ctx))
    {
        MPS_DEBUG_WITH_INFO(ctx, "Returning since some errors have been detected");
        ctx->exit_required = true;
    }
    else
    {
        MPS_DEBUG_WITH_INFO(ctx, "Validating the inclusions");
        if (ctx->active_poly->prec > 0)
        {
            mps_validate_inclusions(ctx);
        }

        mps_copy_roots(ctx);
    }

    MPS_DEBUG_WITH_INFO(ctx, "MPSolve reached the final stage, preparing for exiting.");

    if (ctx->exit_required)
    {
#ifndef DISABLE_DEBUG
        MPS_DEBUG_WITH_INFO(ctx, "Time used from MPSolve: %ld ms", mps_stop_timer(total_clock));
#endif
        return;
    }

    /* Eventually improve the roots if approximation is required */
    if (ctx->output_config->goal == MPS_OUTPUT_GOAL_APPROXIMATE)
    {
#ifdef NICE_DEBUG
        clock_t* my_timer = mps_start_timer();
#endif
        mps_improve(ctx);

#ifdef NICE_DEBUG
        long improve_time = mps_stop_timer(my_timer);
        MPS_DEBUG(ctx, "mps_improve took %lu ms", improve_time);
#endif
    }

    /* Recheck the inclusions before exiting. */
    mps_mupdate_inclusions(ctx);

    /* Debug total time taken but only if debug is enabled */
#ifndef DISABLE_DEBUG
    long tot_time = mps_stop_timer(total_clock);
    if (ctx->debug_level & MPS_DEBUG_TIMINGS)
    {
        MPS_DEBUG(ctx, "Time used for regeneration: %ld ms",
            ctx->regeneration_time);
        MPS_DEBUG(ctx, "Time used in floating point iterations: %ld ms",
            ctx->fp_iteration_time);
        MPS_DEBUG(ctx, "Time used in DPE iterations: %ld ms",
            ctx->dpe_iteration_time);
        MPS_DEBUG(ctx, "Time used in multiprecision iterations: %ld ms",
            ctx->mp_iteration_time);
        MPS_DEBUG(ctx, "Total time using MPSolve: %ld ms",
            tot_time);
    }
#endif
}
