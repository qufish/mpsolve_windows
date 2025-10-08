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


#include <stdio.h>
#include <float.h>
#include <mps/mps.h>
#include <string.h>

MPS_PRIVATE void
    mps_set_default_values(mps_context* ctx)
{
    /* flags */
    ctx->skip_float = false;        /* set to true to skip float phase     */
    ctx->resume = false;            /* resume from pre-computed roots      */
    ctx->chkrad = false;            /* check radii after completion        */

    /* I/O flags */
    ctx->DOLOG = false;             /* if nonzero enables the logging      */
    ctx->DOWARN = true;             /* if nonzero enables the warnings     */
    ctx->DOSORT = true;             /* if nonzero enables root sorting     */
    ctx->debug_level = 0;

    /* I/O streams */
    ctx->instr = stdin;              /* input stream                        */
    ctx->outstr = stdout;             /* output stream                       */
    ctx->logstr = stderr;             /* log stream                          */
    ctx->rtstr = NULL;              /* root stream                         */

    /* constants/parameters */
    ctx->max_pack = 100000;           /* number of max packets of iterations */
    ctx->max_it = 20;                /* number of max iterations per packet */
    ctx->max_newt_it = 15;           /* number of max newton iterations for */
    ctx->jacobi_iterations = false;

    /* Set number of threads to 1.5 * number_of_cores, if this is
    * computable. Set it to 12 otherwise.                     */
    ctx->n_threads = (int)1.5 * mps_thread_get_core_number(ctx);
    if (!ctx->n_threads)
        ctx->n_threads = 12;

    ctx->clusterization = NULL;

    ctx->mpwp_max = 100000000;     /* maximum allowed bits for mp         */

    ctx->zero_roots = 0;

    /* soution related variables */
    ctx->lastphase = no_phase;      /* store last computed phase           */
    ctx->order = NULL;              /* output index order: ord[0],..,ord[n] */
    ctx->fppc1 = NULL;              /* standard complex coefficients       */
    ctx->dpc1 = NULL;               /* dpe complex coefficients            */
    ctx->dpc2 = NULL;               /* dpe complex coefficients            */
    ctx->mfpc1 = NULL;              /* temp multiprec. complex coeff. of p' */
    ctx->mfppc1 = NULL;             /* temp multiprec. complex coeff.      */

    ctx->spar1 = NULL;              /* temp sparsity structure of poly     */
    ctx->oldpunt = NULL;            /* stores the previous value of punt   */
    ctx->fap1 = NULL;               /* moduli of the coefficients as double */
    ctx->fap2 = NULL;               /* temp. log of the coeffs as double   */
    ctx->dap1 = NULL;               /* temp moduli of the coeffs as dpe    */

    ctx->again_old = NULL;          /* temp flag vector: true where more   */

    ctx->random_seed = 0;
    ctx->newtis = 0;
    ctx->last_sigma = 0.1;

    ctx->secular_equation = NULL;
    ctx->active_poly = NULL;

    /* Input */
    ctx->input_config->starting_phase = no_phase;

    /* Output */
    ctx->output_config->format = MPS_OUTPUT_FORMAT_COMPACT;
    ctx->output_config->prec = (long int)(0.8 * DBL_DIG * (int)LOG2_10);
    ctx->output_config->goal = MPS_OUTPUT_GOAL_ISOLATE;
    ctx->output_config->multiplicity = false;
    ctx->output_config->root_properties = MPS_OUTPUT_PROPERTY_NONE;
    ctx->output_config->search_set = MPS_SEARCH_SET_COMPLEX_PLANE;

    ctx->data_prec_max.value = 53;

    ctx->mpwp = DBL_DIG * (int)LOG2_10;

    /* Setup for the default algorithm */
    ctx->mpsolve_ptr = MPS_MPSOLVE_PTR(mps_standard_mpsolve);
    ctx->algorithm = MPS_ALGORITHM_STANDARD_MPSOLVE;
    ctx->starting_strategy = MPS_STARTING_STRATEGY_DEFAULT;

    /* Allocate the thread_pool used in computations. */
    ctx->pool = mps_thread_pool_new(ctx, 0);

    /* Callbacks for async version */
    ctx->callback = NULL;
    ctx->user_data = NULL;

    /* Error handling */
    ctx->error_state = false;
    ctx->last_error = NULL;

    ctx->over_max = false;

    ctx->bmpc = NULL;

    /* Set a sensible gnuplot_format by default so the user won't see (null) in the
    * output. */
    ctx->gnuplot_format = "points";

    ctx->self_thread_pool = NULL;
    ctx->avoid_multiprecision = false;
    ctx->crude_approximation_mode = false;

    /* Get a standard regeneration driver to use for MPSolve. Note that this functions
    * does not actually allocate anything, but simply provides a pointer to the internal
    * MPSolve implementation of the standard regeneration_driver. This instance must
    * not be freed, since it is statically declared inside
    * secsolve/standard-regeneration-driver.cpp */
    ctx->regeneration_driver = mps_regeneration_driver_new_standard(ctx);
}
