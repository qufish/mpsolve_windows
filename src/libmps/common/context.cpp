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

#define _GNU_SOURCE

#include <mps/mps.h>
#include <string.h>

long int
    mps_context_get_minimum_precision(mps_context* ctx)
{
    return ctx->minimum_gmp_precision;
}

long int
    mps_context_get_data_prec_max(mps_context* ctx)
{
    long int ret;

    mps_mutex_lock(ctx->data_prec_max.value_mutex);
    ret = ctx->data_prec_max.value;
    mps_mutex_unlock(ctx->data_prec_max.value_mutex);
    return ret;
}

int
    mps_context_get_degree(mps_context* ctx)
{
    return ctx->n;
}

/**
* @brief Select algorithm to use for computation.
*
* Valid values for this field are
* - MPS_ALGORITHM_STANDARD_MPSOLVE for the standard MPSolve algorithm;
* - MPS_ALGORITHM_SECULAR_GA for the algorithm based on coefficient regeneration
*   applied to secular equations;
*/
void
    mps_context_select_algorithm(mps_context* ctx, mps_algorithm algorithm)
{
    /* First set algorithm in the mps_context */
    ctx->algorithm = algorithm;

    switch (algorithm)
    {
    case MPS_ALGORITHM_STANDARD_MPSOLVE:
        ctx->mpsolve_ptr = MPS_MPSOLVE_PTR(mps_standard_mpsolve);
        break;

    case MPS_ALGORITHM_SECULAR_GA:
        ctx->mpsolve_ptr = MPS_MPSOLVE_PTR(mps_secular_ga_mpsolve);
        break;
    }
}

/**
* @brief Select the starting strategy used to dispose the initial approximations.
*/
void
    mps_context_select_starting_strategy(mps_context* ctx, mps_starting_strategy strategy)
{
    ctx->starting_strategy = strategy;
}

static void
    mps_context_init(mps_context* ctx)
{
    mpf_t test;

    /* Set default streams */
    ctx->instr = stdin;
    ctx->outstr = stdout;
    ctx->logstr = stdout;

    /* Allocate space for the configurations */
    ctx->input_config = (mps_input_configuration*)mps_malloc(sizeof(mps_input_configuration));
    ctx->output_config = (mps_output_configuration*)mps_malloc(sizeof(mps_output_configuration));

    mps_set_default_values(ctx);

    /* Find minimum GMP supported precision */
    mpf_init2(test, 1);
    ctx->minimum_gmp_precision = (int)mpf_get_prec(test);
    mpf_clear(test);

    /* Set standard precision */
    ctx->output_config->prec = (int)(0.9 * DBL_DIG * (int)LOG2_10);
    MPS_DEBUG(ctx, "Setting prec_out to %ld digits", ctx->output_config->prec);

    mps_mp_set_prec(ctx, DBL_DIG * (int)LOG2_10 + 1);

    ctx->initialized = false;
    ctx->exit_required = false;
}

/**
* @brief Allocate a new mps_context struct with default
* options.
*/
mps_context*
    mps_context_new()
{
    mps_context* ctx = NULL;

    mps_new_defined_obj(mps_context, ctx, sizeof(mps_context));
    mps_context_init(ctx);

    return ctx;
}

/**
* @brief Free a not more useful mps_context.
*
* @param ctx the mps_context struct pointer to free.
*/
void
    mps_context_free(mps_context* ctx)
{
    /* Close input and output streams if they're not stdin, stdout and
    * stderr. For the case in which this context will re-used, set them
    * to their default values. */
    if (ctx->instr != stdin && ctx->instr != NULL)
        fclose(ctx->instr);
    if (ctx->logstr != stderr && ctx->logstr != stdout && ctx->logstr != NULL)
        fclose(ctx->logstr);

    ctx->instr = stdin;
    ctx->logstr = stderr;

    /* There's no need to resize bmpc since they will be allocated on demand.
    * We free them here to correct bad assumptions on the size of this
    * vector. */
    mps_free(ctx->bmpc);
    ctx->bmpc = NULL;

    if (ctx->initialized)
        mps_free_data(ctx);

    mps_thread_pool_free(ctx, ctx->pool);

    mps_free(ctx->input_config);
    mps_free(ctx->output_config);

    ctx->active_poly = NULL;

    if (ctx->secular_equation)
        mps_secular_equation_free(ctx, MPS_POLYNOMIAL(ctx->secular_equation));

    mps_del_obj(ctx);
}

void
    mps_context_abort(mps_context* ctx)
{
    ctx->exit_required = true;
}

static void
    mps_context_shrink(mps_context* ctx, int n)
{
    int i;

    for (i = n; i < ctx->n - ctx->zero_roots; i++)
    {
        mps_approximation_free(ctx, ctx->approx_root[i]);
    }

    ctx->approx_root = (mps_approximation**)mps_realloc(ctx->approx_root, sizeof(mps_approximation*) * n);

    ctx->order = (int*)mps_realloc(ctx->order, sizeof(int) * n);

    ctx->fppc1 = (cplx_t*)mps_realloc(ctx->fppc1, sizeof(cplx_t) * (n + 1));

    for (i = n + 1; i <= ctx->n - ctx->zero_roots; i++)
        mpc_clear(ctx->mfpc1[i]);

    ctx->mfpc1 = (mpc_t*)mps_realloc(ctx->mfpc1, sizeof(mpc_t) * (n + 1));

    for (i = n + 1; i <= ctx->n - ctx->zero_roots; i++)
        mpc_clear(ctx->mfppc1[i]);

    ctx->mfppc1 = (mpc_t*)mps_realloc(ctx->mfppc1, sizeof(mpc_t) * (n + 1));

    /* temporary vectors */
    ctx->spar1 = (mps_boolean*)mps_realloc(ctx->spar1, sizeof(mps_boolean) * (n + 2));
    ctx->again_old = (mps_boolean*)mps_realloc(ctx->again_old, sizeof(mps_boolean) * (n));

    ctx->fap1 = (double*)mps_realloc(ctx->fap1, sizeof(double) * (n + 1));
    ctx->fap2 = (double*)mps_realloc(ctx->fap2, sizeof(double) * (n + 1));

    ctx->dap1 = (rdpe_t*)mps_realloc(ctx->dap1, sizeof(rdpe_t) * (n + 1));
    ctx->dpc1 = (cdpe_t*)mps_realloc(ctx->dpc1, sizeof(cdpe_t) * (n + 1));
    ctx->dpc2 = (cdpe_t*)mps_realloc(ctx->dpc2, sizeof(cdpe_t) * (n + 1));

    /* Setting some default here, that were not settable because we didn't know
    * the degree of the polynomial */
    for (i = 0; i < n; i++)
        ctx->approx_root[i]->wp = DBL_DIG * (int)LOG2_10;
}

static void
    mps_context_expand(mps_context* ctx, int n)
{
    int i;
    long int previous_prec = mpc_get_prec(ctx->mfpc1[0]);

    ctx->approx_root = (mps_approximation**)mps_realloc(ctx->approx_root, sizeof(mps_approximation*) * n);
    for (i = ctx->n - ctx->zero_roots; i < n; i++)
    {
        ctx->approx_root[i] = mps_approximation_new(ctx);
    }

    ctx->order = (int*)mps_realloc(ctx->order, sizeof(int) * n);

    ctx->fppc1 = (cplx_t*)mps_realloc(ctx->fppc1, sizeof(cplx_t) * (n + 1));
    ctx->mfpc1 = (mpc_t*)mps_realloc(ctx->mfpc1, sizeof(mpc_t) * (n + 1));

    for (i = ctx->n + 1 - ctx->zero_roots; i < n + 1; i++)
        mpc_init2(ctx->mfpc1[i], previous_prec);

    ctx->mfppc1 = (mpc_t*)mps_realloc(ctx->mfppc1, sizeof(mpc_t) * (n + 1));
    for (i = ctx->n + 1 - ctx->zero_roots; i <= n; i++)
        mpc_init2(ctx->mfppc1[i], previous_prec);

    /* temporary vectors */
    ctx->spar1 = (mps_boolean*)mps_realloc(ctx->spar1, sizeof(mps_boolean) * (n + 2));
    ctx->again_old = (mps_boolean*)mps_realloc(ctx->again_old, sizeof(mps_boolean) * (n));

    ctx->fap1 = (double*)mps_realloc(ctx->fap1, sizeof(double) * (n + 1));
    ctx->fap2 = (double*)mps_realloc(ctx->fap2, sizeof(double) * (n + 1));

    ctx->dap1 = (rdpe_t*)mps_realloc(ctx->dap1, sizeof(rdpe_t) * (n + 1));
    ctx->dpc1 = (cdpe_t*)mps_realloc(ctx->dpc1, sizeof(cdpe_t) * (n + 1));
    ctx->dpc2 = (cdpe_t*)mps_realloc(ctx->dpc2, sizeof(cdpe_t) * (n + 1));

    /* Setting some default here, that were not settable because we didn't know
    * the degree of the polynomial */
    for (i = 0; i < n; i++)
        ctx->approx_root[i]->wp = DBL_DIG * (int)LOG2_10;
}

void
    mps_context_resize(mps_context* ctx, int n)
{
    /* We're excluding the case n == ctx->n that, clearly, doesn't
    * need any operation at all. */
    if (n > ctx->n)
        mps_context_expand(ctx, n);
    if (n < ctx->n)
        mps_context_shrink(ctx, n);
}

void
    mps_context_set_degree(mps_context* ctx, int n)
{
    if (ctx->initialized)
    {
        if (ctx->secular_equation != NULL)
        {
            mps_secular_equation_free(ctx, MPS_POLYNOMIAL(ctx->secular_equation));
            ctx->secular_equation = NULL;
        }

        mps_context_resize(ctx, n);
    }

    ctx->deg = ctx->n = n;

    /* Check if the numer of thread is greater of the number of roots,
    and in that case decrease it */
    if (ctx->n_threads > ctx->deg)
    {
        MPS_DEBUG_WITH_INFO(ctx, "Adjusting concurrency limit to %d", ctx->deg);
        mps_thread_pool_set_concurrency_limit(ctx, ctx->pool, ctx->deg);
    }

    /* If a secular equation is present in the old context we should free it
    * now so it will be reallocated on the first call to the algorithm. */
    if (ctx->secular_equation && MPS_POLYNOMIAL(ctx->secular_equation)->degree < n)
        mps_secular_equation_free(ctx, MPS_POLYNOMIAL(ctx->secular_equation));
    ctx->secular_equation = NULL;

}

/**
* @brief Set the monomial poly p as the input polynomial for
* the current equation.
*
* @param ctx The mps_context to set the monomial_poly into.
* @param p The mps_monomial_poly to solve.
*/
void
    mps_context_set_input_poly(mps_context* ctx, mps_polynomial* p)
{
    MPS_DEBUG_THIS_CALL(ctx);

    MPS_DEBUG(ctx, "Setting input poly");

    if (p->degree < 0)
    {
        mps_error(ctx, "Polynomial degree should be positive");
        return;
    }

    int i;
    ctx->active_poly = p;

    if (!p->thread_safe)
        mps_thread_pool_set_concurrency_limit(ctx, ctx->pool, 1);

    /* Set the density or sparsity of the polynomial, if it's not
    * a user polynomial */
    if (MPS_IS_MONOMIAL_POLY(p))
    {
        int original_degree = p->degree;
        mps_monomial_poly* mp = MPS_MONOMIAL_POLY(p);

        /* Deflate the polynomial if necessary */
        mps_monomial_poly_deflate(ctx, p);
        ctx->zero_roots = original_degree - p->degree;

        MPS_DEBUG_WITH_INFO(ctx, "Degree = %d", p->degree);

        /* Check if the input polynomial is sparse or not. We can simply check if
        * the again vector is all of true values */
        p->density = MPS_DENSITY_DENSE;
        for (i = 0; i <= MPS_POLYNOMIAL(mp)->degree; ++i)
        {
            if (!mp->spar[i])
            {
                p->density = MPS_DENSITY_SPARSE;
                break;
            }
        }
    }

    mps_context_set_degree(ctx, p->degree);
}

/**
* @brief Set active polynomial as a real floating point coefficient
* polynomial of degree <code>n</code> with coefficient exactly
* determined by components of vector coeff.
*
* Precisely, if \f${\mathrm coeff}\f$ is a vector of \f$n+1\f$ components,
* \f[
*   p(x) = \sum_{i = 0}^{n} {\mathrm coeff}_i  x^i
* \f]
*/
int
    mps_context_set_poly_d(mps_context* ctx, cplx_t* coeff, long unsigned int n)
{
    long unsigned int i;

    /* Allocate space for a polynomial of degree n */
    mps_monomial_poly* p = mps_monomial_poly_new(ctx, n);

    /* Fill polynomial coefficients */
    for (i = 0; i <= n; i++)
    {
        mps_monomial_poly_set_coefficient_d(ctx, p, i, cplx_Re(coeff[i]),
            cplx_Im(coeff[i]));
    }

    mps_context_set_input_poly(ctx, MPS_POLYNOMIAL(p));

    return 0;
}

/**
* @brief Set active polynomial as a integer coefficient
* polynomial of degree <code>n</code> with coefficient exactly
* determined by components of vector coeff.
*
* Precisely, if \f${\mathrm coeff}\f$ is a vector of \f$n+1\f$ components,
* \f[
*   p(x) = \sum_{i = 0}^{n} {\mathrm coeff}_i  x^i
* \f]
*/
int
    mps_context_set_poly_i(mps_context* ctx, int* coeff, long unsigned int n)
{
    long unsigned int i;

    /* Allocate data in mps_context to hold the polynomial of degree n */
    mps_monomial_poly* p = mps_monomial_poly_new(ctx, n);

    /* Fill polynomial */
    for (i = 0; i <= n; i++)
    {
        mpq_set_si(p->initial_mqp_r[i], coeff[i], 1U);
    }

    mps_context_set_input_poly(ctx, MPS_POLYNOMIAL(p));

    return 0;
}

/**
* @brief Set <code>roots[i]</code> to the i-th root of the polynomial
* and (if it is not <code>NULL</code>) <code>radius[i]</code>
* to the i-th inclusion radius.
*
* @param ctx The current mps_context.
* @param roots A pointer to an array of cplx_t variables. if *roots == NULL,
* MPSolve will take care of allocating these for you. You are in charge to free
* them when you don't need them anymore.
*
* @param radius A pointer to an array of double where MPSolve should store the
* inclusion radii. If *radius == NULL MPSolve will allocate those radii for you.
* If radius == NULL no radii will be returned.
*/
int
    mps_context_get_roots_d(mps_context* ctx, cplx_t** roots, double** radius)
{
    int i;

    if (*roots == NULL)
        *roots = cplx_valloc(ctx->n);

    if (radius && !*radius)
        *radius = double_valloc(ctx->n);

    for (i = 0; i < ctx->n; i++)
    {
        if (radius && *radius != NULL)
        {
            if (ctx->lastphase == float_phase || ctx->lastphase == dpe_phase)
            {
                (*radius)[i] = ctx->approx_root[i]->frad;
            }
            else
            {
                (*radius)[i] = rdpe_get_d(ctx->approx_root[i]->drad);
            }
        }

        if (ctx->lastphase == mp_phase)
        {
            mpc_get_cplx((*roots)[i], ctx->approx_root[i]->mvalue);
        }
        else if (ctx->lastphase == float_phase)
        {
            cplx_set((*roots)[i], ctx->approx_root[i]->fvalue);
        }
        else if (ctx->lastphase == dpe_phase)
        {
            cdpe_get_x((*roots)[i], ctx->approx_root[i]->dvalue);
        }
    }
    return 0;
}

/**
* @brief Get the roots computed as multiprecision complex numbers.
*
* @param roots A pointer to an array of mpc_t variables. if *roots == NULL,
* MPSolve will take care of allocate and init those for you. You are in charge to free
* and clear them when you don't need them anymore.
*
* @param radius A pointer to an array of rdpe_t where MPSolve should store the
* inclusion radii. If *radius == NULL MPSolve will allocate those radii for you.
* If radius == NULL no radii will be returned.
*/
int
    mps_context_get_roots_m(mps_context* ctx, mpc_t** roots, rdpe_t** radius)
{
    int i;

    if (!*roots)
    {
        *roots = mpc_valloc(ctx->n);
        mpc_vinit2(*roots, ctx->n, 0);
    }

    if (radius && !*radius)
    {
        *radius = rdpe_valloc(ctx->n);
    }

    {
        mpc_t* local_roots = *roots;
        rdpe_t* local_radius = radius ? *radius : NULL;

        for (i = 0; i < ctx->n; i++)
        {
            mpc_set_prec(local_roots[i], mpc_get_prec(ctx->approx_root[i]->mvalue));
            mpc_set(local_roots[i], ctx->approx_root[i]->mvalue);

            if (radius)
                rdpe_set(local_radius[i], ctx->approx_root[i]->drad);
        }
    }

    return 0;
}

/**
* @brief Retrieve a pointer to the current approximations in the context.
*
* @param ctx The current mps_context.
* @return A vector of n mps_approximation pointer,
* where n is the degree of the current polynomial
* that can be retrieve with mps_context_get_degree().
*
* Note that the value returned is a copy of the original approximations, so
* it should be freed when not needed anymore.
*
* Also note that mps_approximation_free() need the context where the approximations
* were created to proceed, so you should free those approximaitions _before_
* freeing the context.
*/
mps_approximation**
    mps_context_get_approximations(mps_context* ctx)
{
    mps_approximation** approximations = NULL;
    int i;

    if (!ctx->approx_root)
    {
        /* This means that an error occurred in the previous computation and the
        * polynomial has not been solved (or mps_mpsolve has not been called at all). */
        return NULL;
    }
    else
        approximations = mps_newv(mps_approximation*, ctx->n + ctx->zero_roots);

    for (i = 0; i < ctx->n; i++)
    {
        approximations[i] = mps_approximation_copy(ctx, ctx->approx_root[i]);

        /* Copy relevant data from the multiprecision values */
        mpc_get_cdpe(approximations[i]->dvalue, approximations[i]->mvalue);
        mpc_get_cplx(approximations[i]->fvalue, approximations[i]->mvalue);
        approximations[i]->frad = rdpe_get_d(approximations[i]->drad);
    }

    for (i = ctx->n; i < ctx->n + ctx->zero_roots; i++)
    {
        approximations[i] = mps_approximation_new(ctx);
        approximations[i]->status = MPS_ROOT_STATUS_APPROXIMATED;

        mpc_set_ui(approximations[i]->mvalue, 0U, 0U);
        cdpe_set(approximations[i]->dvalue, cdpe_zero);
        cplx_set(approximations[i]->fvalue, cplx_zero);

        rdpe_set(approximations[i]->drad, rdpe_zero);
        approximations[i]->frad = 0.0;
    }

    return approximations;
}


/**
* @brief Set the output precision for the roots.
*
* This has different meaning based on the output goal.
* If the goal is <code>MPS_OUTPUT_GOAL_ISOLATE</code>, this
* is the maximum precision used to try to isolate the roots,
* but roots won't be approximated at this precision if they
* are isolated with less precision.
*
* If the goal is <code>MPS_OUTPUT_GOAL_APPROXIMATE</code>,
* this is the minimum precision required for the roots in
* output.
*
* @param ctx The <code>mps_context</code> of the computation.
* @param prec The desired output precision.
*/
void
    mps_context_set_output_prec(mps_context* ctx, long int prec)
{
    ctx->output_config->prec = prec;
    rdpe_set_2dl(ctx->eps_out, 1.0, -prec);
}

/**
* @brief Set the bits of precision of the input coefficients.
*
* This has meaning only for fp coefficients, and the special
* value 0 means infinite precision.
*
* @param ctx The mps_context of the current computation.
* @param prec The precisision to be set.
*/
void
    mps_context_set_input_prec(mps_context* ctx, long int prec)
{
    if (!ctx->active_poly)
        return;
    else
        mps_polynomial_set_input_prec(ctx, ctx->active_poly, prec);
}

/**
* @brief Set the desired output format that will be used when
* calling mps_output().
*
* @param ctx The <code>mps_context</code> of the current computation.
* @param format The format chosen for the output.
*/
void
    mps_context_set_output_format(mps_context* ctx, mps_output_format format)
{
    ctx->output_config->format = format;

    /* Special handling of GNUPLOT format */
    if (format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
        ctx->gnuplot_format = "xyerrorbars";
    }
}

/**
* @brief Set the output goal for the computation.
*
* @param ctx The <code>mps_context</code> of the computation.
* @param goal The goal that will be reached before stopping.
*/
void
    mps_context_set_output_goal(mps_context* ctx, mps_output_goal goal)
{
    ctx->output_config->goal = goal;
}

/**
 * @brief Restrict the search set for the roots.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param set The search set for the roots.
 */
void
mps_context_set_search_set(mps_context *ctx, mps_search_set set)
{
  ctx->output_config->search_set = set;
}

/**
* @brief Set the value of the jacobi iterations switch in the MPSolve context.
*
* If jacobi_iterations is true then the Ehrlich-Aberth iterations will be carried
* out in a Jacobi fashion, otherwise Gauss-Seidel style will be employed.
*
* @param ctx The mps_context where the value will be set
* @param jacobi_iterations The desired value for the jacobi_iterations switch.
*/
void
    mps_context_set_jacobi_iterations(mps_context* ctx, mps_boolean jacobi_iterations)
{
    ctx->jacobi_iterations = jacobi_iterations;
}


/**
* @brief Set the debug level in MPSolve.
*
* @param ctx The <code>mps_context</code> of the current computation.
* @param level The debug_level to set.
*/
void
    mps_context_set_debug_level(mps_context* ctx, mps_debug_level level)
{
    ctx->debug_level = level;
    if (level)
    {
        ctx->DOLOG = true;
        if (!ctx->logstr)
            ctx->logstr = stderr;
    }
}

/**
* @brief Add another debug domain to the ones displayed. This will
* enable debug if disabled and show message from the given region.
*
* @param ctx The <code>mps_context</code> of the current computation.
* @param level The domain to add to the already set debug_level.
*/
void
    mps_context_add_debug_domain(mps_context* ctx, mps_debug_level level)
{
    mps_context_set_debug_level(ctx, ctx->debug_level | level);
}

/**
* @brief Get a pointer to the input config stored in the
* given mps_context.
*
* @param ctx The <code>mps_context</code> of the current computation.
*/
mps_input_configuration*
    mps_context_get_input_config(mps_context* ctx)
{
    return ctx->input_config;
}

/**
* @brief Get a pointer to the output config stored in the
* given mps_context.
*
* @param ctx The <code>mps_context</code> of the current computation.
*/
mps_output_configuration*
    mps_context_get_output_config(mps_context* ctx)
{
    return ctx->output_config;
}


/**
* @brief Set logstr as the default output for logging.
*
* @param ctx The <code>mps_context</code> of the current computation.
* @param logstr The desired stream to be used for logging.
*/
void
    mps_context_set_log_stream(mps_context* ctx, FILE* logstr)
{
    ctx->logstr = logstr;
}

/**
 * @brief Set root stream.
 *
 * @param s The <code>mps_context</code> of the current computation.
 * @param rtstr The stream used to resume an interrupted computation
 * or to load the approximations from a custom file.
 */
void
mps_context_set_root_stream(mps_context* ctx, FILE* rtstr)
{
    ctx->rtstr = rtstr;
}

/**
* @brief Set the phase from which the computation should start.
*
* @param ctx The <code>mps_context</code> of the current computation.
* @param phase The phase which should be chosen at the start of the
* computation.
*/
void
    mps_context_set_starting_phase(mps_context* ctx, mps_phase phase)
{
    ctx->input_config->starting_phase = phase;
}

/**
* @brief Get the number of zero roots in the output. Must be called
* after mps_mpsolve() has complete.
*
* @param ctx The <code>mps_context</code> of the current computation.
*/
int
    mps_context_get_zero_roots(mps_context* ctx)
{
    return ctx->zero_roots;
}

/**
* @brief Return true of the computation has passed the maximum
* admitted precision, and so was unable to reach desired output
* precision without any further information on the input
* coefficients.
*
* @param ctx The <code>mps_context</code> of the current computation.
*/
mps_boolean
    mps_context_get_over_max(mps_context* ctx)
{
    return ctx->over_max;
}

/**
* @brief Return true if mpsolve has encountered an error
* in the computation.
*
* @param ctx The <code>mps_context</code> of the current computation.
*/
mps_boolean mps_context_has_errors(mps_context* ctx)
{
    return ctx->error_state;
}

/**
* @brief Return a copy of the string describing the error msg, or
* NULL if no error has been encountered.
*
* This char array should be freed by the user.
*
* @param ctx The <code>mps_context</code> of the current computation.
*/
char* mps_context_error_msg(mps_context* ctx)
{
    if (ctx->last_error)
        return mps_strdup(ctx->last_error);
    else
        return NULL;
}

/**
* @brief Retrieve a pointer to the active polynomial being solved.
*
* @return A pointer to the requested active polynomial.
*/
mps_polynomial*
    mps_context_get_active_poly(mps_context* ctx)
{
    return ctx->active_poly;
}

/**
* @brief Retrieve the status of the root in position i.
*
* This method can be used to obtain more insight on the status of
* the approximations previously obtained by a call to mps_context_get_roots_m()
* or mps_context_get_roots_d().
*
* @return A copy to the mps_root_status of the approximation.
*/
mps_root_status
    mps_context_get_root_status(mps_context* ctx, int i)
{
    return ctx->approx_root[i]->status;
}

/**
* @brief Set the internal flag "avoid_multiprecision" to the specified value.
*
* If avoid_multiprecision is true MPSolve will not enter a multiprecision stage
* thus making impossible the computation of more digits than the one that are
* representable in standard floating point.
*
* This may be a useful flag to approximate roots of very ill-conditioned polynomials
* of high degree when no strict isolation is required. In this case it's possible to
* obtain very good approximations that are not Newton-isolated but are still satisfactory.
*/
void
    mps_context_set_avoid_multiprecision(mps_context* ctx, mps_boolean avoid_multiprecision)
{
    ctx->avoid_multiprecision = avoid_multiprecision;
}

/**
* @brief Enable the "crude" only approximation mode of MPSolve.
*
* If this mode is activated MPSolve will only perform a basic Aberth iteration
* in floating point and then exit. Note that the output result will still be
* guaranteed but in general it will not be possible to reach arbitrary precision
* and the results may be quite far from the roots for bad conditioned polynomials.
*/
void
    mps_context_set_crude_approximation_mode(mps_context* ctx, mps_boolean crude_approximation_mode)
{
    ctx->crude_approximation_mode = crude_approximation_mode;
}

/**
* @brief Select the regeneration driver implementation to use.
*
* The user of the library can change the default regeneration driver by providing
* a custom implementation of mps_regeneration_driver and calling this function.
*
* The custom implementation is obtained by hooking the function pointers defined
* in mps_regeneration_driver to custom routines that update a secular equation
* given new nodes.
*
* As usual, three routines for every data type supported by MPSolve must be provided.
*
* @param ctx The context where the change will have effect.
* @param rd The custom regeneration driver. Optionally this field may be NULL to restore
* the default implementation.
*/
void
    mps_context_set_regeneration_driver(mps_context* ctx, mps_regeneration_driver* rd)
{
    ctx->regeneration_driver = rd;
}

/**
 * @brief Set the number of threads.
 *
 * @param ctx The context where the change will have effect.
 * @param n_threads The number of thread to be spawned.
 */
void mps_context_set_n_threads(mps_context* ctx, int n_threads)
{
    ctx->n_threads = n_threads;
}

/**
* @brief Define which properties of the roots must be determined by MPSolve.
 *
 * @param ctx The context where the change will have effect.
 * @param root_properties Possible values:
 *  -# MPS_OUTPUT_PROPERTY_NONE
 *  -# MPS_OUTPUT_PROPERTY_REAL
 *  -# MPS_OUTPUT_PROPERTY_IMAGINARY
 *  -# MPS_OUTPUT_PROPERTY_REAL | MPS_OUTPUT_PROPERTY_IMAGINARY
 */
void mps_context_set_root_properties(mps_context* ctx, char root_properties)
{
    ctx->output_config->root_properties = root_properties;
}
