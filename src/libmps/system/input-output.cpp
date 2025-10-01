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

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <mps/mps.h>
#include <ctype.h>
#include <locale.h>

#define ISZERO -1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef __WINDOWS
#include <unistd.h>
#else
#include <io.h>
#endif

/**
* @brief Read the approximations from the file opened in the
* rtstr member of the mps_context.
*
* @param ctx A pointer to the current mps_context.
*
* NOTE: This function is not used anywhere in the code, atm.
*/
void
    mps_readroots(mps_context* ctx)
{
    long digits;
    int i, read_elements;

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "Reading roots...\n");

    read_elements = fscanf(ctx->rtstr, "%ld", &digits);
    if (!read_elements)
    {
        mps_error(ctx, "Error while reading roots, aborting.");
    }

    /* precision setup code goes here */

    for (i = 0; i < ctx->n; i++)
        mpc_inp_str_u(ctx->approx_root[i]->mvalue, ctx->rtstr, 10);
}

/**
* @brief Count the roots that are included in the search set, excluded
* form it, or have an undetermined inclusion state.
*
* @param ctx A pointer to the current mps_context.
*/
void
    mps_countroots(mps_context* ctx)
{
    int k;

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "Counting roots...\n");

    ctx->count[0] = ctx->count[1] = ctx->count[2] = 0;

    for (k = 0; k < ctx->n; k++)
        switch (ctx->approx_root[k]->inclusion)
        {
        case MPS_ROOT_INCLUSION_IN:
            ctx->count[0]++;
            break;

        case MPS_ROOT_INCLUSION_OUT:
            ctx->count[1]++;
            break;

        default:
            ctx->count[2]++;
            break;
        }

    if (ctx->output_config->search_set == MPS_SEARCH_SET_UNITARY_DISC_COMPL)
        ctx->count[1] += ctx->zero_roots;
    else
        ctx->count[0] += ctx->zero_roots;
}

/**
* @brief Print a summary of the count of the roots that can be obtained
* through a call to mps_countroots() to the outstr member of the
* current mps_context.
*
* @param ctx A pointer to the current mps_context.
*/
void
    mps_outcount(mps_context* ctx)
{
    mps_countroots(ctx);

    fprintf(ctx->outstr, "%d roots are inside;\n", ctx->count[0]);
    fprintf(ctx->outstr, "%d roots are outside;\n", ctx->count[1]);
    fprintf(ctx->outstr, "%d roots are uncertain.\n", ctx->count[2]);
    if (ctx->DOLOG)
    {
        fprintf(ctx->logstr, "%d roots are inside;\n", ctx->count[0]);
        fprintf(ctx->logstr, "%d roots are outside;\n", ctx->count[1]);
        fprintf(ctx->logstr, "%d roots are uncertain.\n", ctx->count[2]);
    }
}

/**
* @brief Print a float to stdout (or whatever the output stream is
* atm) respecting the given options, and only with the significant
* digits.
*
* @param ctx A pointer to the current mps_context.
* @param f The float approximation that should be printed.
* @param rad The current inclusion radius for that approximation.
* @param out_digit The number of output digits required.
* @param sign The sign of the approximation.
*/
MPS_PRIVATE void
    mps_outfloat(mps_context* ctx, mpf_t f, rdpe_t rad, long out_digit,
        mps_boolean sign)
{
    mpf_t t;
    rdpe_t r, ro;
    double d;
    long l, digit, true_digit;

    if (ctx->output_config->format == MPS_OUTPUT_FORMAT_FULL)
    {
        mpf_init2(t, mpf_get_prec(f));
        mpf_set(t, f);
        mpf_out_str(ctx->outstr, 10, 0, t);
        mpf_clear(t);
        return;
    }

    mpf_init2(t, ctx->output_config->prec);

    mpf_get_rdpe(ro, f);
    if (ctx->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT ||
        ctx->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
        rdpe_out_str_u(ctx->outstr, ro);
    else
    {
        rdpe_abs_eq(ro);
        if (rdpe_ne(ro, rdpe_zero))
            rdpe_div(r, rad, ro);
        else
            rdpe_set_d(r, 1.0e-10);
        digit = (long)(-rdpe_log10(r) + 1.5);
        if (digit <= 0)
        {
            rdpe_get_dl(&d, &l, ro);
            fprintf(ctx->outstr, "0.e%ld", l);
        }
        else
        {
            true_digit = (long)(LOG10_2 * mpf_get_prec(f)) + 1;
            true_digit = MIN(digit, true_digit);
            true_digit = MIN(true_digit, out_digit);
            if (sign)
                mpf_set(t, f);
            else
                mpf_abs(t, f);
            mpf_out_str(ctx->outstr, 10, true_digit, t);
        }
    }

    mpf_clear(t);
}

/**
* @brief Print an approximation to stdout (or whatever the output
* stream currently selected in the mps_context is).
*
* @param ctx A pointer to the current mps_context.
* @param i The index of the approxiomation that shall be printed.
* @param num The number of zero roots.
*/
MPS_PRIVATE void
mps_outroot(mps_context* ctx, int i, int num)
{
    long out_digit;

    out_digit = (long)(LOG10_2 * ctx->output_config->prec) + 10;

    /* print format header */
    switch (ctx->output_config->format)
    {
    case MPS_OUTPUT_FORMAT_COMPACT:
    case MPS_OUTPUT_FORMAT_FULL:
        fprintf(ctx->outstr, "(");
        break;

    case MPS_OUTPUT_FORMAT_VERBOSE:
        fprintf(ctx->outstr, "Root(%d) = ", num);
        break;

    default:
        break;
    }

    /* print real part */
    if (i == ISZERO || ctx->approx_root[i]->attrs == MPS_ROOT_ATTRS_IMAG)
        fprintf(ctx->outstr, "0");
    else
        mps_outfloat(ctx, mpc_Re(ctx->approx_root[i]->mvalue), ctx->approx_root[i]->drad, out_digit, true);

    /* print format middle part */
    switch (ctx->output_config->format)
    {
    case MPS_OUTPUT_FORMAT_BARE:
        fprintf(ctx->outstr, " ");
        break;

    case MPS_OUTPUT_FORMAT_GNUPLOT:
    case MPS_OUTPUT_FORMAT_GNUPLOT_FULL:
        fprintf(ctx->outstr, "\t");
        break;

    case MPS_OUTPUT_FORMAT_COMPACT:
    case MPS_OUTPUT_FORMAT_FULL:
        fprintf(ctx->outstr, ", ");
        break;

    case MPS_OUTPUT_FORMAT_VERBOSE:
        if (i == ISZERO || mpf_sgn(mpc_Im(ctx->approx_root[i]->mvalue)) >= 0)
            fprintf(ctx->outstr, " + I * ");
        else
            fprintf(ctx->outstr, " - I * ");
        break;

    default:
        break;
    }

    /* print imaginary part */
    if (i == ISZERO || ctx->approx_root[i]->attrs == MPS_ROOT_ATTRS_REAL)
        fprintf(ctx->outstr, "0");
    else
        mps_outfloat(ctx, mpc_Im(ctx->approx_root[i]->mvalue), ctx->approx_root[i]->drad, out_digit,
            ctx->output_config->format != MPS_OUTPUT_FORMAT_VERBOSE);

    /* If the output format is GNUPLOT_FORMAT_FULL, print out also the radius */
    if (ctx->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
        fprintf(ctx->outstr, "\t");
        rdpe_out_str_u(ctx->outstr, ctx->approx_root[i]->drad);
        fprintf(ctx->outstr, "\t");
        rdpe_out_str_u(ctx->outstr, ctx->approx_root[i]->drad);
    }

    /* print format ending */
    switch (ctx->output_config->format)
    {
    case MPS_OUTPUT_FORMAT_COMPACT:
        fprintf(ctx->outstr, ")");
        break;

    case MPS_OUTPUT_FORMAT_FULL:
        fprintf(ctx->outstr, ")\n");
        if (i != ISZERO)
        {
            rdpe_outln_str(ctx->outstr, ctx->approx_root[i]->drad);
            fprintf(ctx->outstr, "Status: %s, %s, %s\n",
                MPS_ROOT_STATUS_TO_STRING(ctx->approx_root[i]->status),
                MPS_ROOT_ATTRS_TO_STRING(ctx->approx_root[i]->attrs),
                MPS_ROOT_INCLUSION_TO_STRING(ctx->approx_root[i]->inclusion));
        }
        else
            fprintf(ctx->outstr, " 0\n ---\n");
        break;

    default:
        break;
    }
    fprintf(ctx->outstr, "\n");

    /* debug info */
    if (ctx->DOLOG)
    {
        if (i == ISZERO)
            fprintf(ctx->logstr, "zero root %-4d = 0", num);
        else
        {
            fprintf(ctx->logstr, "Root %-4d = ", i);
            mpc_out_str_2(ctx->logstr, 10, 0, 0, ctx->approx_root[i]->mvalue);
            fprintf(ctx->logstr, "\n");
            fprintf(ctx->logstr, "  Radius = ");
            rdpe_outln_str(ctx->logstr, ctx->approx_root[i]->drad);
            fprintf(ctx->logstr, "  Prec = %ld\n",
                (long)(mpc_get_prec(ctx->approx_root[i]->mvalue) / LOG2_10));
            fprintf(ctx->logstr, "  Approximation = %s\n",
                MPS_ROOT_STATUS_TO_STRING(ctx->approx_root[i]->status));
            fprintf(ctx->logstr, "  Attributes = %s\n",
                MPS_ROOT_ATTRS_TO_STRING(ctx->approx_root[i]->attrs));
            fprintf(ctx->logstr, "  Inclusion = %s\n",
                MPS_ROOT_INCLUSION_TO_STRING(ctx->approx_root[i]->inclusion));
            fprintf(ctx->logstr, "--------------------\n");
        }
    }
}

/**
* @brief Print the approximations to stdout (or whatever the output
* stream currently selected in the mps_context is).
*
* @param ctx A pointer to the current mps_context.
*/
void
    mps_output(mps_context* ctx)
{
    int i, ind, num = 0;

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "--------------------\n");

    if (ctx->output_config->format != MPS_OUTPUT_FORMAT_GNUPLOT &&
        ctx->output_config->format != MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
        if (ctx->over_max)
        {
            mps_warn(ctx, "Warning: Input precision has been reached during computation, "
                "so not all the required digits may have been computed.");
        }
    }

    /* Start with plotting instructions in the case of
    * MPS_OUTPUT_GNUPLOT_FULL, so the output can be
    * piped directly to gnuplot */
    if (ctx->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      fprintf (ctx->outstr, "# MPSolve output for GNUPLOT\n");
      fprintf (ctx->outstr, "# Make user that this output is piped into gnuplot using a command like\n");
      fprintf (ctx->outstr, "# mpsolve -Ogf | gnuplot \n");
      fprintf (ctx->outstr, "set pointsize 0.3\n");
      fprintf (ctx->outstr, "plot '-' title 'Computed roots' with %s\n", ctx->gnuplot_format);
    }

    if (ctx->output_config->goal == MPS_OUTPUT_GOAL_COUNT)
        mps_outcount(ctx);
    else
    {
        if (ctx->output_config->search_set != MPS_SEARCH_SET_UNITARY_DISC_COMPL)
            for (i = 0; i < ctx->zero_roots; i++)
                mps_outroot(ctx, ISZERO, num++);
        for (ind = 0; ind < ctx->n; ind++)
        {
            i = ctx->order[ind];
            if (ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_OUT)
                continue;
            mps_outroot(ctx, i, num++);
        }
    }

    if (ctx->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      fprintf (ctx->outstr, "e\n");
      fprintf (ctx->outstr, "pause mouse close\n");
      fprintf (ctx->outstr, "# End of MPSolve GNUPLOT output. If you are seeing this maybe\n");
      fprintf (ctx->outstr, "# you forgot to pipe the ***solve command into gnuplot?\n");
    }
}

/**
* @brief Update the MP version of the roots to the latest and greatest approximations.
*
* @param ctx A pointer to the current mps_context.
*/
MPS_PRIVATE void
    mps_copy_roots(mps_context* ctx)
{
    int i;

    MPS_DEBUG_THIS_CALL(ctx);

    switch (ctx->lastphase)
    {
    case no_phase:
        mps_error(ctx, "Nothing to copy");
        break;

    case float_phase:
        if (ctx->DOSORT)
            mps_fsort(ctx);
        for (i = 0; i < ctx->n; i++)
        {
            mpc_set_prec(ctx->approx_root[i]->mvalue, DBL_MANT_DIG);
            mpc_set_cplx(ctx->approx_root[i]->mvalue, ctx->approx_root[i]->fvalue);
            rdpe_set_d(ctx->approx_root[i]->drad, ctx->approx_root[i]->frad);
        }
        break;

    case dpe_phase:
        if (ctx->DOSORT)
            mps_dsort(ctx);
        for (i = 0; i < ctx->n; i++)
        {
            mpc_set_prec(ctx->approx_root[i]->mvalue, DBL_MANT_DIG);
            mpc_set_cdpe(ctx->approx_root[i]->mvalue, ctx->approx_root[i]->dvalue);
        }
        break;

    case mp_phase:
        if (ctx->DOSORT)
            mps_msort(ctx);
        break;
    }
}

/**
* @brief Dump all the current approximation to the logstr selected
* in the current mps_context.
*
* @param ctx A pointer to the current mps_context.
*
* This function is tipically used when encountering some errors.
*/
void
    mps_dump(mps_context* ctx)
{
    int i;
    FILE* dmpstr = ctx->logstr;

    MPS_DEBUG(ctx, "Dumping the approximations:");

    /* output current status */
    /* fprintf (dmpstr, */
    /*          "Phase=%d, In=%d, Out=%d, Uncertain=%d, Zero=%d, Clusters=%ld\n", */
    /*          ctx->lastphase, ctx->count[0], ctx->count[1], ctx->count[2], ctx->zero_roots, */
    /*          ctx->clusterization->n); */

    MPS_DEBUG(ctx, "Phase = %s, In = %d, Out = %d, Uncertain = %d, Zero = %d, Cluster = %ld",
        MPS_PHASE_TO_STRING(ctx->lastphase), ctx->count[0], ctx->count[1], ctx->count[2],
        ctx->zero_roots, ctx->clusterization->n);

    /* output current approximations */
    /* fprintf (dmpstr, "\nCurrent approximations:\n"); */
    MPS_DEBUG(ctx, "Current approximations:");
    for (i = 0; i < ctx->n; i++)
    {
        switch (ctx->lastphase)
        {
        case no_phase:
        case float_phase:
            MPS_DEBUG_CPLX(ctx, ctx->approx_root[i]->fvalue, "Approximation  %4d", i);
            break;

        case dpe_phase:
            MPS_DEBUG_CDPE(ctx, ctx->approx_root[i]->dvalue, "Approximation  %4d", i);
            break;

        case mp_phase:
            MPS_DEBUG_MPC(ctx, ctx->mpwp, ctx->approx_root[i]->mvalue, "Approximation  %4d", i);
            break;
        }
    }

    /* output radii */
    MPS_DEBUG(ctx, "Current radii:");
    for (i = 0; i < ctx->n; i++)
    {
        switch (ctx->lastphase)
        {
        case no_phase:
        case float_phase:
            MPS_DEBUG(ctx, "Radius of root %4d = %e", i, ctx->approx_root[i]->frad);
            break;

        case dpe_phase:
        case mp_phase:
            MPS_DEBUG_RDPE(ctx, ctx->approx_root[i]->drad, "Radius of root %4d", i);
            break;
        }
    }

    MPS_DEBUG(ctx, " ");
    mps_dump_status(ctx, dmpstr);
}

/**
* @brief Dump cluster structure to <code>outstr</code>.
*
* @param ctx the mps_context struct pointer.
* @param outstr The output stream where the cluster structure
*  will be dumped.
*/
MPS_PRIVATE void
    mps_dump_cluster_structure(mps_context* ctx, FILE* outstr)
{
  fprintf (outstr,
           "    MPS_DUMP_CLUSTER_STRUCTURE: Dumping cluster structure\n");

  mps_cluster_item * cluster_item;
  mps_cluster * cluster;
  mps_root * root;

  for (cluster_item = ctx->clusterization->first_cluster_item; cluster_item != NULL;
       cluster_item = cluster_item->next_cluster_item)
    {
      cluster = cluster_item->cluster;
      fprintf (outstr, "     Cluster contains %ld roots:\n", cluster->n);

      /* Dump cluster roots, but not more than 15 for line, to make
       * the output readable. */
      int j = 0;
      for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
          /* Go to a newlint if 15 roots are printed out */
          if (j % 15 == 0)
            {
              fprintf (outstr, "\n       ");
            }

          fprintf (outstr, " %4ld", root->k);
          j++;
        }

      /* Make space untile the next cluster */
      fprintf (outstr, "\n\n");
    }
}

/**
* @brief Dump status of all the root approximations
*/
MPS_PRIVATE void
    mps_dump_status(mps_context* ctx, FILE* outstr)
{
    int i;

    MPS_DEBUG(ctx, "              Approximation              Attributes       Inclusion");
    for (i = 0; i < ctx->n; i++)
    {
        MPS_DEBUG(ctx, "Status  %4d: %-25s  %-15s  %-15s", i,
            MPS_ROOT_STATUS_TO_STRING(ctx->approx_root[i]->status),
            MPS_ROOT_ATTRS_TO_STRING(ctx->approx_root[i]->attrs),
            MPS_ROOT_INCLUSION_TO_STRING(ctx->approx_root[i]->inclusion));
    }
}

/**
* @brief Print a warning to the user.
*
* @param st The current mps_context
* @param format The printf-like format for the data to print
*/
void
    mps_warn(mps_context* st, char* format, ...)
{
    char* exclamation_mark = "";

    va_list ap;

    va_start(ap, format);

    if (mps_is_a_tty(st->logstr))
        exclamation_mark = "\033[33;1m!\033[0m";

    if (st->DOWARN)
    {
        if (format[strlen(format)] == '\n')
        {
            fprintf(stderr, "%s ", exclamation_mark);
            vfprintf(stderr, format, ap);
        }
        else
        {
            fprintf(stderr, "%s ", exclamation_mark);
            vfprintf(stderr, format, ap);
            fprintf(stderr, "\n");
        }
    }

    va_end(ap);
}

/**
* @brief Check if the file descriptor associated to stream
* is bounded to a tty.
*
* @param stream the stream to check
*/
MPS_PRIVATE mps_boolean
    mps_is_a_tty(FILE* stream)
{
#ifndef __WINDOWS
    return isatty(fileno(stream));
#else
    return _isatty(_fileno(stream));
#endif
}

/* Prepare a implicit definition if not provided by the compiler but
* available as a non-conformant extension. */
#ifndef vsnprintf
#ifdef HAVE_VSNPRINTF
int snprintf(char* str, size_t size, const char* format, ...);
#endif
#endif

/**
* @brief Record an error happened during the computation.
* This will set the internal status error to on, and the actual
* errors will printed with the first call to mps_print_errors().
* @param ctx A pointer to the current mps_context.
*/
void
    mps_error(mps_context* ctx, const char* format, ...)
{
    va_list ap;
    int buffer_size = 32;
    int missing_characters = 0;

    va_start(ap, format);

    ctx->error_state = true;
    if (ctx->last_error == NULL)
        ctx->last_error = mps_newv(char, buffer_size);

    /* Measure space needed for the string, if our initial guess for the space neede
    * is not enough */
    while ((missing_characters = vsnprintf(ctx->last_error, buffer_size, format, ap)) > buffer_size)
    {
        buffer_size += missing_characters + 1;
        ctx->last_error = (char*)mps_realloc(ctx->last_error, buffer_size);
    }

    va_end(ap);
}

/**
* @brief Print all the errors that have been recorded up to now.
* This function should be called only if mps_context_has_errors()
* returns true.
*
* @param ctx A pointer to the current mps_context.
*/
void
    mps_print_errors(mps_context* ctx)
{
  const char *error = ctx->last_error;
  size_t length = (int)strlen (error);

  if (ctx->logstr == NULL)
    ctx->logstr = stderr;

  const char *exclamation_mark = "!";
  if (mps_is_a_tty (ctx->logstr))
    exclamation_mark = "\033[31;1m!\033[0m";

  if (error[length] == '\n')
    {
      fprintf (stderr, "%s %s %s", exclamation_mark, "MPSolve encountered an error:", error);
    }
  else
    {
      fprintf (stderr, "%s %s %s\n", exclamation_mark, "MPSolve encountered an error:", error);
    }

  /* Dump approximations, but only if they are present */
  if (ctx->approx_root && ctx->lastphase)
    mps_dump (ctx);
}
