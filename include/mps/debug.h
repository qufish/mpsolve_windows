/*
* This file is part of MPSolve 3.2.2
*
* Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
* License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
*
* Authors:
*   Leonardo Robol <leonardo.robol@unipi.it>
*/

/**
* @file
* @brief Debugging functions, that honor the status os <code>ctx->DOLOG</code> and
* autodetect if the output stream is or not a tty to select the best output method.
*/

#ifndef  MPS_DEBUG_H_
#define  MPS_DEBUG_H_

#ifndef __WINDOWS
#include <mps_unistd.h>
#else
#include <io.h>
#endif

#include <mps/gmptools.h>
#include <gmp.h>
#include <time.h>

MPS_BEGIN_DECLS

    /* Timer functions */
    clock_t* mps_start_timer(void);
unsigned long int mps_stop_timer(clock_t* my_timer);

/* Keep away assert() when compiling without debug */
#ifdef DISABLE_DEBUG
#define NDEBUG
#endif

#ifdef NICE_DEBUG

/**
* @brief Shorthand for debugging with the MPS_DEBUG_INFO level.
*/
#define MPS_DEBUG_WITH_INFO(ctx, templ ...) {      \
if (ctx->debug_level & MPS_DEBUG_INFO) {      \
    MPS_DEBUG (ctx, templ);                      \
    }                                           \
}

/**
* @brief Shorthand for debugging with the MPS_DEBUG_IO level.
*/
#define MPS_DEBUG_WITH_IO(ctx, templ ...) {      \
if (ctx->debug_level & MPS_DEBUG_IO) {      \
    MPS_DEBUG (ctx, templ);                      \
    }                                           \
}

#define MPS_DEBUG(ctx, templ ...) {                \
__MPS_DEBUG (ctx, templ);                     \
if (ctx->DOLOG) {                             \
    fprintf (ctx->logstr, "\n");                 \
    }                                           \
}

/**
* @brief Debug the value of a complex multiprecision
* variable.
*/
#define MPS_DEBUG_MPC(ctx, digits, c, name ...) {  \
__MPS_DEBUG_EQ (ctx, name);                    \
if (ctx->DOLOG) {                             \
    mpc_outln_str (ctx->logstr, 10, digits, c);  \
    }                                           \
}

/**
* @brief Debug the value of a real multiprecision
* variable.
*/
#define MPS_DEBUG_MPF(ctx, digits, c, name ...) {  \
__MPS_DEBUG_EQ (ctx, name);                    \
if (ctx->DOLOG) {                             \
    mpf_out_str (ctx->logstr, 10, digits, c);    \
    fprintf (ctx->logstr, "\n");                \
    }                                           \
}

#define MPS_DEBUG_MPC2(ctx, radius, c, name ...) {         \
__MPS_DEBUG_EQ (ctx, name);                              \
if (ctx->DOLOG) {                                       \
    int t = -rdpe_log10 (radius) - 1;                     \
    cdpe_t ctmp;                                          \
    rdpe_t mm;                                            \
    mpc_get_cdpe (ctmp, c);                               \
    cdpe_mod (mm, ctmp);                                  \
    t += rdpe_log (mm);                                   \
    mpc_outln_str (ctx->logstr, 10, t, c);                  \
    }                                                     \
}


/**
    * @brief Debug the value of a rdpe variable.
    */
#define MPS_DEBUG_RDPE(ctx, r, name ...) {         \
__MPS_DEBUG_EQ (ctx, name);                    \
if (ctx->DOLOG) {                             \
    rdpe_outln_str (ctx->logstr, r);             \
    }                                           \
}

    /**
    * @brief Debug the value of a cdpe variable.
    */
#define MPS_DEBUG_CDPE(ctx, c, name ...) {         \
__MPS_DEBUG_EQ (ctx, name);                    \
if (ctx->DOLOG) {                             \
    cdpe_outln_str (ctx->logstr, c);             \
    }                                           \
}

    /**
    * @brief Debug the values of a cplx_t variable
    */
#define MPS_DEBUG_CPLX(ctx, c, name ...) {         \
__MPS_DEBUG_EQ (ctx, name);                    \
if (ctx->DOLOG) {                             \
    cplx_outln_str (ctx->logstr, c);             \
    }                                           \
}

    /**
    * @brief Make some space in the debug stream to make clean that
    * another section is starting.
    */
#define MPS_DEBUG_BREAK(ctx) if (ctx->DOLOG) {      \
    fprintf (ctx->logstr, "\n");                   \
}

    /**
        * @brief Low-level debug print that appends an " = " sign to the
        * output (useful for debugging values of variables).
        */
#define __MPS_DEBUG_EQ(ctx, templ ...)             \
__MPS_DEBUG (ctx, templ);                        \
if (ctx->DOLOG) {                               \
    fprintf (ctx->logstr, " = ");                  \
}

        /**
        * @brief Debug that a function is going to be called.
        */
#ifndef __WINDOWS
#define MPS_DEBUG_CALL(ctx, function) if (ctx->DOLOG && (ctx->debug_level & MPS_DEBUG_FUNCTION_CALLS)) { \
    if (isatty (fileno (ctx->logstr))) {                                   \
        __MPS_DEBUG (ctx, "Called \033[31;1m");                              \
    }                                                                   \
    else {                                                              \
        __MPS_DEBUG (ctx, "Called ");                                        \
    }                                                                   \
    if (isatty (fileno (ctx->logstr))) {                                   \
        fprintf (ctx->logstr, function); fprintf (ctx->logstr, "()\033[0m\n");  \
    }                                                                   \
    else                                                                \
    {                                                                 \
        fprintf (ctx->logstr, function); fprintf (ctx->logstr, "()\n");       \
    }                                                                 \
}
#else
#define MPS_DEBUG_CALL(ctx, function) if (ctx->DOLOG && (ctx->debug_level & MPS_DEBUG_FUNCTION_CALLS)) { \
    if (_isatty (_fileno (ctx->logstr))) {                                  \
        __MPS_DEBUG (ctx, "Called \033[31;1m");                              \
    }                                                                   \
    else {                                                              \
        __MPS_DEBUG (ctx, "Called ");                                        \
    }                                                                   \
    if (_isatty (_fileno (ctx->logstr))) {                                  \
        fprintf (ctx->logstr, function); fprintf (ctx->logstr, "()\033[0m\n");  \
    }                                                                   \
    else                                                                \
    {                                                                 \
        fprintf (ctx->logstr, function); fprintf (ctx->logstr, "()\n");       \
    }                                                                 \
}
#endif

#define MPS_DEBUG_THIS_CALL(ctx) MPS_DEBUG_CALL (ctx, __FUNCTION__)


        /**
        * @brief Low-level DEBUG() used by other MPS_DEBUG_* statements.
        */
#if __STDC_VERSION__ >= 199901L
#ifndef __WINDOWS
#define __MPS_DEBUG(ctx, templ ...) if (ctx->DOLOG) {                \
    if (isatty (fileno (ctx->logstr))) {                           \
        fprintf (ctx->logstr, "%s:%d \033[32;1m%s()\033[;0m ",       \
                __FILE__, __LINE__, __FUNCTION__);                \
    }                                                           \
    else {                                                      \
        fprintf (ctx->logstr, "%s:%d %s() ",                         \
                __FILE__, __LINE__, __FUNCTION__);                \
    }                                                           \
    fprintf (ctx->logstr, templ);                                  \
}
#else
#define __MPS_DEBUG(ctx, templ ...) if (ctx->DOLOG) {                \
    if (_isatty (_fileno (ctx->logstr))) {                          \
        fprintf (ctx->logstr, "%s:%d \033[32;1m%s()\033[;0m ",       \
                __FILE__, __LINE__, __FUNCTION__);                \
    }                                                           \
    else {                                                      \
        fprintf (ctx->logstr, "%s:%d %s() ",                         \
                __FILE__, __LINE__, __FUNCTION__);                \
    }                                                           \
    fprintf (ctx->logstr, templ);                                  \
}
#endif

#endif
#else
#define MPS_DEBUG(args, ...)
#define MPS_DEBUG_MPC(args, ...)
#define MPS_DEBUG_CPLX(args, ...)
#define MPS_DEBUG_RDPE(args, ...)
#define MPS_DEBUG_CDPE(args, ...)
#define MPS_DEBUG_CALL(args, ...)
#define MPS_DEBUG_MCLUSTER_ROOTS(args, ...)
#define MPS_DEBUG_THIS_CALL(ctx)
#define MPS_DEBUG_WITH_INFO(args, ...)
#define MPS_DEBUG_WITH_IO(args, ...)
#define __MPS_DEBUG(args, ...)
#endif

/*
* Ansi C version implemented without using variadic macros, only for compatibility
* since could became imprecise when using complex debug instructions.
*/
#ifndef DISABLE_DEBUG
#if __STDC_VERSION__ < 199901L
#include <mps/mps.h>
#ifdef MPS_DEBUG
#undef MPS_DEBUG
#endif
#ifdef __MPS_DEBUG
#undef __MPS_DEBUG
#endif
#define MPS_DEBUG __c_impl__MPS_DEBUG
#define __MPS_DEBUG __c_impl____MPS_DEBUG
void __c_impl__MPS_DEBUG(mps_context* ctx, const char* templ, ...);
void __c_impl____MPS_DEBUG(mps_context* ctx, const char* templ, ...);
#endif
#endif

/**
* @brief This is the flag that enables debugging of general
* information about the flow of the program.
*/
#define MPS_DEBUG_INFO     (0x0001 << 0)

/**
* @brief This is the flag used to debug cluster-related
* informations.
*/
#define MPS_DEBUG_CLUSTER  (0x0001 << 1)

/**
    * @brief This is the flag used to debug informations on
    * approximations during the execution of MPSolve.
    */
#define MPS_DEBUG_APPROXIMATIONS (0x0001 << 2)

    /**
    * @brief Flag used to show debug about final approximation
    * of the roots using newton iterations.
    */
#define MPS_DEBUG_IMPROVEMENT (0x0001 << 3)

    /**
    * @brief Flag used to obtain information about timings in
    * the algorithm.
    */
#define MPS_DEBUG_TIMINGS (0x0001 << 4)

    /**
    * @brief Debug function calls
    */
#define MPS_DEBUG_FUNCTION_CALLS (0x0001 << 5)

    /**
        * @brief Debug I/O informations
        */
#define MPS_DEBUG_IO (0x0001 << 6)

        /**
        * @brief Debug memory management
        */
#define MPS_DEBUG_MEMORY (0x0001 << 7)

        /**
        * @brief Debug checks for the convergence
        * of the various iterations packets.
        */
#define MPS_DEBUG_PACKETS (0x0001 << 8)

        /**
        * @brief Debug the regenration of the coefficients-
        */
#define MPS_DEBUG_REGENERATION (0x0001 << 9)

        /**
            * @brief This is the flag used to debug informations
            * about virtually every step in the program, it enables
            * every debug level.
            */
#define MPS_DEBUG_TRACE    (0xFFFF)

#define MPS_DEBUG_IF(ctx, debug_level, debug_instruction) if (debug_level & ctx->debug_level) { \
    debug_instruction;                                                  \
}

MPS_END_DECLS

#endif                          /* DEBUG_H */
