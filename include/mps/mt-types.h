/*
* This file is part of MPSolve 3.2.2
*
* Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
* License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
*
* Authors:
*   Leonardo Robol <leonardo.robol@unipi.it>
*/

#include <mps/mps.h>

/**
* @file
* @brief Implementation of some thread-safe types that can be easily used
*/

#ifndef MPS_MT_TYPES_
#define MPS_MT_TYPES_

/**
 * @brief A thread safe version of mps_boolean.
 *
* Must be accessed locking/unlocking a mutex.
 */
struct mps_boolean_mt {
	mps_boolean value;
	mps_mutex_t mutex;
};


/**
 * @brief A thread safe version of mps_boolean.
 *
 * Must be accessed using the macro MPS_LOCK (x) and
 * MPS_UNLOCK (x).
 */
struct mps_long_int_mt {
	long int value;
	mps_mutex_t value_mutex;
};

#endif
