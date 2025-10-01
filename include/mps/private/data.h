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
*
* @brief Implementation of the routines that handle the management of data inside
* mps_context objects.
*/

#ifndef MPS_DATA_H_
#define MPS_DATA_H_

#include <mps/mps.h>

MPS_BEGIN_DECLS

/* functions in data.cpp */
void mps_mp_set_prec(mps_context* ctx, long int prec);
void mps_allocate_data(mps_context* ctx);
void mps_prepare_data(mps_context* ctx, long int prec);
void mps_restore_data(mps_context* ctx);
void mps_free_data(mps_context* ctx);
long int mps_raise_data(mps_context* ctx, long int prec);
void mps_raise_data_raw(mps_context* ctx, long int prec);

/* functions in main.cpp */
void mps_setup(mps_context* ctx);
void mps_check_data(mps_context* ctx, char* which_case);
void mps_compute_sep(mps_context* ctx);

MPS_END_DECLS

#endif /* endif _MPS_DATA_H */
