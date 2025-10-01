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
* @brief General routines ported from MPSolve 2.2
*/

#ifndef MPS_SOLVE_H_
#define MPS_SOLVE_H_

MPS_BEGIN_DECLS

	/* functions in solve.cpp */
	void mps_update(mps_context* ctx);
void mps_fsrad(mps_context* ctx, mps_cluster* cluster, cplx_t sc, double* sr);
void mps_dsrad(mps_context* ctx, mps_cluster* cluster, cdpe_t sc, rdpe_t sr);
void mps_msrad(mps_context* ctx, mps_cluster* cluster, mpc_t sc, rdpe_t sr);

mps_boolean mps_check_stop(mps_context* ctx);
void mps_fsolve(mps_context* ctx, mps_boolean* d_after_f);
void mps_dsolve(mps_context* ctx, mps_boolean d_after_f);
void mps_msolve(mps_context* ctx);
void mps_fpolzer(mps_context* ctx, int* it, mps_boolean* excep);
void mps_dpolzer(mps_context* ctx, int* it, mps_boolean* excep);
void mps_mpolzer(mps_context* ctx, int* it, mps_boolean* excep);

/* Functions in modify.cpp */
void mps_fmodify(mps_context* ctx, mps_boolean track_new_cluster);
void mps_dmodify(mps_context* ctx, mps_boolean track_new_cluster);
void mps_mmodify(mps_context* ctx, mps_boolean track_new_cluster);

/* Functions in inclusion.cpp */
void mps_fupdate_inclusions(mps_context* ctx);
void mps_dupdate_inclusions(mps_context* ctx);
void mps_mupdate_inclusions(mps_context* ctx);

/* Functions in test.cpp */
mps_boolean mps_inclusion(mps_context* ctx);

MPS_END_DECLS

#endif /* MPS_SOLVE_H_ */
