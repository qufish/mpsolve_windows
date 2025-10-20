/*
* This file is part of MPSolve 3.2.2
*
* Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
* License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
*
* Authors:
*   Leonardo Robol <leonardo.robol@unipi.it>
*/

#ifndef MPS_VERSION_H_
#define MPS_VERSION_H_

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "3"
#endif


#ifndef MPS_MAJOR_VERSION
#define MPS_MAJOR_VERSION 3
#endif

#ifndef MPS_MINOR_VERSION
#define MPS_MINOR_VERSION 2
#endif

#ifndef MPS_PATCH_VERSION
#define MPS_PATCH_VERSION 2
#endif

const char* mps_get_version();

unsigned int mps_get_major_version();

unsigned int mps_get_minor_version();

unsigned int mps_get_patch_version();


#endif /* MPS_VERSION_H_ */
