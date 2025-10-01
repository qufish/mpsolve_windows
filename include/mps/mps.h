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
* This file is the header for the libmps library. Including
* this file is needed to access all the MPSolve routines by
* MPSolve internals.
*
* @brief Header file for libmps
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef MPS_CORE_H_
#define MPS_CORE_H_

#define __MPS_NOT_DEFINE_BOOL

#ifdef __MPS_MATLAB_MODE
#define __MPS_NOT_DEFINE_BOOL
#endif

#ifdef HAVE_HIDDEN_VISIBILITY_ATTRIBUTE
  #ifdef MPS_PUBLISH_PRIVATE_METHODS
    #define MPS_PRIVATE
  #else
    #define MPS_PRIVATE __attribute__((visibility ("hidden")))
  #endif
#else
  #define MPS_PRIVATE
#endif
#define _MPS_PRIVATE

#define MPS_BEGIN_DECLS
#define MPS_END_DECLS

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//#define MPS_USE_PTHREADS // if toggling back to pthreads, will need pthrerads.h and several lib files to be made available
#ifndef MPS_USE_PTHREADS

//  implement threading and mutexes using the Windows std lib methods
#include <condition_variable>
#include <mutex>
#include <thread>

#define mps_thread_join(threadptr) \
    (*(threadptr)).join(); \
    mps_del_obj(threadptr)
#define mps_thread_processing_exit(thread) 
#define mps_thread_t std::thread
#define mps_start_new_thread(threadptr,method,arg) \
    threadptr = new std::thread(method, arg)
#define mps_thread_id_t std::thread::id
#define mps_thread_self() std::this_thread::get_id()
#define mps_thread_equal(thread, self) ((thread).get_id() == self)
#define mps_thread_create(threadptr,method,start_arg) 
#define mps_thread_yield() std::this_thread::yield()

#define mps_mutex_t std::mutex
#define mps_initialized_mutex(mutexvar,type) std::mutex mutexvar
#define mps_mutex_lock(mutexvar) (mutexvar).lock()
#define mps_mutex_unlock(mutexvar) (mutexvar).unlock()
#define mps_mutex_init(mutexvar)
#define mps_mutex_destroy(mutexvar)
#define mps_once_t std::once_flag
#define mps_static_initialized_once(var,type) static std::once_flag var
#define mps_once(once_key,create_key_method) std::call_once(once_key,create_key_method)

#define mps_cond_t std::condition_variable
#define mps_unique_mutex_t std::unique_mutex<std::mutex>
#define mps_unique_mutex_ptr(assocmutex) mps_unique_mutex_t * assocmutex

#define mps_cond_t std::condition_variable
#define mps_cond_init(cond)
#define mps_cond_destroy(cond)
#define mps_cond_broadcast(cond,assocmutex) \
{ \
    std::unique_lock<std::mutex> signal_lock(assocmutex); \
    (cond).notify_all(); \
}
#define mps_cond_signal(cond,assocmutex) \
{ \
    std::unique_lock<std::mutex> signal_lock(assocmutex); \
    (cond).notify_one(); \
}
#define mps_cond_wait(cond,assocmutex,mutexvar,donecond) \
{ \
    mps_mutex_unlock(mutexvar); \
    std::unique_lock<std::mutex> signal_lock(assocmutex); \
    (cond).wait(signal_lock, [&] { return donecond; }); \
} \
mps_mutex_lock(mutexvar)
#define mps_static_tls_struct_init(keyinit) static thread_local bool keyinit = false
#define mps_tls_struct(key) static thread_local mps_tls key
#define mps_tls_create(key,keytype,keyinit,ptr) \
    keytype* ptr = &key; \
    keyinit = true
#define mps_tls_setspecific(key,keyinit,ptr)
#define mps_tls_getspecific(key,keyinit) \
keyinit? &key:NULL
#define mps_tls_key_create(key,clean_method)

#else

//  implement threading and mutexes using unix pthread
#include <pthread.h>
#define mps_mutex_t pthread_mutex_t
#define mps_initialized_mutex(mutexvar,type) pthread_mutex_t mutexvar = type
#define mps_thread_join(threadptr) \
  pthread_join((*(threadptr)),NULL); \
  mps_del_obj(threadptr)
#define mps_thread_processing_exit(thread) pthread_exit(NULL)
#define mps_thread_t pthread_t
#define mps_start_new_thread(threadptr,method,arg) \
  threadptr = new mps_thread_t; \
  pthread_create(threadptr,NULL,method,arg)
#define mps_thread_id_t pthread_t
#define mps_thread_self() pthread_self()
#define mps_thread_equal(thread, self) pthread_equal (thread, self)
#define mps_thread_create(threadptr,method,start_arg) pthread_create(&threadptr,NULL,method,start_arg)
#define mps_thread_yield() sched_yield()

#define mps_unique_mutex(mutexvar) pthread_mutex_t mutexvar
#define mps_static_tls_struct_init(keyinit)
#define mps_tls_struct(key) static pthread_key_t key
#define mps_tls_create(key,keytype,keyinit,ptr) mps_new_obj(keytype, ptr, sizeof(keytype))
#define mps_tls_setspecific(key,keyinit,ptr) pthread_setspecific(key,ptr)
#define mps_tls_getspecific(key,keyinit) pthread_getspecific(key)
#define mps_tls_key_create(key,clean_method) pthread_key_create(key,clean_method)

#define mps_mutex_lock(mutexvar) pthread_mutex_lock(&(mutexvar))
#define mps_mutex_unlock(mutexvar) pthread_mutex_unlock(&(mutexvar))
#define mps_mutex_init(mutexvar) pthread_mutex_init(&(mutexvar),NULL)
#define mps_mutex_destroy(mutexvar) pthread_mutex_destroy(&mutexvar)
#define mps_once_t pthread_once_t
#define mps_static_initialized_once(var,type) static pthread_once_t var = type
#define mps_once(once_key,create_key_method) pthread_once(&once_key,create_key_method)

#define mps_cond_t pthread_cond_t
#define mps_unique_mutex_t pthread_t
#define mps_unique_mutex_ptr(assocmutex) mps_unique_mutex_t * assocmutex
#define mps_new_unique_mutex(name,uniquemutex) new mps_unique_mutex_t
#define mps_del_unique_mutex(name,assocmutex) 
#define mps_cond_init(cond) pthread_cond_init(&(cond), NULL)
#define mps_cond_destroy(cond) pthread_cond_destroy(&(cond), NULL)
#define mps_cond_broadcast(cond,assocmutex) pthread_cond_broadcast(&(cond))
#define mps_cond_signal(cond,assocmutex) pthread_cond_signal(&(cond))
#define mps_cond_wait(cond,assocmutex,mutexvar,donecond) pthread_cond_wait(&(cond),&(mutexvar))
#endif

// this is the same for both mutex handlers
#define mps_mutex_tracked_lock(mutexvar) mps_tracked_lock(mutexvar)
#define mps_mutex_tracked_unlock(mutexvar) mps_tracked_unlock(mutexvar)

#include <mps_unistd.h> 

/* This header should be included first since it contains all the forward
* declaration that must be available in the others. */
#include <mps/types.h>

/* Include more types that are needed in the declarations of the generic
* functions, such as DPE and MP. */
#include <mps/mt.h>
#include <mps/gmptools.h>
#include <mps/mpc.h>
#include <mps/link.h>
#include <mps/polynomial.h>

/* Public types, mostly custom polynomial types such as Chebyshev, Monomial and
* Secular equations. */
#include <mps/matrix.h>
#include <mps/chebyshev.h>
#include <mps/monomial-matrix-poly.h>
#include <mps/monomial-poly.h>
#include <mps/secular-equation.h>
#include <mps/nroots-polynomial.h>
#include <mps/regeneration-driver.h>

/* Public interface functions for MPSolve */
#include <mps/approximation.h>
#include <mps/context.h>
#include <mps/debug.h>
#include <mps/interface.h>
#include <mps/parser.h>
#include <mps/version.h>

/* Private inclusions. Please note that these header files may not be distributed with
* MPSolve, so it's safe to use them only for internal functions. */
#ifdef _MPS_PRIVATE

#ifndef getline
MPS_BEGIN_DECLS
    ssize_t getline(char** lineptr, size_t* n, FILE* stream);
MPS_END_DECLS
#endif
#include <mps/private/system/abstract-input-stream.h>
#include <mps/private/system/file-input-stream.h>
#include <mps/private/system/memory-file-stream.h>
#include <mps/private/aberth.h>
#include <mps/private/algorithms.h>
#include <mps/private/cluster.h>
#include <mps/private/convex.h>
#include <mps/private/data.h>
#include <mps/private/hessenberg-determinant.h>
#include <mps/private/horner.h>
#include <mps/private/jacobi-aberth.h>
#include <mps/private/improve.h>
#include <mps/private/input-buffer.h>
#include <mps/private/input-output.h>
#include <mps/private/list.h>
#include <mps/private/mandelbrot-user.h>
#include <mps/private/newton.h>
#include <mps/private/options.h>
#include <mps/private/radii.h>
#include <mps/private/secular-evaluation.h>
#include <mps/private/solve.h>
#include <mps/private/sort.h>
#include <mps/private/starting.h>
#include <mps/private/starting-configuration.h>
#include <mps/private/threading.h>
#include <mps/private/tools.h>
#include <mps/private/touch.h>
#include <mps/private/utils.h>
#include <mps/private/formal/formal-monomial.h>
#include <mps/private/formal/formal-polynomial.h>
#include <mps/private/secular-regeneration.h>
#endif

#define ERRMSG_SIZE 5000
void mps_fatal_exit(char* errmsg);
bool mps_tracked_lock(mps_mutex_t& target_mutex);
bool mps_tracked_unlock(mps_mutex_t& target_mutex);

#define mps_new_obj(objname,varname,varsize) \
objname * varname = new objname

#define mps_new_defined_obj(objname,varname,varsize) \
varname = new objname

#pragma warning( disable : 4150 )
#define mps_del_obj(objname) \
delete objname; \
objname = NULL

#define mps_new_array_obj(objname,varname,varsize,count) \
objname * varname = new objname[count]

#define mps_new_defined_array_obj(objname,varname,varsize,count) \
varname = new objname[count]

#define mps_del_array_obj(objname) \
delete[] objname; \
objname = NULL

#define mps_malloc(size) malloc(size)
#define mps_realloc(ptr, size) realloc(ptr, size);
#define mps_free(ptr) free(ptr)

#endif                          /* endef MPSCORE_H */
