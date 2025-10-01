/*
 * This file was originally is part of MPSolve 3.2.1
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>

/*
	These methods have been revised to integrate with the num math library,
	in particular the use of realnum and compnum

	The num math library has Copyright(C) 2025 Drew Meyer ddmeyer@attglobal.net
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

#define MPS_BEGIN_DECLS
#define MPS_END_DECLS

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//#define MPS_USE_PTHREADS
#ifndef MPS_USE_PTHREADS

//std::unique_lock<std::mutex> lock(my_mutex)

//std::mutex mymutex;
//std::unique_lock<std::mutex> lock(mymutex, std::defer_lock);

//  implement threading and mutexes using the Windows std lib methods
#include <condition_variable>
#include <mutex>
#include <thread>
#define mps_mutex_t std::mutex
#define mps_initialized_mutex(mutexvar,type) std::mutex mutexvar
#define mps_initialized_mutex_ptr_t std::mutex
#define mps_unique_mutex(mutexvar) std::mutex mutexvar
#define mps_once_t std::once_flag
#define mps_initialized_once(var,type) std::once_flag var
#define mps_with_lock(pmutex, code) \
{ \
    mps_mutex_lock (&pmutex); \
    code \
    mps_mutex_unlock (&pmutex);  \
}
#define mps_thread_once(once_key,create_key_method) std::call_once(once_key,create_key_method)
#define mps_mutex_lock(mutex) (*mutex).lock()
#define mps_mutex_unlock(mutex) (*mutex).unlock()
#define mps_mutex_init(mutex)
#define mps_mutex_destroy(mutex)
#define mps_thread_cond_t std::condition_variable
#define mps_thread_cond_wait(cond,mutex) (cond).wait(mutex)
#define mps_thread_cond_init(cond)
#define mps_thread_cond_broadcast(arg) (arg).notify_all()
#define mps_thread_cond_signal(cond) (*cond).notify_one()
#define mps_thread_join(thread) (thread).join()
#define mps_thread_exit(thread) \
delete thread;\
mps_unreg(thread)
#define mps_thread_t std::thread
#define mps_start_new_thread(thread_ptr,method,arg)\
  *thread_ptr = std::thread(method, arg);\
  mps_reg("thread_ptr", thread_ptr, sizeof(mps_thread_t))
#define mps_thread_id_t std::thread::id
#define mps_thread_self() std::this_thread::get_id()
#define mps_thread_equal(thread, self) ((thread).get_id() == self)
#define mps_thread_create(thread_ptr,method,start_arg) 
#define mps_thread_static_struct_init(keyinit) static thread_local bool keyinit = false
#define mps_thread_struct(key)\
  thread_local mps_tls key
#define mps_thread_yield() std::this_thread::yield()
/*
#define mps_thread_setspecific(key,keyinit,ptr)\
  memcpy(ptr,key,sizeof(mps_tls));\
  mps_free(ptr);\
  keyinit = true;\
  ptr = key
#define mps_thread_getspecific(keyptr,keyinit)\
  ((keyinit) ? keyptr : NULL)
#define mps_thread_key_create(key,clean_method)
*/

#else

//  implement threading and mutexes using unix pthread
#define mps_mutex_t pthread_mutex_t
#define mps_initialized_mutex(mutexvar,type) pthread_mutex_t mutexvar = type
#define mps_initialized_mutex_ptr_t pthread_mutex_t
#define mps_unique_mutex(mutexvar) pthread_mutex_t mutexvar
#define mps_once_t pthread_once_t
#define mps_initialized_once(var,type) pthread_once_t var = type
#define mps_with_lock(pmutex, code) { \
    mps_mutex_lock (&pmutex); \
    code \
    mps_mutex_unlock (&pmutex); \
  }
#define mps_thread_once(once_key,create_key_method) pthread_once(&once_key,create_key_method)
#define mps_mutex_lock(mutex) pthread_mutex_lock(mutex)
#define mps_mutex_unlock(mutex) pthread_mutex_unlock(mutex)
#define mps_mutex_init(mutex) pthread_mutex_init(mutex,NULL)
#define mps_mutex_destroy(mutex) pthread_mutex_destroy(mutex)
#define mps_thread_cond_t pthread_cond_t
#define mps_thread_cond_wait(cond,mutex) pthread_cond_wait(&(cond),&(mutex))
#define mps_thread_cond_init(cond) pthread_cond_init(&cond, NULL)
#define mps_thread_cond_broadcast(arg) pthread_cond_broadcast(&arg)
#define mps_thread_cond_signal(cond) pthread_cond_signal(cond)
#define mps_thread_join(thread) pthread_join((thread),NULL)
#define mps_thread_exit(thread)\
pthread_exit(NULL);\
mps_free(thread)
#define mps_thread_t pthread_t
#define mps_start_new_thread(thread_ptr,method,arg)\
thread_ptr = mps_new (mps_thread_t);\
pthread_create(thread_ptr,NULL,method,arg)
#define mps_thread_id_t pthread_t
#define mps_thread_self() pthread_self()
#define mps_thread_equal(thread, self) pthread_equal (thread, self)
#define mps_thread_create(thread_ptr,method,start_arg) pthread_create(&thread_ptr,NULL,method,start_arg)
#define mps_thread_static_struct_init(keyinit) 
#define mps_thread_struct(key) pthread_key_t key
#define mps_thread_yield() sched_yield()
#endif

#include <pthread.h>
#define mps_thread_setspecific(keyptr,keyinit,ptr) pthread_setspecific(*keyptr,ptr)
#define mps_thread_getspecific(keyptr,keyinit) pthread_getspecific(*keyptr)
#define mps_thread_key_create(key,clean_method) pthread_key_create(key,clean_method)

#include <mps/mps_unistd.h> 

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

#ifndef getline
MPS_BEGIN_DECLS
ssize_t getline (char **lineptr, size_t *n, FILE *stream);
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

#include <num_defines.h>
#define ERRMSG_SIZE 5000
extern void mpsolve_fatal_exit(char* errmsg);

#ifdef USE_NUM_MEM
#define mps_malloc(size) NUM_MEM_MALLOC_MPSOLVE(size)
#define mps_realloc(ptr,size) NUM_MEM_REALLOC_MPSOLVE(ptr, size);
#define mps_free(ptr) NUM_MEM_FREE_MPSOLVE(ptr)
#define mps_reg(desc,ptr,size) NUM_MEM_REG_MPSOLVE(desc,ptr,size)
#define mps_reg_mem_leak(desc,ptr,size) num_register_mem_leak(desc,ptr,size,false)
#define mps_unreg(ptr) NUM_MEM_UNREG_MPSOLVE(ptr)
#else
#define mps_malloc(size) malloc(size)
#define mps_realloc(ptr, size) realloc(ptr, size);
#define mps_free(ptr) NUM_MEM_FREEONLY_MPSOLVE(ptr)   // do a free allowing ptr to be null and clearing it to NULL aterwards
#define mps_reg(desc,ptr,size)
#define mps_reg_mem_leak(desc,ptr,size)
#define mps_unreg(ptr)
#endif

#endif                          /* endef MPSCORE_H */
