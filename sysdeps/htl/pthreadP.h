/* Declarations of internal pthread functions used by libc.  Hurd version.
   Copyright (C) 2016-2025 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#ifndef _PTHREADP_H
#define _PTHREADP_H	1

#define __PTHREAD_HTL

#include <pthread.h>
#include <link.h>
#include <bits/cancelation.h>

/* Attribute to indicate thread creation was issued from C11 thrd_create.  */
#define ATTR_C11_THREAD ((void*)(uintptr_t)-1)

extern void __pthread_init_static_tls (struct link_map *) attribute_hidden;

/* These represent the interface used by glibc itself.  */

extern int __pthread_barrier_destroy (pthread_barrier_t *__barrier);
libc_hidden_proto (__pthread_barrier_destroy)
extern int __pthread_barrier_init (pthread_barrier_t *__barrier,
				const pthread_barrierattr_t *__attr,
				unsigned __count);
libc_hidden_proto (__pthread_barrier_init)
extern int __pthread_barrier_wait (pthread_barrier_t *__barrier);
libc_hidden_proto (__pthread_barrier_wait)
extern int __pthread_barrierattr_destroy (pthread_barrierattr_t *__attr);
libc_hidden_proto (__pthread_barrierattr_destroy)
extern int __pthread_barrierattr_getpshared (const pthread_barrierattr_t *__attr,
					   int *__pshared);
libc_hidden_proto (__pthread_barrierattr_getpshared)
extern int __pthread_barrierattr_init (pthread_barrierattr_t *__attr);
libc_hidden_proto (__pthread_barrierattr_init)
extern int __pthread_barrierattr_setpshared (pthread_barrierattr_t *__attr,
					   int __pshared);
libc_hidden_proto (__pthread_barrierattr_setpshared)
extern int __pthread_mutex_init (pthread_mutex_t *__mutex, const pthread_mutexattr_t *__attr);
libc_hidden_proto (__pthread_mutex_init)
extern int __pthread_mutex_clocklock (pthread_mutex_t *__mutex, clockid_t __clockid,
				    const struct timespec *__abstime);
libc_hidden_proto (__pthread_mutex_clocklock)
extern int __pthread_mutex_consistent (pthread_mutex_t *__mtxp);
libc_hidden_proto (__pthread_mutex_consistent)
extern int __pthread_mutex_lock (pthread_mutex_t *__mutex);
libc_hidden_proto (__pthread_mutex_lock)
extern int __pthread_mutex_getprioceiling (const pthread_mutex_t *__mutex,
					 int *__prioceiling);
libc_hidden_proto (__pthread_mutex_getprioceiling)
extern int __pthread_mutex_setprioceiling (pthread_mutex_t *__mutex,
					 int __prio, int *__oldprio);
libc_hidden_proto (__pthread_mutex_setprioceiling)
extern int __pthread_mutex_trylock (pthread_mutex_t *_mutex);
libc_hidden_proto (__pthread_mutex_trylock)
extern int __pthread_mutex_timedlock (pthread_mutex_t *__mutex,
     const struct timespec *__abstime);
libc_hidden_proto (__pthread_mutex_timedlock)
extern int __pthread_mutex_unlock (pthread_mutex_t *__mutex);
libc_hidden_proto (__pthread_mutex_unlock)
extern int __pthread_mutexattr_getpshared(const pthread_mutexattr_t *__restrict __attr,
					int *__restrict __pshared);
libc_hidden_proto (__pthread_mutexattr_getpshared)
extern int __pthread_mutexattr_setpshared(pthread_mutexattr_t *__attr,
					int __pshared);
libc_hidden_proto (__pthread_mutexattr_setpshared)

/* Get the robustness flag of the mutex attribute ATTR.  */
extern int __pthread_mutexattr_getrobust (const pthread_mutexattr_t *__attr,
					int *__robustness);
libc_hidden_proto (__pthread_mutexattr_getrobust)
extern int __pthread_mutexattr_getrobust_np (const pthread_mutexattr_t *__attr,
					   int *__robustness);
libc_hidden_proto (__pthread_mutexattr_getrobust_np)
extern int __pthread_mutexattr_setrobust_np (pthread_mutexattr_t *__attr,
					   int __robustness);
libc_hidden_proto (__pthread_mutexattr_setrobust_np)

/* Set the robustness flag of the mutex attribute ATTR.  */
extern int __pthread_mutexattr_setrobust (pthread_mutexattr_t *__attr,
					int __robustness);
libc_hidden_proto (__pthread_mutexattr_setrobust)

extern int __pthread_mutexattr_getprioceiling(const pthread_mutexattr_t *__restrict __attr,
					    int *__restrict __prioceiling);
libc_hidden_proto (__pthread_mutexattr_getprioceiling)
extern int __pthread_mutexattr_setprioceiling(pthread_mutexattr_t *__attr,
					    int __prioceiling);
libc_hidden_proto (__pthread_mutexattr_setprioceiling)

extern int __pthread_mutexattr_getprotocol(const pthread_mutexattr_t *__restrict __attr,
					 int *__restrict __protocol);
libc_hidden_proto (__pthread_mutexattr_getprotocol)
extern int __pthread_mutexattr_setprotocol(pthread_mutexattr_t *__attr,
					 int __protocol);
libc_hidden_proto (__pthread_mutexattr_setprotocol)

extern int __pthread_mutexattr_gettype(const pthread_mutexattr_t *__restrict __attr,
				     int *__restrict __type);
libc_hidden_proto (__pthread_mutexattr_gettype)
extern int __pthread_mutexattr_settype(pthread_mutexattr_t *__attr,
				     int __type);
libc_hidden_proto (__pthread_mutexattr_settype)
extern int __pthread_rwlock_clockrdlock (pthread_rwlock_t *__rwlock,
				       clockid_t __clockid, const struct timespec *__abstime);
libc_hidden_proto (__pthread_rwlock_clockrdlock)
extern int __pthread_rwlock_clockwrlock (pthread_rwlock_t *__rwlock,
				       clockid_t __clockid, const struct timespec *__abstime);
libc_hidden_proto (__pthread_rwlock_clockwrlock)
extern int __pthread_rwlock_timedrdlock (struct __pthread_rwlock *__rwlock,
				       const struct timespec *__abstime);
libc_hidden_proto (__pthread_rwlock_timedrdlock)
extern int __pthread_rwlock_timedwrlock (struct __pthread_rwlock *__rwlock,
				       const struct timespec *__abstime);
libc_hidden_proto (__pthread_rwlock_timedwrlock)
extern int __pthread_rwlockattr_destroy (pthread_rwlockattr_t *__attr);
libc_hidden_proto (__pthread_rwlockattr_destroy)
extern int __pthread_rwlockattr_getpshared (const pthread_rwlockattr_t *__attr,
					  int *__pshared);
libc_hidden_proto (__pthread_rwlockattr_getpshared)
extern int __pthread_rwlockattr_init (pthread_rwlockattr_t *__attr);
libc_hidden_proto (__pthread_rwlockattr_init)
extern int __pthread_rwlockattr_setpshared (pthread_rwlockattr_t *__attr,
					  int __pshared);
libc_hidden_proto (__pthread_rwlockattr_setpshared)

extern int __pthread_cond_init (pthread_cond_t *cond,
				const pthread_condattr_t *cond_attr);
libc_hidden_proto (__pthread_cond_init)
extern int __pthread_cond_signal (pthread_cond_t *cond);
libc_hidden_proto (__pthread_cond_signal);
extern int __pthread_cond_broadcast (pthread_cond_t *cond);
libc_hidden_proto (__pthread_cond_broadcast);
extern int __pthread_cond_wait (pthread_cond_t *cond, pthread_mutex_t *mutex);
libc_hidden_proto (__pthread_cond_wait);
extern int __pthread_cond_timedwait (pthread_cond_t *cond,
				     pthread_mutex_t *mutex,
				     const struct timespec *abstime);
libc_hidden_proto (__pthread_cond_timedwait);
extern int __pthread_cond_clockwait (pthread_cond_t *cond,
				     pthread_mutex_t *mutex,
				     clockid_t clockid,
				     const struct timespec *abstime)
  __nonnull ((1, 2, 4));
libc_hidden_proto (__pthread_cond_clockwait);
extern int __pthread_cond_destroy (pthread_cond_t *cond);
libc_hidden_proto (__pthread_cond_destroy);

extern int __pthread_setcanceltype (int __type, int *__oldtype);
libc_hidden_proto (__pthread_setcanceltype);
extern int __pthread_sigmask (int, const sigset_t *, sigset_t *);
libc_hidden_proto (__pthread_sigmask);

typedef struct __cthread *__cthread_t;
typedef int __cthread_key_t;
typedef void *	(*__cthread_fn_t)(void *__arg);

__cthread_t __cthread_fork (__cthread_fn_t, void *);
int __pthread_create (pthread_t *newthread,
		      const pthread_attr_t *attr,
		      void *(*start_routine) (void *), void *arg);

void __cthread_detach (__cthread_t);
int __pthread_detach (pthread_t __threadp);
void __pthread_exit (void *value) __attribute__ ((__noreturn__));
int __pthread_join (pthread_t, void **);
int __cthread_keycreate (__cthread_key_t *);
int __cthread_getspecific (__cthread_key_t, void **);
int __cthread_setspecific (__cthread_key_t, void *);
int __pthread_key_create (pthread_key_t *key, void (*destr) (void *));
void *__pthread_getspecific (pthread_key_t key);
int __pthread_setspecific (pthread_key_t key, const void *value);
int __pthread_key_delete (pthread_key_t key);
int __pthread_once (pthread_once_t *once_control, void (*init_routine) (void));

int __pthread_getattr_np (pthread_t, pthread_attr_t *);
int __pthread_attr_getstackaddr (const pthread_attr_t *__restrict __attr,
				 void **__restrict __stackaddr);
libc_hidden_proto (__pthread_attr_getstackaddr)
int __pthread_attr_setstackaddr (pthread_attr_t *__attr, void *__stackaddr);
libc_hidden_proto (__pthread_attr_setstackaddr)
int __pthread_attr_getstacksize (const pthread_attr_t *__restrict __attr,
				 size_t *__restrict __stacksize);
libc_hidden_proto (__pthread_attr_getstacksize)
int __pthread_attr_setstacksize (pthread_attr_t *__attr, size_t __stacksize);
libc_hidden_proto (__pthread_attr_setstacksize)
int __pthread_attr_setstack (pthread_attr_t *__attr, void *__stackaddr,
			     size_t __stacksize);
libc_hidden_proto (__pthread_attr_setstack)
int __pthread_attr_getstack (const pthread_attr_t *, void **, size_t *);
libc_hidden_proto (__pthread_attr_getstack)
void __pthread_testcancel (void);
int __pthread_attr_init (pthread_attr_t *attr);
int __pthread_condattr_init (pthread_condattr_t *attr);

#define __pthread_raise_internal(__sig) raise (__sig)

libc_hidden_proto (__pthread_self)
libc_hidden_proto (__pthread_attr_init)
libc_hidden_proto (__pthread_condattr_init)

#if IS_IN (libpthread)
hidden_proto (__pthread_create)
hidden_proto (__pthread_detach)
hidden_proto (__pthread_key_create)
hidden_proto (__pthread_getspecific)
hidden_proto (__pthread_setspecific)
hidden_proto (__pthread_get_cleanup_stack)
#endif

#if !defined(__NO_WEAK_PTHREAD_ALIASES) && !IS_IN (libpthread)
# ifdef weak_extern
weak_extern (__pthread_exit)
# else
#  pragma weak __pthread_exit
# endif
#endif

#define ASSERT_TYPE_SIZE(type, size) 					\
  _Static_assert (sizeof (type) == size,				\
		  "sizeof (" #type ") != " #size)

 /* Special cleanup macros which register cleanup both using
    __pthread_cleanup_{push,pop} and using cleanup attribute.  This is needed
    for qsort, so that it supports both throwing exceptions from the caller
    sort function callback (only cleanup  attribute works there) and
    cancellation of the thread running the callback if the callback or some
    routines it calls don't have unwind information.
    TODO: add support for cleanup routines.  */
#ifndef pthread_cleanup_combined_push
# define pthread_cleanup_combined_push  __pthread_cleanup_push
#endif
#ifndef pthread_cleanup_combined_pop
# define pthread_cleanup_combined_pop   __pthread_cleanup_pop
#endif

#endif	/* pthreadP.h */
