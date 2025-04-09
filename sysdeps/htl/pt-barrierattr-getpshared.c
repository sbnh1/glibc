/* pthread_barrierattr_getpshared.  Generic version.
   Copyright (C) 2002-2025 Free Software Foundation, Inc.
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
   License along with the GNU C Library;  if not, see
   <https://www.gnu.org/licenses/>.  */

#include <pthread.h>
#include <pt-internal.h>
#include <shlib-compat.h>

int
__pthread_barrierattr_getpshared (const pthread_barrierattr_t *attr,
				int *pshared)
{
  *pshared = attr->__pshared;
  return 0;
}
libc_hidden_def (__pthread_barrierattr_getpshared)
versioned_symbol (libc, __pthread_barrierattr_getpshared, pthread_barrierattr_getpshared, GLIBC_2_42);

#if OTHER_SHLIB_COMPAT (libpthread, GLIBC_2_12, GLIBC_2_42)
compat_symbol (libpthread, __pthread_barrierattr_getpshared, pthread_barrierattr_getpshared, GLIBC_2_12);
#endif
