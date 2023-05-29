
/*
Copyright (C) 2008 Peter Langfelder; parts based on R by R Development team

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/


#ifndef __conditionalThreading_h__
#define __conditionalThreading_h__

#define MxThreads      128

#ifdef WITH_THREADS

  // #warning Including pthread headers.

  #include <unistd.h>
  #include <pthread.h>

#else

  // define fake pthread functions so we don't have to put a #ifdef everywhere
  //
  // This prevents competing definitions of pthread types to be included
  #define _BITS_PTHREADTYPES_H

  typedef int pthread_mutex_ref; // replaced pthread_mutex_t to avoid clash with CLANG interpreter
  typedef int pthread_ref;
  // in the original code this was called "pthread_attr_t",
  // which causes issues with arm64 macs, as the OS uses this variable internally already, hence the change
  typedef int pthread_attr;

  #define PTHREAD_MUTEX_INITIALIZER 0

  static inline void pthread_mutex_lock ( pthread_mutex_ref * lock ) { }
  static inline void pthread_mutex_unlock ( pthread_mutex_ref * lock ) { }

  static inline int pthread_join ( pthread_ref t, void ** p) { return 0; }

#endif


// Conditional pthread routines

static inline void pthread_mutex_lock_c( pthread_mutex_ref * lock, int threaded)
{
  if (threaded) pthread_mutex_lock(lock);
}

static inline void pthread_mutex_unlock_c(pthread_mutex_ref * lock, int threaded)
{
  if (threaded) pthread_mutex_unlock(lock);
}

static inline int pthread_create_c(pthread_ref *thread, const pthread_attr *attr,
    void *(*start_routine)(void*), void *arg, int threaded)
{
  #ifdef WITH_THREADS
  if (threaded)
    return pthread_create(thread, attr, start_routine, arg);
  else
  #endif
    (*start_routine)(arg);
  return 0;
}

static inline int pthread_join_c(pthread_ref thread, void * * value_ptr, int threaded)
{
  if (threaded) return pthread_join(thread, (void * *) value_ptr);
  return 0;
}


#endif