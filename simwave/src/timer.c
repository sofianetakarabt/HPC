///
/// \copyright Copyright 2015 Issam Said. All rights reserved.
/// This file is part of \b simwave.
///
/// \b simwave is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// \b simwave is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with \b simwave.  If not, see <http://www.gnu.org/licenses/>.
///
/// \author Issam Said
/// \file timer.c
/// \version $Id: timer.c 1855 2012-12-29 00:45:28Z isaid $
/// \brief This file contains the implementation of the timer.
///

#include "timer.h"
#include "config.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#include <sys/types.h>
#include <sys/time.h>
#endif //  _WIN32

#include "stdio.h"

static timing_t _elapsed;

#ifdef _WIN32
static timing_t _freq=0;
#elif _GETTIMEOFDAY
static struct timeval  _tv_start;
static struct timeval  _tv_end;
#elif _CLOCK_GETTIME 
static struct timespec _ts_start;
static struct timespec _ts_end;
#else 
#error unknown timer
#endif // _WIN32

char *timing_unit_string(timing_unit_t unit) {
    switch(unit) {
    case SECONDS:
        return "s";
    case MILLISECONDS:
        return "ms";
    case MICROSECONDS:
        return "us";
    case NANOSECONDS:
        return "ns";
    default:
        return "unknown";
    }
}

double timing_to_seconds(timing_unit_t unit) {
    return 1.e-9*((double)unit);
}

void init_timer(){}

INLINE void start_timer() {
#ifdef _WIN32
    QueryPerformanceCounter((LARGE_INTEGER*)&_tstart);
#elif _GETTIMEOFDAY
    gettimeofday(&_tv_start, NULL);
#elif _CLOCK_GETTIME
    clock_gettime(CLOCK_REALTIME, &_ts_start);
#else
#error unknown timer
#endif // _WIN32
}

INLINE double stop_timer(timing_unit_t unit) {
#ifdef _WIN32
    QueryPerformanceCounter((LARGE_INTEGER *)&_tend);	
#elif _GETTIMEOFDAY
    gettimeofday(&_tv_end, NULL);
    _elapsed  = (_tv_end.tv_sec  -  _tv_start.tv_sec) * 1e9; 
    _elapsed += (_tv_end.tv_usec - _tv_start.tv_usec) * 1e3;
#elif _CLOCK_GETTIME
    clock_gettime( CLOCK_REALTIME, &_ts_end);
    _elapsed  = (_ts_end.tv_sec  -  _ts_start.tv_sec) * 1e9;
    _elapsed += (_ts_end.tv_nsec - _ts_start.tv_nsec); 
#else 
#error unknown timer
#endif // _WIN32

#ifdef _WIN32
    return ((double)(_elapsed)/
           (double)_freq)*1.e9/(double)unit;
#else
    return (double)(_elapsed)/(double)unit;
#endif // _WIN32
}
