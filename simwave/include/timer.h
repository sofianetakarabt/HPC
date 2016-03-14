#ifndef _TIMER_H_
#define _TIMER_H_
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
/// along with \b simwave. If not, see <http://www.gnu.org/licenses/>.
///
/// \author Issam Said
/// \file timer.h
/// \version $Id: timer.h 2767 2015-10-26 08:02:36Z isaid $
/// \brief Timing utilities.
///
/// Contains a cross-pltaform timing set of routines. In order to time an execution
/// call start_timer before execution and stop_timer after
///

/// \brief Defines a type to describe a timing units
typedef enum {
    SECONDS      = 1000000000,
    MILLISECONDS = 1000000   ,
    MICROSECONDS = 1000      ,
    NANOSECONDS  = 1         ,
} timing_unit_t;


/// \brief Converts a timing unit to a string 
/// \param unit the timing unit
/// \return A string that describes the timing unit 
char *timing_unit_string(timing_unit_t unit);

/// \brief Calculates a scale to convert to seconds
/// \param unit the timing unit
/// \return A double that contains the converting factor
double timing_to_seconds(timing_unit_t unit);


#ifdef _WIN32
/// \brief Defines a type that contains cpu ticks number
typedef __int64 timing_t;
#else
/// \brief Defines a type that contains cpu ticks number
typedef long long timing_t;
#endif // _WIN32

/// \brief Detects the timer resolution latency
void init_timer();

/// \brief Starts the timer
void start_timer();

/// \brief Stops the timer and returns the elapsed time since
/// \ref start_timer was last called
/// \param unit the unit of the measured time
/// \return the elapsed time
double stop_timer(timing_unit_t unit);

#endif // _TIMER_H_
