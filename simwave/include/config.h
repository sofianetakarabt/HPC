#ifndef _CONFIG_H_
#define _CONFIG_H_
///
/// \copyright Copyright 2015 Issam Said. All rights reserved.
/// This file is part of \b simwave.
///
/// \b simwave is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// simwave is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
///
/// \author Issam Said
/// \file config.h
/// \version $Id: config.h 2781 2015-10-26 09:17:02Z isaid $
/// \brief Contains the major settings of project.
///

#ifdef _OPENMP
#include <omp.h>
#endif 

#ifdef _WIN32
#define INLINE __forceinline
#else
#define INLINE inline
#endif // _WIN32

/// \def DATATYPE
/// \brief Abstarcts the type of data.
/// It can be either float or double.
#define DATATYPE float

/// \def _GETTIMEOFDAY
/// \brief Set linux gettimeofday wall-clock timer
/// to be used by the timing routines \see timer.h. 
#define _GETTIMEOFDAY 1

#if defined(__APPLE__) || defined(__MACOSX)
#define _GETTIMEOFDAY 1
#endif

#endif //  _CONFIG_H_
