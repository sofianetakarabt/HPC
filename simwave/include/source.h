#ifndef _SOURCE_H_
#define _SOURCE_H_
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
/// \file source.h
/// \version $Id: source.h 2767 2015-10-26 08:02:36Z isaid $
/// \brief A wave source generator.
///
/// Create a periodic signal that simulates waves pressure impulses.
///

#include "config.h"

/// \brief Creates a Ricker Wavelet to be used as a wave source.
/// \see http://subsurfwiki.org/wiki/Ricker_wavelet
/// \param source an array that will contains the source discretized terms
/// \param dt is the duration of a time step
/// \param time_steps the number of the simulation time steps
/// \param fmax is the max frequency of the input signal
void source_ricker_wavelet(DATATYPE* source,
                           DATATYPE dt, int time_steps, DATATYPE fmax);

#endif
