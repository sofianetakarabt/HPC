#ifndef _VELOCITY_H_
#define _VELOCITY_H_
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
/// \file velocity.h
/// \version $Id: velocity.h 2767 2015-10-26 08:02:36Z isaid $
/// \brief A velocity model generator.
///
/// Generates a phase velocity model for the simulation.
///

#include "config.h"

/// \brief Generates a three dimensional velocity model in m/s
/// \param vtab is an array that will contain the model terms
/// \param dimx is the width of the domain
/// \param dimy is the height of the domain
/// \param dimz is the depth of the domain
/// \param vmin is the minimum velocity required
/// \param vmax is the maximum velocity tolerated
/// \param layers is the maximum velocity layers number
void velocity_generate_model(DATATYPE* vtab,
                             int dimx, int dimy, int dimz,
                             DATATYPE vmin, DATATYPE vmax, unsigned int layers);

/// \brief Calculates the maximum velocity of a three dimensional model
/// \param vtab is an array that contains the model terms
/// \param dimx is the width of the domain
/// \param dimy is the height of the domain
/// \param dimz is the depth of the domain
/// \return the maximum velocity
DATATYPE velocity_get_max(DATATYPE* vtab,
                          int dimx, int dimy, int dimz);

#endif
