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
/// \file velocity.c
/// \version $Id: velocity.c 1864 2012-12-29 23:09:23Z isaid $
/// \brief This file contains the implementation of the velocity generator.
///

#include "velocity.h"
void velocity_generate_model(DATATYPE* vtab,
                             int dimx, int dimy, int dimz,
                             DATATYPE vmin, DATATYPE vmax, unsigned int layers) {
    int x, y, z, i;
    DATATYPE delta = layers == 1 ? (vmax - vmin) : (vmax - vmin)/(layers-1);
    for (i=0; i<layers; ++i) {
        for (y=(i*dimy)/layers; y<((i+1)*dimy)/layers; y++)
            for (z=0; z<dimz; z++)
                for (x=0; x<dimx; x++) {
                    vtab[dimx*(dimy*z + y) + x] = vmin + (i*delta);
                }
    }
}

#define MAX(a,b) (a>b?a:b)

DATATYPE velocity_get_max(DATATYPE* vtab,
                          int dimx, int dimy, int dimz) {
    DATATYPE max_velocity = 0.0;
    int x,y,z;
    for (z=0; z<dimz; z++ )
        for (y=0; y<dimy; y++ )
            for (x=0; x<dimx; x++ )
                max_velocity = MAX(vtab[dimx*(dimy*z + y) + x], max_velocity); 
    return max_velocity;
}
