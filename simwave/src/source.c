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
/// \file source.c
/// \version $Id: source.c 1807 2012-12-19 13:13:36Z isaid $
/// \brief This file contains the implementation of the source similator.
///

#include "source.h"
#include "math.h"
#include <stdio.h>

void source_ricker_wavelet(DATATYPE* source, DATATYPE dt, int time_steps, DATATYPE fmax) {
    int      it;
    DATATYPE sigma, tau, t, scale;
    FILE* source_file = fopen("data/source.dat", "wb");
    sigma = 0.6*fmax;
    tau   = 1.0;
    scale = 8.0;
    for(it=0; it < time_steps; ++it) {
        t = dt*(DATATYPE)(it-1);
        source[it] = - 2*scale*sigma*
            (sigma-2*sigma* scale*(sigma*t-tau)*(sigma*t-tau))*
            exp(-scale*(sigma*t-tau)*(sigma*t-tau));
        fprintf(source_file, "%d %f\n", it, source[it]);
    }
    fclose(source_file);
}
