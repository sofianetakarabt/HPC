#ifndef _WAVE_H_
#define _WAVE_H_
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
/// \file include/wave.h
/// \version $Id: wave.h 2783 2015-10-26 09:43:28Z isaid $
/// \brief A \em reflectionless wave (\b pwave) handler on the CPU.
///
/// Contains all the routines necessary for a wave simulation with boundary
/// conditions (using the <b>Perfectly Matched Layer</b> method) on the CPU side.
/// We call the new wave \b pwave for sake of clarity. 
///

#include "config.h"
#include "parser.h"
#include <stdio.h>

/// \brief Typedef \ref _pwave for convenience
typedef struct _pwave pwave_t;

/// \brief A descriptor for a reflectionless wave (pwave) for the CPU.
/// i.e. The wave is absorbed on the domain edges
struct _pwave {
    /// The domain width
    unsigned int dimx;
    
    /// The domain height
    unsigned int dimy;
    
    /// The domain depth
    unsigned int dimz;
    
    /// The size of a grid step along the \b x axis
    unsigned int dx;
    
    /// The size of a grid step along the \b y axis
    unsigned int dy;

    /// The size of a grid step along the \b z axis
    unsigned int dz;
    
    /// The half size of the stencil along the \b x axis
    unsigned int sx;

    /// The half size of the stencil along the \b y axis 
    unsigned int sy;

    /// The half size of the stencil along the \b z axis
    unsigned int sz;

    /// The length of the PML layer on the \b x axis
    unsigned int pmlx;
    
    /// The length of the PML layer on the \b y axis
    unsigned int pmly;

    /// The length of the PML layer on the \b z axis
    unsigned int pmlz;
   
    /// The position of the source in the \b x axis
    unsigned int sourcex;

    /// The position of the source in the \b y axis
    unsigned int sourcey;

    /// The position of the source in the \b z axis
    unsigned int sourcez;

    /// The number of the simulation time steps
    unsigned int time_steps;
    
    /// The size of the fields in number of elements
    size_t size;

    /// Contains the fields pressure value at time step \em t 
    DATATYPE* u0;
    
    /// Contains the fields pressure value at time step \em t+1 
    DATATYPE* u1;

    /// Contains the velocity values of the traversed mediums 
    DATATYPE* roc2;

    /// Contains the terms of the source
    DATATYPE* source;

    /// Contains the stencil coefficients along the \b x axis 
    DATATYPE* coefx;

    /// Contains the stencil coefficients along the \b y axis 
    DATATYPE* coefy;

    /// Contains the stencil coefficients along the \b x axis 
    DATATYPE* coefz;

    /// Contains the PML initial terms 
    DATATYPE* eta;

    /// Contains the altered PML terms 
    DATATYPE* phi;

    /// The coefficient of the stencil central point
    DATATYPE  coef0;

    /// The PML absorbing factor on the \b x axis
    DATATYPE hdx2;

    /// The PML absorbing factor on the \b y axis
    DATATYPE hdy2;

    /// The PML absorbing factor on the \b z axis
    DATATYPE hdz2;

    /// The duration of a time step
    DATATYPE  dt;

    /// The minimum points number per wavelength
    DATATYPE  lambda;

    /// The descriptor for the snapshots file (produced by the CPU)
    FILE* snapshot_file;
};

/// \brief Creates a pwave descriptor for the CPU
/// \param p is a pointer to an option parser 
/// \return a pointer to the created pwave 
pwave_t *pwave_create(parser *p,int nproc);

/// \brief Deletes a pwave descriptor
/// \param w is a pointer to the pwave descriptor to be deleted 
/// \param p is a pointer to an option parser 
void pwave_delete(pwave_t *w, parser *p);

/// \brief Prints informations about a pwave descriptor
/// \param w is a pointer to the pwave descriptor to be printed 
void pwave_print(pwave_t *w);

/// \brief Computes the pressure on each element of the domain
/// \param w is a pointer to the pwave descriptor to be updated
void pwave_update_fields(pwave_t *w);
void pwave_update_fields_mp(pwave_t *w);
/// \brief Adds an impulse to the grid point situated at the
/// source location 
/// \param w is a pointer to the pwave descriptor to be updated 
/// \param time_step is the current time step
void pwave_update_source(pwave_t *w, unsigned int time_step);

/// \brief Switches the \ref wave_t::u0 and \ref wave_t::u1 pointers 
/// \param w is a pointer to the pwave descriptor to be updated 
void pwave_swap_pointers(pwave_t *w);

/// \brief Saves a snapshot of the pwave fields
/// \param w is a pointer to the pwave descriptor to be saved 
void pwave_snapshot(pwave_t* w);

/// \brief Computes the PML initial terms
/// \param w is a pointer to the pwave descriptor to be updated
/// \param p is a pointer to an option parser 
void pwave_compute_pml(pwave_t* w, parser *p);

/// \brief Generates a Makefile rule to visualize the pwave propagation
/// \see make display 
/// \param w is a pointer to the pwave descriptor to be saved
/// \param p is a pointer to an option parser 
void pwave_make_display_cmd(pwave_t* w, parser *p);

#endif // _WAVE_H_
