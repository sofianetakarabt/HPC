#ifndef _MEMORY_H_
#define _MEMORY_H_
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
/// \file memory.h
/// \version $Id: memory.h 2767 2015-10-26 08:02:36Z isaid $
/// \brief Routines for memory manipulation.
///
/// A set of routines that help allocating, manipulating and releasing
/// memory buffers.
///

#include <stdlib.h>
#include "config.h"

/// \brief Creates a memory buffer
/// \param size is the number of items that the buffer should contain
/// \return a pointer to the created buffer
DATATYPE* create_buffer(size_t size);

/// \brief Deletes a memory buffer
/// \param buffer is a pointer to the buffer to be deleted
void delete_buffer(DATATYPE * buffer);

/// \brief Sets all the buffer entries to zero
/// \param buffer is a pointer to the buffer to be nullified
/// \param size is the number of the buffer elements
void zero_buffer(DATATYPE* buffer, size_t size);

/// \brief Sets all the buffer entries to given value
/// \param buffer is a pointer to the buffer to be initialized
/// \param size is the number of the buffer elements
/// \param value is the value to be duplicated along the buffer
void init_buffer(DATATYPE* buffer, size_t size, DATATYPE value);

/// \brief Randomly initializes the buffer 
/// \param buffer is a pointer to the buffer to be initialized
/// \param size is the number of the buffer elements
void rand_buffer(DATATYPE* buffer, size_t size);

/// \brief Copies a source buffer to a destination buffer
/// \param dest_buffer is the destination buffer
/// \param src_buffer is the source buffer
/// \param size is the number of the source buffer elements
void copy_buffer(DATATYPE* dest_buffer,
                 DATATYPE* src_buffer, size_t size);

/// \brief Checks if tow buffers have the same values
/// \param buffer is a pointer to the first buffer
/// \param reference is a pointer to the second buffer
/// \param size is the number of the buffers elements
/// \param e is the floating point error margin
int check_buffer(DATATYPE* buffer,
                 DATATYPE* reference, size_t size, DATATYPE e);

#endif //  _MEMORY_H_
