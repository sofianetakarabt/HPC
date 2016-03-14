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
/// \file memory.c
/// \version $Id: memory.c 1855 2012-12-29 00:45:28Z isaid $
/// \brief This file contains the implementatins of the memory manipulation routines
///

#include "memory.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

DATATYPE* create_buffer(size_t size) {
    DATATYPE *buffer;
    int status = posix_memalign((void**)&buffer, 4096, size*sizeof(DATATYPE));
    if(!(buffer) || status) {
	    fprintf(stderr,"failed to allocate heap memory\n");
	    exit(EXIT_FAILURE);
	}
    return buffer;
}

void delete_buffer(DATATYPE* buffer) {
    free(buffer);
}

void zero_buffer(DATATYPE* buffer, size_t size) {
        memset(buffer, 0, size*sizeof(DATATYPE));
}

void init_buffer(DATATYPE* buffer,
                 size_t size, DATATYPE value) {
    int i;
#pragma omp parallel for shared(buffer) private(i)
    for(i=0; i<size; i++)
        *(buffer+i) = value;
}

void rand_buffer(DATATYPE* buffer, size_t size) {
    int i;
    srand(time(NULL)); 
#pragma omp parallel for shared(buffer) private(i)
    for(i=0; i<size; i++)
        *(buffer+i) = rand() % 10 + 1;
}

void copy_buffer(DATATYPE* dest_buffer,
                 DATATYPE* src_buffer,
                 size_t size) {
    int i;
#pragma omp parallel for shared(dest_buffer, src_buffer) private(i)
    for(i=0; i<size; i++)
        *(dest_buffer+i)=*(src_buffer+i);
}

int check_buffer(DATATYPE* buffer,
                 DATATYPE* reference, size_t size, DATATYPE e) {
    size_t i, not_valid=0;
    for(i = 0; i < size; i++) {
//    printf(" %lf\t%lf\n",buffer[i],reference[i]);
        if(fabs(buffer[i] - reference[i])>e){
            not_valid ++;
	if (not_valid < 3)printf("err: %d -> %3.45f != %3.45f\n",i,buffer[i],reference[i]);
	}
    }
    return not_valid;
}
