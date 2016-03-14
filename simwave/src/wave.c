#include "wave.h"
#include "memory.h"
#include "math.h"
#include "velocity.h"
#include "source.h"
#include <stdio.h>

#define U0(z,y,x)   (w->u0[(x+w->sx) + (2*w->sx + w->dimx) * \
                           ((2*w->sy + w->dimy) * (z+w->sz) + (y+w->sy))])
#define U1(z,y,x)   (w->u1[(x+w->sx) + (2*w->sx + w->dimx) * \
                           ((2*w->sy + w->dimy) * (z+w->sz) + (y+w->sy))])
#define ROC2(z,y,x) (w->roc2[x + (w->dimx * (w->dimy*(z) + y))])
#define PHI(z,y,x)  ( w->phi[x + (w->dimx * (w->dimy*(z) + y))])
#define ETA(z,y,x)  (w->eta[(2+w->dimx)*((2+w->dimy)*(z+1) + (y+1)) + (x+1)])
#define MAX(a,b)    (a>b?a:b)
#define MIN(a,b)    (a<b?a:b)

#define SETUP_STENCIL_COEFFS(t) \
    t[0] = -205. /  72.;        \
    t[1] =    8. /   5.;        \
    t[2] =   -1. /   5.;        \
    t[3] =    8. / 315.;        \
    t[4] =   -1. / 560.;

#define SCALE_STENCIL_COEFFS(delta, t) \
    t[0] /= pow(delta,2);              \
    t[1] /= pow(delta,2);              \
    t[2] /= pow(delta,2);              \
    t[3] /= pow(delta,2);              \
    t[4] /= pow(delta,2);   

pwave_t* pwave_create(parser* p,int nproc) {
    DATATYPE courant_number;
    unsigned int i, z, y, x;
    int grid[3], du[3], source_loc[3], taper[3];
    pwave_t *w               = (pwave_t*) malloc(sizeof(pwave_t));
    DATATYPE cfl            = parser_get_float(p, "cfl");
    DATATYPE max_frequency  = parser_get_float(p, "fmax");
    DATATYPE vmax           = parser_get_float(p, "vmax");
    DATATYPE vmin           = parser_get_float(p, "vmin");
    parser_get_vect_int(p, "grid",       grid,       3);
    parser_get_vect_int(p, "dgrid",      du,         3);
    parser_get_vect_int(p, "source_loc", source_loc, 3);
    parser_get_vect_int(p, "taper",      taper,      3);

    w->dimx = grid[0];
    w->dimy = grid[1];
    w->dimz = grid[2]/nproc;

    w->sx = w->sy = w->sz = 4;

    w->sourcex = source_loc[0] == -1 ? w->dimx/4 : source_loc[0];
    w->sourcey = source_loc[1] == -1 ? w->dimy/4 : source_loc[1];
    w->sourcez = source_loc[2] == -1 ? w->dimz/2 : source_loc[2];

    w->size   = (2*w->sx + w->dimx)*(2*w->sy + w->dimy)*(2*w->sz + w->dimz);//*sizeof(float);
    w->u0     = create_buffer(w->size);
    w->u1     = create_buffer(w->size);
    w->roc2   = create_buffer(w->dimx * w->dimy * w->dimz);
    w->coefx  = create_buffer(2*w->sx + 1);
    w->coefy  = create_buffer(2*w->sy + 1);
    w->coefz  = create_buffer(2*w->sz + 1);

    w->eta    = create_buffer((w->dimx+2)*(w->dimy+2)*(w->dimz+2));
    w->phi    = create_buffer(w->dimx * w->dimy * w->dimz);

    zero_buffer(w->u0, w->size);
    zero_buffer(w->u1, w->size);
    zero_buffer(w->phi, w->dimx * w->dimy * w->dimz);
    zero_buffer(w->eta, (w->dimx+2)*(w->dimy+2)*(w->dimz+2));
  
    w->lambda = ((unsigned int)(vmax/max_frequency));
    w->dx     = MAX(w->lambda/16, du[0]);
    w->dy     = MAX(w->lambda/16, du[1]);
    w->dz     = MAX(w->lambda/16, du[2]);
    
    w->hdx2   = 1. / (4. * pow(w->dx, 2));
    w->hdy2   = 1. / (4. * pow(w->dy, 2));
    w->hdz2   = 1. / (4. * pow(w->dz, 2));
    
    w->pmlx = taper[0]*w->lambda/w->dx;
    w->pmly = taper[1]*w->lambda/w->dy;
    w->pmlz = taper[2]*w->lambda/w->dz;
        
    SETUP_STENCIL_COEFFS(w->coefx);
    SETUP_STENCIL_COEFFS(w->coefy);
    SETUP_STENCIL_COEFFS(w->coefz);
    
    SCALE_STENCIL_COEFFS(w->dx, w->coefx);
    SCALE_STENCIL_COEFFS(w->dy, w->coefy);
    SCALE_STENCIL_COEFFS(w->dz, w->coefz);
    w->coef0 = w->coefx[0] + w->coefy[0] + w->coefz[0];

    courant_number = fabs(w->coefx[0]) + fabs(w->coefy[0]) + fabs(w->coefz[0]);
    for(i=1; i<w->sx+1; i++)
        courant_number += 2.*fabs(w->coefx[i]);
    for(i=1; i<w->sy+1; i++) 
        courant_number += 2.*fabs(w->coefy[i]);
    for(i=1; i<w->sz+1; i++)
        courant_number += 2.*fabs(w->coefz[i]);   
    courant_number = 2. / sqrt(courant_number);

    velocity_generate_model(w->roc2, w->dimx, w->dimy, w->dimz, vmin, vmax, 2);

    w->dt = cfl * courant_number / velocity_get_max(w->roc2, w->dimx, w->dimy, w->dimz);
    w->time_steps = parser_get_float(p, "iter");
    
    for (z=0; z < w->dimz; z++)
        for (y=0; y < w->dimy; y++)
            for (x=0; x < w->dimx; x++)
                ROC2(z,y,x) = pow(w->dt, 2) * pow(ROC2(z,y,x), 2);

    pwave_compute_pml(w, p);

    w->source = create_buffer(w->time_steps); 
    source_ricker_wavelet(w->source, w->dt, w->time_steps, max_frequency);

    if(!parser_get_bool(p, "nosnap"))
        w->snapshot_file = fopen(parser_get_string(p, "output"), "wb");

    pwave_make_display_cmd(w, p);
    return w;
}

void pwave_delete(pwave_t* w, parser *p) {
    delete_buffer(w->u0);
    delete_buffer(w->u1);
    delete_buffer(w->roc2);
    delete_buffer(w->eta);
    delete_buffer(w->phi); 
    delete_buffer(w->coefx);
    delete_buffer(w->coefy);
    delete_buffer(w->coefz);
    delete_buffer(w->source);
    if(!parser_get_bool(p, "nosnap")) fclose(w->snapshot_file);
    free(w);
}

void pwave_print(pwave_t* w) {
    int i;
    fprintf(stdout, "\n");
    fprintf(stdout, "... Wave information:\n");
    fprintf(stdout, "... domain size      = %u x %u x %u (%f MB)\n",
            w->dimx, w->dimy, w->dimz,
            w->size/1024./1024.);
    fprintf(stdout, "... Lambda           = %f\n", w->lambda);
    fprintf(stdout, "... dS               = %u x %u x %u\n",
            w->dx, w->dy, w->dz); 
    fprintf(stdout, "... dt               = %g\n", w->dt);
    fprintf(stdout, "... stencil size     = %u x %u x %u\n",
            2*w->sx+1, 2*w->sy+1, 2*w->sz+1);
    fprintf(stdout, "... pml damping      = %u x %u x %u\n",
            w->pmlx, w->pmly, w->pmlz);
    fprintf(stdout, "... source location  = %u x %u x %u\n",
            w->sourcex, w->sourcey, w->sourcez);
    fprintf(stdout, "... coefx            = ");
    for(i=w->sx; i>0; i--)
        fprintf(stdout, "%g ", w->coefx[i]);
    fprintf(stdout, "%g ", w->coefx[0]);
    for(i=1; i<w->sx+1; i++)
        fprintf(stdout, "%g ", w->coefx[i]); 
    fprintf(stdout, "\n... coefy            = ");
    for(i=w->sy; i>0; i--)
        fprintf(stdout, "%g ", w->coefy[i]);
    fprintf(stdout, "%g ", w->coefy[0]);
    for(i=1; i<w->sy+1; i++)
        fprintf(stdout, "%g ", w->coefy[i]); 
    fprintf(stdout, "\n... coefz            = ");
    for(i=w->sz; i>0; i--)
        fprintf(stdout, "%g ", w->coefz[i]);
    fprintf(stdout, "%g ", w->coefz[0]);
    for(i=1; i<w->sz+1; i++)
        fprintf(stdout, "%g ", w->coefz[i]);
    fprintf(stdout, "\n... hdx2, hdy2, hdz2 = %f, %f, %f\n",
            w->hdx2, w->hdy2, w->hdz2);
    fprintf(stdout, "... time steps       = %u\n", w->time_steps);
    fprintf(stdout, "\n");
}

#define PWAVE_COMPUTE_LAPLACIAN()                               \
    laplacian = w->coef0 * U0(z, y, x)                          \
        + w->coefx[1]*( U0(z,   y,   x+1) + U0(z,   y,   x-1))  \
        + w->coefy[1]*( U0(z,   y+1, x  ) + U0(z,   y-1, x  ))  \
        + w->coefz[1]*( U0(z+1, y,   x  ) + U0(z-1, y,   x  ))  \
        + w->coefx[2]*( U0(z,   y,   x+2) + U0(z,   y,   x-2))  \
        + w->coefy[2]*( U0(z,   y+2, x  ) + U0(z,   y-2, x  ))  \
        + w->coefz[2]*( U0(z+2, y,   x  ) + U0(z-2, y,   x  ))  \
        + w->coefx[3]*( U0(z,   y,   x+3) + U0(z,   y,   x-3))  \
        + w->coefy[3]*( U0(z,   y+3, x  ) + U0(z,   y-3, x  ))  \
        + w->coefz[3]*( U0(z+3, y,   x  ) + U0(z-3, y,   x  ))  \
        + w->coefx[4]*( U0(z,   y,   x+4) + U0(z,   y,   x-4))  \
        + w->coefy[4]*( U0(z,   y+4, x  ) + U0(z,   y-4, x  ))  \
        + w->coefz[4]*( U0(z+4, y,   x  ) + U0(z-4, y,   x  ));

#define PWAVE_UPDATE_INNER_FIELDS()                                  \
    U1(z,y,x) = 2.*U0(z,y,x) - U1(z,y,x) + ROC2(z,y,x) * laplacian;

#define PWAVE_UPDATE_PML_FIELDS()                                                        \
    U1(z,y,x) = ((2.-ETA(z,y,x)*ETA(z,y,x) + 2.*ETA(z,y,x))*U0(z,y,x)                    \
                 - U1(z,y,x) + ROC2(z,y,x)*(laplacian + PHI(z,y,x)))/(1.+2.*ETA(z,y,x)); \
    PHI(z,y,x)= (PHI(z,y,x)-                                                             \
                 (( ETA(z,   y,   x+1) - ETA(z,   y,   x-1))                             \
                  *( U0(z,   y,   x+1) -  U0(z,   y,   x-1))*w->hdx2                     \
                  +(ETA(z,   y+1, x  ) - ETA(z,   y-1, x  ))                             \
                  *( U0(z,   y+1, x  ) -  U0(z,   y-1, x  ))*w->hdy2                     \
                  +(ETA(z+1, y,   x  ) - ETA(z-1, y,   x  ))                             \
                  *( U0(z+1, y,   x  ) -  U0(z-1, y,   x  ))*w->hdz2))/(1.+ETA(z,y,x));



INLINE void pwave_update_fields(pwave_t *w) {
    unsigned int z, y, x;
    DATATYPE laplacian;
//    printf("tat");
//#pragma omp parallel for private(laplacian, z,y,x)
    for(z = 0; z < w->dimz ; z++) {
        for(y = 0; y < w->dimy; y++) {
            for(x = 0; x < w->dimx; x++) {
//                PWAVE_COMPUTE_LAPLACIAN();
		float l1,l2,l3,l4,l5;
		l1 = w->coef0 * U0(z, y, x);
		l2 = w->coefx[1]*( U0(z,   y,   x+1) + U0(z,   y,   x-1));
		l3 = w->coefy[1]*( U0(z,   y+1, x  ) + U0(z,   y-1, x  )); 
		l4 = w->coefz[1]*( U0(z+1, y,   x  ) + U0(z-1, y,   x  ));
		l5 = w->coefx[2]*( U0(z,   y,   x+2) + U0(z,   y,   x-2));
		float l6 = w->coefy[2]*( U0(z,   y+2, x  ) + U0(z,   y-2, x  ));
        float l7 = w->coefz[2]*( U0(z+2, y,   x  ) + U0(z-2, y,   x  ));
		float l8 = w->coefx[3]*( U0(z,   y,   x+3) + U0(z,   y,   x-3));
        float l9 = w->coefy[3]*( U0(z,   y+3, x  ) + U0(z,   y-3, x  ));
        float l10 = w->coefz[3]*( U0(z+3, y,   x  ) + U0(z-3, y,   x  ));
        float l11 = w->coefx[4]*( U0(z,   y,   x+4) + U0(z,   y,   x-4));
        float l12 = w->coefy[4]*( U0(z,   y+4, x  ) + U0(z,   y-4, x  ));
        float l13 = w->coefz[4]*( U0(z+4, y,   x  ) + U0(z-4, y,   x  ));

	l1 = l1+l2+l3;
	l2 = l4+l5+l6;
	l3 = l7+l8+l9;
	l4 = l10+l11+l12+l13;
		laplacian = l1+l2+l3+l4;
/*		if(x==25 && y == 25 && z == 44){
			printf("lap g == %3.45f\n",laplacian);
			printf("------------- l12 = %3.45f + %3.45f\n +%3.45f= %3.45f\n------------------\n",l1,l2,l3,l1+l2+l3);

		}
*/               if((z >= w->pmlz) && (z < w->dimz - w->pmlz) &&
                   (y >= w->pmly) && (y < w->dimy - w->pmly) && 
                   (x >= w->pmlx) && (x < w->dimx - w->pmlx)) {
                    PWAVE_UPDATE_INNER_FIELDS();
                } else {
                   PWAVE_UPDATE_PML_FIELDS();
                }
            }
        }
    }
}



/**fonction orig + optimisation**/
/*
INLINE void pwave_update_fields(pwave_t *w) {
    unsigned int z, y, x,zpml,ypml,xpml;
    DATATYPE laplacian;

//	#pragma omp parallel
//	#pragma omp for private(y,x)

	xpml=w->dimx - w->pmlx ;ypml=w->dimy - w->pmly ;zpml = w->dimz - w->pmlz ;
#pragma omp parallel for private(laplacian, y,x)
    for(z = w->pmlz; z < zpml ; z++) {
        for(y = w->pmly; y < ypml; y++) {
            for(x = w->pmlx; x < xpml; x++) {
                PWAVE_COMPUTE_LAPLACIAN();
                PWAVE_UPDATE_INNER_FIELDS();
            }
        }
    }
#pragma omp parallel for private(laplacian, y,x)
    for(z = 0; z < w->pmlz ; z++) {
        for(y = 0; y < w->dimy; y++) {
            for(x = 0; x < w->dimx; x++) {
                PWAVE_COMPUTE_LAPLACIAN();
                PWAVE_UPDATE_PML_FIELDS();
            }
        }
    }
#pragma omp parallel for private(laplacian, y,x)
    for(z = zpml; z < w->dimz ; z++) {
        for(y = 0; y < w->dimy; y++) {
            for(x = 0; x < w->dimx; x++) {
                PWAVE_COMPUTE_LAPLACIAN();
                PWAVE_UPDATE_PML_FIELDS();
            }
        }
    }
#pragma omp parallel for private(laplacian, y,x)
    for(z = w->pmlz; z < zpml ; z++) {
        for(y = 0; y < w->pmly; y++) {
            for(x = 0; x < w->dimx; x++) {
                PWAVE_COMPUTE_LAPLACIAN();
                PWAVE_UPDATE_PML_FIELDS();
            }
        }
    }
#pragma omp parallel for private(laplacian, y,x)
    for(z = w->pmlz; z < zpml ; z++) {
        for(y = ypml; y < w->dimy; y++) {
            for(x = 0; x < w->dimx; x++) {
                PWAVE_COMPUTE_LAPLACIAN();
                PWAVE_UPDATE_PML_FIELDS();
            }
        }
    }
#pragma omp parallel for private(laplacian, y,x)
    for(z = w->pmlz; z < zpml ; z++) {
        for(y = w->pmly; y < ypml; y++) {
            for(x = 0; x < w->pmlx; x++) {
                PWAVE_COMPUTE_LAPLACIAN();
                PWAVE_UPDATE_PML_FIELDS();
            }
        }
    }
#pragma omp parallel for private(laplacian, y,x)
    for(z = w->pmlz; z < zpml ; z++) {
        for(y = w->pmly; y < ypml; y++) {
            for(x = xpml; x < w->dimx; x++) {
                PWAVE_COMPUTE_LAPLACIAN();
                PWAVE_UPDATE_PML_FIELDS();
            }
        }
    }
}
*/
INLINE void pwave_update_source(pwave_t* w, unsigned int time_step) {
    U1(w->sourcez, w->sourcey, w->sourcex) += w->source[time_step];
}

INLINE void pwave_swap_pointers(pwave_t* w) {
    DATATYPE* tmp = w->u0;
    w->u0 = w->u1;
    w->u1 = tmp;
}

INLINE void pwave_snapshot(pwave_t* w) {
    unsigned int x, y;
    size_t size = sizeof(DATATYPE);
    for (x = 0; x < w->dimx; x++)
        for (y = 0; y < w->dimy; y++)
            fwrite(&U0(w->sourcez, y, x), size, 1, w->snapshot_file);
}

#define PML_SCALE_FACTOR(taper, f) 3.*f*log(1000.)/(2.*taper);

INLINE void pwave_compute_pml(pwave_t *w, parser *p) {
    int x, y, z;
    DATATYPE vx, vy, vz;
    DATATYPE scalex, scaley, scalez;
    DATATYPE fmax = parser_get_float(p, "fmax");
    int taper[3];
    parser_get_vect_int(p, "taper", taper, 3);
    scalex = PML_SCALE_FACTOR(taper[0], fmax);
    scaley = PML_SCALE_FACTOR(taper[1], fmax);
    scalez = PML_SCALE_FACTOR(taper[2], fmax);
    for(z = 0; z < w->dimz+2; ++z) {
        for(y = 0; y < w->dimy+2; ++y) {
            for(x = 0; x < w->dimx+2; ++x) {
                if(z < w->pmlz)
                    vz = pow(1.*(w->pmlz - z)/w->pmlz, 2)*scalez;
                else if(z >= ((w->dimz+2) - w->pmlz))
                    vz = pow(1.*(z - ((w->dimz+2) - w->pmlz - 1))/w->pmlz, 2)*scalez;
                else
                    vz = 0.;
                if(y < w->pmly)
                    vy = pow(1.*(w->pmly - y)/w->pmly, 2)*scaley;
                else if(y >= ((w->dimy+2) - w->pmly))
                    vy = pow(1.*(y - ((w->dimy+2) - w->pmly - 1))/w->pmly, 2)*scaley;
                else
                    vy = 0.;
                if(x < w->pmlx)
                    vx = pow(1.*(w->pmlx - x)/w->pmlx, 2)*scalex;
                else if(x >= ((w->dimx+2) - w->pmlx))
                    vx = pow(1.*(x - ((w->dimx+2) - w->pmlx - 1))/w->pmlx, 2)*scalex;
                else
                    vx = 0.;
                ETA(z-1,y-1,x-1) = MAX(vx, MAX(vy, vz))*w->dt;
            }
        }
    }
}

void pwave_make_display_cmd(pwave_t* w, parser *p) {
    FILE* f=fopen("make/display.mk", "w");
    float ratio = (float)w->dimy/(float)w->dimx;
    if (f == NULL) {
        fprintf(stderr,"Can not open output file\n");
        exit(EXIT_FAILURE);
    }
    fprintf(f, "display:\n");
    fprintf(f, "ifndef VISU_BINDIR\n");
    fprintf(f, "\t${VERBOSE}echo please set VISU_BINDIR env variable\n");
    fprintf(f, "else\n");
    fprintf(f, "\t${VERBOSE}${VISU} < %s "
            "n1=%d n2=%d width=%d height=%d clip=-1 "
            "title=\"wave simulator\"&\n",
            parser_get_string(p, "output"),
            w->dimy, w->dimx, 768, (int)(ratio*768));
    fprintf(f, "endif\n");
    fclose(f);
}
