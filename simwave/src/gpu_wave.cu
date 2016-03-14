#include "mpi.h"
#include "wave.h"
#include "memory.h"
#include "math.h"
#include "timer.h"
#include "velocity.h"
#include "source.h"
#include <stdio.h>
#include <cuda.h>

#define U0(z,y,x)   (w->u0[(x+w->sx) + (2*w->sx + w->dimx) * \
                           ((2*w->sy + w->dimy) * (z+w->sz) + (y+w->sy))])
#define U1(z,y,x)   (w->u1[(x+w->sx) + (2*w->sx + w->dimx) * \
                           ((2*w->sy + w->dimy) * (z+w->sz) + (y+w->sy))])


#define ROC2(z,y,x) (w->roc2[x + (w->dimx * (w->dimy*(z) + y))])
#define PHI(z,y,x)  ( w->phi[x + (w->dimx * (w->dimy*(z) + y))])
#define ETA(z,y,x)  (w->eta[(2+w->dimx)*((2+w->dimy)*(z+1) + (y+1)) + (x+1)])


#define PWAVE_COMPUTE_LAPLACIAN() \
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

/*
__global__ void gpu_pwave_swap_pointers(pwave_t* w) {
    DATATYPE* tmp = w->u0;
    w->u0 = w->u1;
    w->u1 = tmp;
}
*/

INLINE void pwave_snapshot(pwave_t* w) {
    unsigned int x, y;
    size_t size = sizeof(DATATYPE);
    for (x = 0; x < w->dimx; x++)
        for (y = 0; y < w->dimy; y++)
            fwrite(&U0(w->sourcez, y, x), size, 1, w->snapshot_file);
}
/*
int checkcpy(pwave_t *w,pwave_t *gw){
	if (!(gw->dimx == w->dimx &&
	gw->dimy == w->dimy &&
	gw->dimz == w->dimz &&
	gw->dx == w->dx &&
	gw->dy == w->dy &&
	gw->dz == w->dz &&
	gw->sx == w->sx &&
	gw->sy == w->sy &&
	gw->sz == w->sz &&
	gw->pmlx == w->pmlx &&
	gw->pmly == w->pmly &&
	gw->pmlz == w->pmlz &&
	gw->sourcex == w->sourcex &&
	gw->sourcey == w->sourcey &&
	gw->sourcez == w->sourcez &&
	gw->time_steps == w->time_steps &&

	gw->coef0 == w->coef0 &&
	gw->hdx2 == w->hdx2 &&
	gw->hdy2 == w->hdy2 &&
	gw->hdz2 == w->hdz2 &&
	gw->dt == w->dt &&
	gw->lambda == w->lambda)) return -1;
	size_t i;
	for(i=0;i<w->size;i++){
		if(w->u0[i]!=gw->u0[i] || w->u1[i]!=gw->u1[i] ) return -1;
	}

	for(i=0;i< w->dimx * w->dimy * w->dimz;i++){
		if(w->phi[i]!=gw->phi[i] || w->roc2[i]!=gw->roc2[i] ) return -1;
	}
	return 1;
}
*/

pwave_t* myCudaMemcpy(pwave_t *w, pwave_t **hw){

	pwave_t * gpuw;
	pwave_t * hd = (pwave_t*)malloc( sizeof(pwave_t));

	cudaMalloc((void**)&gpuw,sizeof(pwave_t));

	hd->size = w->size;
	hd->dimx = w->dimx;
	hd->dimy = w->dimy;
	hd->dimz = w->dimz;
	hd->dx = w->dx;
	hd->dy = w->dy;
	hd->dz = w->dz;
	hd->sx = w->sx;
	hd->sy = w->sy;
	hd->sz = w->sz;
	hd->pmlx = w->pmlx;
	hd->pmly = w->pmly;
	hd->pmlz = w->pmlz;
	hd->sourcex = w->sourcex;
	hd->sourcey = w->sourcey;
	hd->sourcez = w->sourcez;
	hd->time_steps = w->time_steps;

	hd->coef0 = w->coef0;
	hd->hdx2 = w->hdx2;
	hd->hdy2 = w->hdy2;
	hd->hdz2 = w->hdz2;
	hd->dt = w->dt;
	hd->lambda = w->lambda;


/************************* allouer et copier les tableaux internes ***************************************/

	size_t taille = sizeof(DATATYPE)*w->size;
//	cudaError_t error;

	cudaMalloc((void**)&(hd->u0),taille);				cudaMemcpy(hd->u0,w->u0,taille,cudaMemcpyHostToDevice);
	cudaMalloc((void**)&(hd->u1),taille);				cudaMemcpy(hd->u1,w->u1,taille,cudaMemcpyHostToDevice);

	taille = w->dimx * w->dimy * w->dimz * sizeof(DATATYPE);
	cudaMalloc((void**)&(hd->roc2),taille);				cudaMemcpy(hd->roc2,w->roc2,taille,cudaMemcpyHostToDevice);
	cudaMalloc((void**)&(hd->phi),taille);				cudaMemcpy(hd->phi,w->phi,taille,cudaMemcpyHostToDevice);

	taille = w->time_steps*sizeof(DATATYPE);
	cudaMalloc((void**)&(hd->source),taille);			cudaMemcpy(hd->source,w->source,taille,cudaMemcpyHostToDevice);

	taille = (2*w->sx + 1)*sizeof(DATATYPE);
	cudaMalloc((void**)&(hd->coefx),taille);			cudaMemcpy(hd->coefx,w->coefx,taille,cudaMemcpyHostToDevice);
	taille = (2*w->sy + 1)*sizeof(DATATYPE);
	cudaMalloc((void**)&(hd->coefy),taille);			cudaMemcpy(hd->coefy,w->coefy,taille,cudaMemcpyHostToDevice);
	taille = (2*w->sz + 1)*sizeof(DATATYPE);
	cudaMalloc((void**)&(hd->coefz),taille);			cudaMemcpy(hd->coefz,w->coefz,taille,cudaMemcpyHostToDevice);

	taille = (w->dimx+2)*(w->dimy+2)*(w->dimz+2)*sizeof(DATATYPE);
	cudaMalloc((void**)&(hd->eta),taille);				cudaMemcpy(hd->eta,w->eta,taille,cudaMemcpyHostToDevice);

//	hd->snapshot_file = w->snapshot_file;

/************************* copier hd sur le gpu pour recuperer les pointers *****************************/
	cudaMemcpy(gpuw,hd,sizeof(pwave_t),cudaMemcpyHostToDevice);

	(*hw)= hd;
return gpuw;
}



__global__ void gpu_pwave_update_source_swap_pointers(pwave_t* w, unsigned int time_step) {

	U1(w->sourcez, w->sourcey, w->sourcex) = U1(w->sourcez, w->sourcey, w->sourcex)  + w->source[time_step];
	DATATYPE* tmp = w->u0;
	w->u0 = w->u1;
	w->u1 = tmp;
}

__global__ void gpu_pwave_update_fields_3D(pwave_t *w,int t) {


	int x = blockDim.x*blockIdx.x+threadIdx.x;
	int y = blockDim.y*blockIdx.y+threadIdx.y;
	int z = blockDim.z*blockIdx.z+threadIdx.z;
	DATATYPE laplacian;

	if((x < w->dimx) && (y < w->dimy) && (z < w->dimz)){
//		PWAVE_COMPUTE_LAPLACIAN();


		float l1 = w->coef0 * U0(z, y, x);
		float l2 = w->coefx[1]*( U0(z,   y,   x+1) + U0(z,   y,   x-1));
		float l3 = w->coefy[1]*( U0(z,   y+1, x  ) + U0(z,   y-1, x  ));
		float l4 = w->coefz[1]*( U0(z+1, y,   x  ) + U0(z-1, y,   x  ));
		float l5 = w->coefx[2]*( U0(z,   y,   x+2) + U0(z,   y,   x-2));
		float l6 = w->coefy[2]*( U0(z,   y+2, x  ) + U0(z,   y-2, x  ));
		float l7 = w->coefz[2]*( U0(z+2, y,   x  ) + U0(z-2, y,   x  ));
		float l8 = w->coefx[3]*( U0(z,   y,   x+3) + U0(z,   y,   x-3));
		float l9 = w->coefy[3]*( U0(z,   y+3, x  ) + U0(z,   y-3, x  ));
		float l10 = w->coefz[3]*( U0(z+3, y,   x  ) + U0(z-3, y,   x  ));
		float l11 = w->coefx[4]*( U0(z,   y,   x+4) + U0(z,   y,   x-4));
		float l12 = w->coefy[4]*( U0(z,   y+4, x  ) + U0(z,   y-4, x  ));
		float l13 = w->coefz[4]*( U0(z+4, y,   x  ) + U0(z-4, y,   x  ));

	int ind = (x+w->sx) + (2*w->sx + w->dimx) * \
                         ((2*w->sy + w->dimy) * (z+w->sz) + (y+w->sy));
		if( ind == 1 ){
//			printf("%3.45f + %3.45f = %3.45f\n",l1,l2,l1+l2);
			printf("check\n");
		}


	l1 = l1+l2+l3;
	l2 = l4+l5+l6;
	l3 = l7+l8+l9;
	l4 = l10+l11+l12+l13;
	laplacian = l1+l2+l3+l4;

		if((z >= w->pmlz) && (z < w->dimz - w->pmlz) &&
		(y >= w->pmly) && (y < w->dimy - w->pmly) &&
		(x >= w->pmlx) && (x < w->dimx - w->pmlx)) {
			PWAVE_UPDATE_INNER_FIELDS();
		} else {
                   PWAVE_UPDATE_PML_FIELDS();
		}
	}

}


__global__ void gpu_pwave_update_fields_2D(pwave_t *w,int t) {

	DATATYPE laplacian;
	int x = blockDim.x*blockIdx.x+threadIdx.x;
	int y = blockDim.y*blockIdx.y+threadIdx.y;
	int z;


	for (z=0;z<w->dimz;z++){
		if((x < w->dimx) && (y < w->dimy)){
//			PWAVE_COMPUTE_LAPLACIAN();
		float l1 = w->coef0 * U0(z, y, x);
		float l2 = w->coefx[1]*( U0(z,   y,   x+1) + U0(z,   y,   x-1));
		float l3 = w->coefy[1]*( U0(z,   y+1, x  ) + U0(z,   y-1, x  ));
		float l4 = w->coefz[1]*( U0(z+1, y,   x  ) + U0(z-1, y,   x  ));
		float l5 = w->coefx[2]*( U0(z,   y,   x+2) + U0(z,   y,   x-2));
		float l6 = w->coefy[2]*( U0(z,   y+2, x  ) + U0(z,   y-2, x  ));
		float l7 = w->coefz[2]*( U0(z+2, y,   x  ) + U0(z-2, y,   x  ));
		float l8 = w->coefx[3]*( U0(z,   y,   x+3) + U0(z,   y,   x-3));
		float l9 = w->coefy[3]*( U0(z,   y+3, x  ) + U0(z,   y-3, x  ));
		float l10 = w->coefz[3]*( U0(z+3, y,   x  ) + U0(z-3, y,   x  ));
		float l11 = w->coefx[4]*( U0(z,   y,   x+4) + U0(z,   y,   x-4));
		float l12 = w->coefy[4]*( U0(z,   y+4, x  ) + U0(z,   y-4, x  ));
		float l13 = w->coefz[4]*( U0(z+4, y,   x  ) + U0(z-4, y,   x  ));

		int ind = (x+w->sx) + (2*w->sx + w->dimx) * \
                         ((2*w->sy + w->dimy) * (z+w->sz) + (y+w->sy));
		if( ind == 1 ){
//			printf("%3.45f + %3.45f = %3.45f\n",l1,l2,l1+l2);
			printf("check\n");
		}
		l1 = l1+l2+l3;
		l2 = l4+l5+l6;
		l3 = l7+l8+l9;
		l4 = l10+l11+l12+l13;
		laplacian = l1+l2+l3+l4;


		if((z >= w->pmlz) && (z < w->dimz - w->pmlz) &&
			(y >= w->pmly) && (y < w->dimy - w->pmly) &&
			(x >= w->pmlx) && (x < w->dimx - w->pmlx)) {
				PWAVE_UPDATE_INNER_FIELDS();
		} else {
				PWAVE_UPDATE_PML_FIELDS();
			}
		}
	}
}


extern "C"
float pwave_multiple_update_fields(pwave_t *w,int nbThreads_x,int nbThreads_y,int nbThreads_z,bool_t snapshot_enabled,int nb_snap) {


    float time;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);


	int tx, ty, tz;
	tx =  w->dimx;
 	ty =  w->dimy;
 	tz =  w->dimz;

	int taille_grille_x = tx%nbThreads_x ? tx/nbThreads_x+1:tx/nbThreads_x;
	int taille_grille_y = ty%nbThreads_y ? ty/nbThreads_y+1:ty/nbThreads_y;
	int taille_grille_z = tz%nbThreads_z ? tz/nbThreads_z+1:tz/nbThreads_z;

	dim3 threadsParBloc3D(nbThreads_x,nbThreads_y,nbThreads_z);
	dim3 threadsParBloc2D(nbThreads_x,nbThreads_y);
	dim3 tailleGrille3D(taille_grille_x,taille_grille_y,taille_grille_z);
	dim3 tailleGrille2D(taille_grille_x,taille_grille_y);

	/* creer w sur le gpu*/
	pwave_t *hw = (pwave_t*)malloc(sizeof(pwave_t));

	pwave_t *gpuw = myCudaMemcpy(w,&hw);

	int t;

//printf("----------------------------------------- DEBUT KERNEL CUDA -------------------------------\n\n");
//	cudaError_t error;
//    printf("Occupancy calculator elapsed time:  %3.3f ms \n", time);

 	for(t = 0; t < w->time_steps; ++t) {

		gpu_pwave_update_fields_3D<<<tailleGrille3D,threadsParBloc3D>>>(gpuw,t);
//		gpu_pwave_update_fields_2D<<<tailleGrille2D,threadsParBloc2D>>>(gpuw,t);
		gpu_pwave_update_source_swap_pointers<<<1,1>>>(gpuw,t);


		if((snapshot_enabled) && (t%nb_snap == 0)){
			cudaMemcpy(w->u0, hw->u0, sizeof(DATATYPE)*w->size ,cudaMemcpyDeviceToHost);
			pwave_snapshot(w);
		}
	}
	cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);

	//  recopier pour le check
    cudaMemcpy(w->u0, hw->u0, sizeof(DATATYPE) *w->size ,cudaMemcpyDeviceToHost);


//printf("--------------------------------------- FIN KERNEL CUDA <<----ERREUR D ARRONDIS---->>-------------------\n");
/*
	cudaFree(gpuw->u0);
	cudaFree(gpuw->u1);
	cudaFree(gpuw->coefx);
	cudaFree(gpuw->coefy);
	cudaFree(gpuw->coefz);
	cudaFree(gpuw->roc2);
	cudaFree(gpuw->phi);
	cudaFree(gpuw->eta);
	cudaFree(gpuw);
	free(hw);
*/
	return time;
}












