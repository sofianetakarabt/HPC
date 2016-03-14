#include "config.h"
#include "parser.h"
#include "timer.h"
#include "memory.h"
#include "wave.h"
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

#define SIZE_H_N 100
#define MPI_Real MPI_FLOAT



#define IND(z,y,x)   ((x+w->sx) + (2*w->sx + w->dimx) *((2*w->sy + w->dimy) * (z+w->sz) + (y+w->sy)))

extern /*DATATYPE*/ float pwave_multiple_update_fields(pwave_t *w,int nbThreads_x,int nbThreads_y,int nbThreads_z,bool_t snapshot_enabled,int nb_snap);
extern float pwave_multiple_update_fields_opt(pwave_t *w,int nbThreads_x,int nbThreads_y,int nbThreads_z,bool_t snapshot_enabled,int nb_snap);

/*
pwave_t* p_create(pwave_t *w, int nproc){

	pwave_t * hd = (pwave_t*)malloc( sizeof(pwave_t));



	hd->size = (2*w->sx + w->dimx)*(2*w->sy + w->dimy)*(2*w->sz + w->dimz/nproc);
	hd->dimx = w->dimx;
	hd->dimy = w->dimy;
	hd->dimz = w->dimz/nproc;
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
/*
	size_t taille = sizeof(DATATYPE)*w->size;

	hd->u0 = calloc(DATATYPE,taille);
	hd->u1 = calloc(DATATYPE,taille);

	taille = w->dimx * w->dimy * w->dimz * sizeof(DATATYPE);
	hd->roc2 = calloc(DATATYPE,taille);
	hd->phi = calloc(DATATYPE,taille);


	taille = w->time_steps*sizeof(DATATYPE);
	hd->source = calloc(DATATYPE,taille);


	taille = (2*w->sx + 1)*sizeof(DATATYPE);
	hd->coefx = calloc(DATATYPE,taille);
	taille = (2*w->sy + 1)*sizeof(DATATYPE);
	hd->coefy = calloc(DATATYPE,taille);
	taille = (2*w->sz + 1)*sizeof(DATATYPE);
	hd->coefz = calloc(DATATYPE,taille);

	taille = (w->dimx+2)*(w->dimy+2)*(w->dimz+2)*sizeof(DATATYPE);
	hd->eta = calloc(DATATYPE,taille);


//	hd->snapshot_file = w->snapshot_file;


return hd;
}
*/

int main(int argc, char* argv[]) {


	printf("BEGIN\n");
	int rank,nproc;
	MPI_Status status;
	MPI_Request request;
	char hostname[SIZE_H_N];
	gethostname(hostname,SIZE_H_N);
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &nproc ); 


    timing_unit_t unit;
    pwave_t *w;
    int t, nb_snap;
    bool_t verbose, snapshot_enabled, check,opt;
    double time = 0.0,epsilon;
	float time1 = 0.0;
    parser *p = parser_create("a synthetic pwave simulator");

    PARSER_BOOTSTRAP(p);
    PARSER_ADD_PML_OPTIONS(p);

    parser_parse(p, argc, argv);

    unit =
        parser_get_bool(p, "ms") ? MILLISECONDS :
        parser_get_bool(p, "us") ? MICROSECONDS :
        parser_get_bool(p, "ns") ? NANOSECONDS  : SECONDS;

    snapshot_enabled = ! parser_get_bool(p, "nosnap");
    nb_snap          = parser_get_int(p, "nbsnap");
    verbose          = parser_get_bool(p, "verbose");

    check            = parser_get_bool(p, "check");
    epsilon          = parser_get_float(p, "epsilon");
	
	opt 			 = !parser_get_bool(p, "opt");

	int bloc[3];
//	int grid[3];
    parser_get_vect_int(p, "bloc", bloc, 3);
//    parser_get_vect_int(p, "grid", grid, 3);


/********** distrib data***********/



 printf("start simulation \n");

if (rank == 0){
	w = pwave_create(p,1);
//    wg = pwave_create(p);
	}
else{
	w = pwave_create(p, nproc);
}

int wsize ;//= (8 + grid[2])*(8 + grid[2])*(8 + grid[2]) ;
if (!rank) wsize= w->size;
//MPI_Bcast (&nproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast (&wsize, 1, MPI_INT, 0, MPI_COMM_WORLD);

//printf("\n ------------	START SIMULATION	------------------ \n");


/*
if (opt){
	printf("optimal mode is OFF\n");
	start_timer();
	time1 = pwave_multiple_update_fields(wg,bloc[0],bloc[1],bloc[2],snapshot_enabled,nb_snap);
}else{
	printf("optimal mode is ON\n");
	start_timer();
	time1 = pwave_multiple_update_fields_opt(wg,bloc[0],bloc[1],bloc[2],snapshot_enabled,nb_snap);
}
*/

double time_cpu = 0.0;



 printf("start 1 \n");                         
//if(check){

printf("EXEC\n");

if (!rank) { init_timer();start_timer();}

	for(t = 0; t < w->time_steps; ++t) {
		pwave_update_fields(w);
        pwave_update_source(w, t);
        pwave_swap_pointers(w);
	printf("time step %d\n",t);
	if (nproc > 1 ){
		if(rank > 0 && rank < nproc-1){

				MPI_Isend(w->u0 + IND(0,0,0),(2*w->sx+w->dimx)*(2*w->sy+w->dimy),MPI_Real,rank-1,0,MPI_COMM_WORLD,&request); /////////////// 4*
				MPI_Isend(w->u0 + IND(w->dimz-4,w->dimy-4,w->dimz-4),(2*w->sx+w->dimx)*(2*w->sy+w->dimy),MPI_Real,rank+1,1,MPI_COMM_WORLD,&request);

				//MPI_Isend(w->u0+ (w->sz+(w->dimz)/nproc*rank)*(2*w->sx+w->dimx)*(2*w->sy+w->dimy),4*(2*w->sx+w->dimx)*(2*w->sy+w->dimy),MPI_Real,rank-1,0,MPI_COMM_WORLD,&request);
				//MPI_Isend(w->u0+ (w->sz+(w->dimz)/nproc*rank)*(2*w->sx+w->dimx)*(2*w->sy+w->dimy),4*(2*w->sx+w->dimx)*(2*w->sy+w->dimy),MPI_Real,rank+1,0,MPI_COMM_WORLD,&request);

				MPI_Irecv(w->u0 + IND(w->dimz-4,w->dimy-4,w->dimx-4),(2*w->sx+w->dimx)*(2*w->sy+w->dimy),MPI_Real,rank+1,0,MPI_COMM_WORLD,&request);
				MPI_Irecv(w->u0 + IND(-4,-4,-4),(2*w->sx+w->dimx)*(2*w->sy+w->dimy),MPI_Real,rank-1,1,MPI_COMM_WORLD,&request);


				MPI_Wait(&request,&status);
			}
		else if (rank==0){

				MPI_Isend(w->u0 + IND(w->dimz-4,w->dimy-4,w->dimz-4),(2*w->sx+w->dimx)*(2*w->sy+w->dimy),MPI_Real,1,1,MPI_COMM_WORLD,&request);
				MPI_Irecv(w->u0 + IND(-4,-4,-4),(2*w->sx+w->dimx)*(2*w->sy+w->dimy),MPI_Real,1,0,MPI_COMM_WORLD,&request);

				MPI_Wait(&request,&status);
		}
		else {

			MPI_Isend(w->u0 + IND(0,0,0),(2*w->sx+w->dimx)*(2*w->sy+w->dimy),MPI_Real,rank-1,0,MPI_COMM_WORLD,&request);
			MPI_Irecv(w->u0 + IND(w->dimx-4,w->dimy-4,w->dimx-4),(2*w->sx+w->dimx)*(2*w->sy+w->dimy),MPI_Real,rank-1,1,MPI_COMM_WORLD,&request);

			MPI_Wait(&request,&status);
		}
	//        if((snapshot_enabled) && (t%nb_snap == 0))  pwave_snapshot(w);
		}
	}
	
	
	/** GATHER ***/
	
	 MPI_Gather(w->u0, wsize/nproc, MPI_Real,w->u0, wsize/nproc, MPI_Real,0, MPI_COMM_WORLD);
	 
	if (!rank)time_cpu = stop_timer(unit);
//fprintf(stdout, "... CPU time: %8.2f.\n", time_cpu);
//}




//printf("------------	END SIMULATION	------------------\n");

unsigned long int taille = w->size;

if (!rank){
    if(verbose) {
        fprintf(stdout, "... Total time: %8.2f %s.\n", time_cpu, timing_unit_string(unit));
        fprintf(stdout, "... check mode is: %s\n", check ? "ON":"OFF");
        if (check) {
        	//vérification
/*			char c[30];
			sprintf(c,"reference1_d%d_i%d", w->dimx,w->time_steps);
			FILE *f1 = fopen(c,"rb");
			printf("Récuperation du fichier référence!!\n");
			DATATYPE *u0 = malloc(sizeof(DATATYPE)*taille);
			int i = fread(u0,sizeof(DATATYPE),taille,f1);
//			if (i){ printf(" Fichier %s introuvable!!\n",c);goto out;}
			printf("Fin de la récuperation!!\n");
			//epsilon = 0.000000001;


			int nberreur = check_buffer(wg->u0,w->u0, taille, epsilon),i;
*/		  	fprintf(stdout, "... errors margin: %g\n", epsilon);
//			fprintf(stdout,"elements total: %d\nelements differents: %d\n", w->size,nberreur);
		}

    } else {
        //fprintf(stdout, "%-4u %-4u %-4u %-8.4g\n", w->dimx, w->dimy, w->dimz, time);
    }
}
//    pwave_delete(w, p);
    parser_delete(p);
    MPI_Finalize();
//    return EXIT_SUCCESS;
	return (int) (time*1000.);
}


