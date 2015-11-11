#ifndef _COMMON_H_
#define _COMMON_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <stdarg.h>
#include <mpi.h>
#include <omp.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#ifndef ALLOC_FACTOR
#define ALLOC_FACTOR 1.4
#endif //ALLOC_FACTOR

#ifdef _LONGIDS
typedef long long lint;
typedef unsigned long long ulint;
#else //_LONGIDS
typedef int lint;
typedef unsigned int ulint;
#endif //_LONGIDS

typedef float complex fcomplex;

extern int NNodes;
extern int NodeThis;

extern int Nx_here; //Number of planes in this node
extern int Ix0_here; //First plane in this node
extern int Ny_here;
extern int Iy0_here;
extern int Ngrid;
extern float Lbox;

extern float *Dens_local;
extern fftwf_complex *Cdens_local;

extern int TaskSmooth;
extern float *Dens_sm_local;
extern fftwf_complex *Cdens_sm_local;

extern int TaskVel;
extern float **Vel_local;

extern int TaskTidal;
extern int TaskTidalDiag;
extern float **Tid_local;
extern fftwf_complex **Ctid_local;

extern int TaskLinvel;
extern float **Lvel_local;
extern fftwf_complex **Clvel_local;

extern ulint Npart_alloc;
extern ulint Npart_saved;
extern ulint Npart_total;
extern float *Pos;
extern float *Vel;
extern ulint *Ids;

//Found in common.c
void report_error(int level,char *fmt,...);
FILE *my_fopen(const char *path,const char *mode);
void *my_malloc(size_t size);
void *my_calloc(size_t nmemb,size_t size);
size_t my_fread(void *ptr,size_t size,size_t nmemb,FILE *stream);
size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream);
void mpi_init(int* p_argc,char*** p_argv);
void domain_decomp(void);
void free_particles(void);
void end_all(void);

//Found in io.c
typedef struct {
  unsigned int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  unsigned int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHigh[6];
  int flag_entropy_instead_u;
  char fill[60];
  // fills to 256 Bytes
} gad_header;
void read_gadget_header(char *prefix,int input,gad_header *head_out);
void read_gadget(char *prefix,int input,ulint *npart,
		 gad_header *head_out,int interp_mode);
void write_output(char *prefix);

//Found in grid_tools.c
void compute_overdensity(ulint np,ulint np_total,float *pos,float *delta,int interp_order);
void compute_velocity_and_overdensity(ulint np,ulint np_total,float *pos,float *vel,
				      float *delta,float **velgrid,int interp_order);
void smooth_density_fourier(float r_smooth);
void get_smoothed_density_real(void);
void get_tidal_field(void);
void get_linearized_velocity(gad_header head);

#endif //_COMMON_H_
