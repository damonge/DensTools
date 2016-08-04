#include "common.h"

int NNodes;
int NodeThis;
int NodeLeft;
int NodeRight;

int UseFD=0;

int Ngrid=128;
int Nx_here,Ix0_here,Ny_here,Iy0_here;
float Lbox=128.;

float *Dens_local;
fftwf_complex *Cdens_local;

int TaskVel=0;
float **Vel_local;

int TaskSmooth=0;
float *Dens_sm_local;
fftwf_complex *Cdens_sm_local;

int TaskTidal=0;
int TaskTidalDiag=0;
float **Tid_local;
fftwf_complex **Ctid_local;

int TaskLinvel=0;
float **Lvel_local;
fftwf_complex **Clvel_local;

int TaskNlvel=0;
float *Nlpot_local;
fftwf_complex *Cnlpot_local;
float **Nlvel_local;
float *SliceLeft_Dens;
float *SliceRight_Dens;
float *SliceLeft_Nlpot;
float *SliceRight_Nlpot;

ulint Npart_alloc,Npart_saved,Npart_total;
float *Pos,*Vel;
ulint *Ids;

void free_particles(void)
{
  free(Pos);
  free(Vel);
  free(Ids);
}

void end_all(void)
{
  int i;

  fftwf_free(Dens_local);
  if(TaskSmooth)
    fftwf_free(Dens_sm_local);
  if(TaskTidal) {
    for(i=0;i<6;i++)
      fftwf_free(Tid_local[i]);
    free(Tid_local);
  }
  if(TaskLinvel) {
    for(i=0;i<3;i++)
      fftwf_free(Lvel_local[i]);
    free(Lvel_local);
  }
  if(TaskNlvel) {
    fftwf_free(Nlpot_local);
    for(i=0;i<3;i++)
      fftwf_free(Nlvel_local[i]);
    free(Nlvel_local);
    free(SliceLeft_Dens);
    free(SliceRight_Dens);
    free(SliceLeft_Nlpot);
    free(SliceRight_Nlpot);
  }
  if(TaskVel) {
    for(i=0;i<3;i++)
      fftwf_free(Vel_local[i]);
    free(Vel_local);
  }
}

void domain_decomp(void)
{
  int i;
  ptrdiff_t dsize,nx,ix0,ny,iy0;

  dsize=fftwf_mpi_local_size_3d_transposed(Ngrid,Ngrid,Ngrid,MPI_COMM_WORLD,
					   &nx,&ix0,&ny,&iy0);
  Nx_here=(int)nx;
  Ix0_here=(int)ix0;
  Ny_here=(int)ny;
  Iy0_here=(int)iy0;

  MPI_Barrier(MPI_COMM_WORLD);
  printf("   Node %d; ix0=%d, nx=%d; iy0=%d, ny=%d;\n",
	 NodeThis,Ix0_here,Nx_here,Iy0_here,Ny_here);
  MPI_Barrier(MPI_COMM_WORLD);

  if(NodeThis==0)
    printf("\n");
  Dens_local=fftwf_alloc_real(2*dsize);
  if(Dens_local==NULL)
    report_error(1,"Couldn't allocate density field\n");
  Cdens_local=(fftwf_complex *)Dens_local;
  if(TaskSmooth) {
    Dens_sm_local=fftwf_alloc_real(2*dsize);
    if(Dens_sm_local==NULL)
      report_error(1,"Couldn't allocate density field\n");
    Cdens_sm_local=(fftwf_complex *)Dens_sm_local;
  }
  if(TaskTidal) {
    Tid_local=my_malloc(6*sizeof(float *));
    Ctid_local=my_malloc(6*sizeof(fftwf_complex *));
    for(i=0;i<6;i++) {
      Tid_local[i]=fftwf_alloc_real(2*dsize);
      if(Tid_local[i]==NULL)
	report_error(1,"Couldn't allocate density field\n");
      Ctid_local[i]=(fftwf_complex *)(Tid_local[i]);
    }
  }
  if(TaskLinvel) {
    Lvel_local=my_malloc(3*sizeof(float *));
    Clvel_local=my_malloc(3*sizeof(fftwf_complex *));
    for(i=0;i<3;i++) {
      Lvel_local[i]=fftwf_alloc_real(2*dsize);
      if(Lvel_local[i]==NULL)
	report_error(1,"Couldn't allocate density field\n");
      Clvel_local[i]=(fftwf_complex *)(Lvel_local[i]);
    }
  }
  if(TaskNlvel) {
    long slice_size=Ngrid*2*(Ngrid/2+1);
    Nlpot_local=fftwf_alloc_real(2*dsize);
    if(Nlpot_local==NULL)
      report_error(1,"Couldn't allocate density field\n");
    Cnlpot_local=(fftwf_complex *)Nlpot_local;
    Nlvel_local=my_malloc(3*sizeof(float *));
    for(i=0;i<3;i++) {
      Nlvel_local[i]=fftwf_alloc_real(2*dsize);
      if(Nlvel_local[i]==NULL)
	report_error(1,"Couldn't allocate density field\n");
    }
    SliceLeft_Dens=my_malloc(slice_size*sizeof(float));
    SliceRight_Dens=my_malloc(slice_size*sizeof(float));
    SliceLeft_Nlpot=my_malloc(slice_size*sizeof(float));
    SliceRight_Nlpot=my_malloc(slice_size*sizeof(float));
  }
  if(TaskVel) {
    Vel_local=my_malloc(3*sizeof(float *));
    for(i=0;i<3;i++) {
      Vel_local[i]=fftwf_alloc_real(2*dsize);
      if(Vel_local[i]==NULL)
	report_error(1,"Couldn't allocate density field\n");
    }
  }

  Pos=my_malloc(3*Npart_alloc*sizeof(float));
  Vel=my_malloc(3*Npart_alloc*sizeof(float));
  Ids=my_malloc(Npart_alloc*sizeof(ulint));
}

void report_error(int level,char *fmt,...)
{
  va_list args;
  char msg[256];

  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);
  
  if(level) {
    fprintf(stderr," Node %d, fatal: %s",NodeThis,msg);
    exit(level);
  }
  else
    fprintf(stderr," Node %d, Warning: %s",NodeThis,msg);
}

void *my_malloc(size_t size)
{
  void *outptr=malloc(size);
  if(outptr==NULL) report_error(1,"Out of memory\n");

  return outptr;
}

size_t my_fread(void *ptr,size_t size,size_t nmemb,FILE *stream)
{
  if(fread(ptr,size,nmemb,stream)!=nmemb)
    report_error(1,"Error freading\n");

  return nmemb;
}

size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  if(fwrite(ptr,size,nmemb,stream)!=nmemb)
    report_error(1,"Error fwriting\n");

  return nmemb;
}

void *my_calloc(size_t nmemb,size_t size)
{
  void *outptr=calloc(nmemb,size);
  if(outptr==NULL)
    report_error(1,"Out of memory\n");

  return outptr;
}

FILE *my_fopen(const char *path,const char *mode)
{
  FILE *fout=fopen(path,mode);
  if(fout==NULL)
    report_error(1,"Couldn't open file %s\n",path);

  return fout;
}

void mpi_init(int* p_argc,char*** p_argv)
{
  MPI_Init(p_argc,p_argv);
  fftwf_mpi_init();

  MPI_Comm_size(MPI_COMM_WORLD,&NNodes);
  MPI_Comm_rank(MPI_COMM_WORLD,&NodeThis);
  if(NodeThis==0)
    NodeLeft=NNodes-1;
  else
    NodeLeft=NodeThis-1;
  if(NodeThis==NNodes-1)
    NodeRight=0;
  else
    NodeRight=NodeThis+1;
}
