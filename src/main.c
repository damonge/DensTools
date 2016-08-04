#include "common.h"

int main(int argc,char **argv)
{
  char fnameIn[256]="i.snap";
  char fnameOut[256]="output";
  char interp_method[256]="NGP";
  float r_smooth=0.;
  gad_header head;
  int input_format=1;
  int ii,interp_order=0;

  char **c;
  mpi_init(&argc,&argv);

  for(c=argv+1;*c;c++) {
    if(!strcmp(*c,"-in")) sprintf(fnameIn,"%s",*++c);
    else if(!strcmp(*c,"-out")) sprintf(fnameOut,"%s",*++c);
    else if(!strcmp(*c,"-ngrid")) Ngrid=atoi(*++c);
    else if(!strcmp(*c,"-smooth")) r_smooth=(float)(atof(*++c));
    else if(!strcmp(*c,"-interp")) sprintf(interp_method,"%s",*++c);
    else if(!strcmp(*c,"-diag_tidal")) TaskTidalDiag=1;
    else if(!strcmp(*c,"-use_finite_differences")) UseFD=1;
    else if(!strcmp(*c,"-do")) {
      int ich=0;
      char d=(*++c)[ich];
      while(d!='\0') {
	if(d=='t')
	  TaskTidal=1;
	else if(d=='v')
	  TaskVel=1;
	else if(d=='l')
	  TaskLinvel=1;
	else if(d=='u')
	  TaskNlvel=1;
	else
	  fprintf(stderr,"Unknown task %c. Supported : v, t, l\n",d);
	d=(*c)[++ich];
      }
    }
    else if(!strcmp(*c,"-in_fmt")) input_format=atoi(*++c);
    else if(!strcmp(*c,"-h")) {
      if(NodeThis==0) {
	fprintf(stderr,"Usage: DensTools -<opt-name> <option>\n");
	fprintf(stderr,"Options:\n");
	fprintf(stderr,"  -in         -> input file/prefix\n");
	fprintf(stderr,"  -in_fmt     -> input GADGET format\n");
	fprintf(stderr,"  -out        -> output prefix\n");
	fprintf(stderr,"  -ngrid      -> #cells per side\n");
	fprintf(stderr,"  -smooth     -> smoothing length (in same units as simulation box)\n");
	fprintf(stderr,"  -interp     -> particle interpolation scheme: NGP, CIC or TSC\n");
	fprintf(stderr,"  -do         -> concatenate tasks: velocity (v), tidal field (t),"
		" linearized velocity (l), non-linear velocity (u)\n");
	fprintf(stderr,"  -diag_tidal -> diagonalize tidal tensor\n");
	fprintf(stderr,"  -use_finite_differences -> use central FDs for derivatives\n");
	fprintf(stderr,"  -h          -> this help\n\n");
	return 0;
      }
    }
    else {
      fprintf(stderr,"Unknown option %s\n",*c);
      exit(1);
    }
  }
  if(!TaskTidal)
    TaskTidalDiag=0;
  if(r_smooth>0 || TaskTidal+TaskLinvel+TaskNlvel)
    TaskSmooth=1;
  if(strcmp(interp_method,"NGP") && strcmp(interp_method,"CIC") && strcmp(interp_method,"TSC"))
    report_error(1,"Unknown interpolation scheme %s. Use NGP, CIC or TSC\n",interp_method);

  if(NodeThis==0) {
    printf("* Run options :\n");
    printf("   Input prefix : %s\n",fnameIn);
    printf("   Output prefix : %s\n",fnameOut);
    printf("   Ngrid : %d\n",Ngrid);
    printf("   R_smooth : %.2lf \n",r_smooth);
    printf("   Interpolation scheme : %s\n",interp_method);
    printf("   Tasks to do:");
    if(TaskVel)
      printf(" V");
    if(TaskSmooth)
      printf(" S");
    if(TaskTidal && !TaskTidalDiag)
      printf(" T");
    if(TaskTidalDiag)
      printf(" TD");
    if(TaskLinvel)
      printf(" LV");
    if(TaskNlvel)
      printf(" NLV");
    printf("\n");
    printf("   Use finite differences : %d\n",UseFD);
    printf("\n");
  }
  if(!strcmp(interp_method,"NGP"))
    interp_order=0;
  else if(!strcmp(interp_method,"CIC"))
    interp_order=1;
  else if(!strcmp(interp_method,"TSC"))
    interp_order=2;
  else
    report_error(1,"Wrong interpolation order\n");

  if(NodeThis==0)
    printf("* Reading header\n");
  read_gadget_header(fnameIn,input_format,&head);
  Lbox=(float)(head.BoxSize);
  Npart_total=0;
  for(ii=0;ii<6;ii++)
    Npart_total+=(head.npartTotal[ii]+((unsigned long long)(head.npartTotalHigh[ii]) << 32));
  Npart_alloc=(unsigned long long)(Npart_total*ALLOC_FACTOR/((float)(NNodes)));

  if(NodeThis==0)
    printf("* Domain decomposition:\n");
  domain_decomp();

  if(NodeThis==0)
    printf("* Reading snapshot\n");
  read_gadget(fnameIn,input_format,&Npart_saved,&head,interp_order);
  if(NodeThis==0)
    printf("\n");

  if(TaskVel) {
    if(NodeThis==0)
      printf("* Computing velocity and density grids\n");
    compute_velocity_and_overdensity(Npart_saved,Npart_total,Pos,Vel,
				     Dens_local,Vel_local,interp_order);
  }
  else {
    if(NodeThis==0)
      printf("* Computing density grid\n");
    compute_overdensity(Npart_saved,Npart_total,Pos,Dens_local,interp_order);
  }
  free_particles();
  if(NodeThis==0)
    printf("\n");

  if(TaskSmooth) {
    if(NodeThis==0)
      printf("* Smoothing density field\n");
    smooth_density_fourier(r_smooth);
    if(NodeThis==0)
      printf("\n");
  }
  if(TaskTidal) {
    if(NodeThis==0)
      printf("* Getting tidal field\n");
    get_tidal_field();
    if(NodeThis==0)
      printf("\n");
  }
  if(TaskLinvel) {
    if(NodeThis==0)
      printf("* Getting linearized velocity\n");
    get_linearized_velocity(head);
    if(NodeThis==0)
      printf("\n");
  }
  if(TaskNlvel) {
    if(NodeThis==0)
      printf("* Getting non-linear velocity\n");
    get_nonlinear_velocity(head);
    if(NodeThis==0)
      printf("\n");
  }
  if(TaskSmooth && !TaskNlvel) {
    if(NodeThis==0)
      printf("* Transforming smoothed density field back\n");
    get_smoothed_density_real();
    if(NodeThis==0)
      printf("\n");
  }
  if(NodeThis==0)
    printf("* Writing output\n");
  write_output(fnameOut);
  
  if(NodeThis==0)
      printf("\nDone!\n");
  end_all();
  MPI_Finalize();

  return 0;
}
