#include "common.h"

typedef struct {
  char label[4];
  int size;
} gad_title;

static void gad_check_block(int b1,int b2)
{
  //////
  // Checks that a block from a snapshot is
  // consistent from its begin/end values
  if(b1!=b2)
    report_error(1,"Corrupted block!\n");
}

static int gad_seek_block(FILE *snap,char name[])
{
  //////
  // Seeks block from title
  gad_title tit;
  int block1,block2;

  rewind(snap);

  while(1>0) {
    if(!(fread(&block1,sizeof(int),1,snap))||
       feof(snap)||ferror(snap)) {
      report_error(1,"Block %s not found!!\n",name);
    }
    my_fread(&tit,sizeof(gad_title),1,snap);
    my_fread(&block2,sizeof(int),1,snap);
    gad_check_block(block1,block2);
    if(strncmp(tit.label,name,2)!=0)
      fseek(snap,tit.size,SEEK_CUR);
    else
      break;
  }

  return 0;
}

static int check_num_files(char *prefix)
{
  FILE *fil;

  fil=fopen(prefix,"r");
  if(fil!=NULL) {
    fclose(fil);
    return 1;
  }
  else {
    int nfils=0;
    while(nfils>=0) {
      char fname[256];
      sprintf(fname,"%s.%d",prefix,nfils);
      fil=fopen(fname,"r");
      if(fil!=NULL) {
	fclose(fil);
	nfils++;
      }
      else {
	if(nfils==0) {
	  report_error(0,"Can't find file %s or %s.x\n",prefix,prefix);
	  return -1;
	}
	else if(nfils==1) {
	  report_error(0,"Can only file %s found. Weird.\n",fname);
	  return -1;
	}
	else {
	  return nfils;
	}
      }
    }
  }
  
  report_error(0,"This shouldn't have happened \n");
  return -1;
}

void read_gadget_header(char *prefix,int input,gad_header *head_out)
{
  int ii,block1,block2;
  char fname[256];
  int nfils=check_num_files(prefix);
  unsigned long long nptot_64[6];
  gad_header head;
  FILE *snap;
  if(nfils<=0) 
    report_error(1,"Can't find snapshot files\n");
  else if(nfils==1)
    sprintf(fname,"%s",prefix);
  else
    sprintf(fname,"%s.0",prefix);

  snap=my_fopen(fname,"r");
  
  //Read header
  if(input==2)
    gad_seek_block(snap,"HEAD");
  my_fread(&block1,sizeof(int),1,snap);
  my_fread(&head,sizeof(gad_header),1,snap);
  my_fread(&block2,sizeof(int),1,snap);
  gad_check_block(block1,block2);
  fclose(snap);

  if(head.num_files!=nfils) {
    report_error(1,
		 "Header and existing files do not match %d != %d.\n"
		 "      There may be some files missing\n",
		 nfils,head.num_files);
  }

  if(NodeThis==0) {
    printf("  The cosmological model is:\n");
    printf("   - Omega_M = %.3lf\n",head.Omega0);
    printf("   - Omega_L = %.3lf\n",head.OmegaLambda);
    printf("   - h = %.3lf\n",head.HubbleParam);
    printf("  This file contains: \n");
    for(ii=0;ii<6;ii++) {
      nptot_64[ii]=head.npartTotal[ii]+((unsigned long long)(head.npartTotalHigh[ii]) << 32);
      printf("   - %llu particles of type %d with mass %.3lE\n",nptot_64[ii],
	     (int)ii,head.mass[ii]);
    }
    printf("  The box size is %.3lf\n",head.BoxSize);
    printf("  Redshift z = %.3lf \n",head.redshift);
    printf("\n");
  }

  *head_out=head;
  for(ii=0;ii<6;ii++) {
    (*head_out).npart[ii]=head.npartTotal[ii];
    (*head_out).num_files=1;
  }
}

void read_gadget(char *prefix,int input,ulint *npart,
		 gad_header *head_out,int interp_mode)
{
  int nfils=check_num_files(prefix);
  if(nfils<=0)
    report_error(1,"Can't find snapshot files\n");
  else {
    lint ii;
    char fname[256];
    gad_header head;
    int block1,block2;
    FILE *snap;
    float norm_vel;
    float agrid=Lbox/Ngrid;
    float x_left_proper,x_right_proper,x_left_bleft,x_right_bleft,x_left_bright,x_right_bright;
    float x_min=0,x_max=Lbox,x_shift_bright=0,x_shift_bleft=0;

    x_left_proper=Ix0_here*agrid;
    x_right_proper=(Ix0_here+Nx_here)*agrid;
    x_left_bleft=2*Lbox;
    x_left_bright=2*Lbox;
    x_right_bleft=-2*Lbox;
    x_right_bright=-2*Lbox;
    if(interp_mode==0) {//NGP
      x_min=-agrid/2;
      x_max=Lbox-agrid/2;
      x_left_proper=(Ix0_here-0.5f)*agrid;
      x_right_proper=(Ix0_here+Nx_here-0.5f)*agrid;
    }
    else if(interp_mode==1) {//CIC
      x_min=0;
      x_max=Lbox;
      x_left_proper=Ix0_here*agrid;
      x_right_proper=(Ix0_here+Nx_here)*agrid;
      x_left_bleft=x_left_proper-agrid;
      x_right_bleft=x_left_proper;
      if(NNodes>1 && NodeThis==0) {
	x_left_bleft=(Ngrid-1)*agrid;
	x_right_bleft=Lbox;
	x_shift_bleft=-Lbox;
      }
    }
    else if(interp_mode==2) {//TSC
      x_min=-agrid/2;
      x_max=Lbox-agrid/2;
      x_left_proper=(Ix0_here-0.5f)*agrid;
      x_right_proper=(Ix0_here+Nx_here-0.5f)*agrid;
      x_left_bleft=x_left_proper-agrid;
      x_left_bright=x_left_proper;
      x_right_bleft=x_right_proper;
      x_right_bright=x_right_proper+agrid;
      if(NNodes>1 && NodeThis==0) {
	x_left_bleft=(Ngrid-1.5f)*agrid;
	x_right_bleft=(Ngrid-0.5f)*agrid;
	x_shift_bleft=-Lbox;
      }
      if(NNodes>1 && NodeThis==NNodes-1) {
	x_left_bright=-0.5f*agrid;
	x_right_bright=0.5f*agrid;
	x_shift_bright=Lbox;
      }
    }
    else
      report_error(1,"Wrong interpolation scheme %d\n",interp_mode);

    if(NodeThis==0)
      printf("  Reading %d snapshot files \n",nfils);

    if(nfils>1)
      sprintf(fname,"%s.0",prefix);
    else
      sprintf(fname,"%s",prefix);
    snap=my_fopen(fname,"r");
  
    //Read header
    if(input==2)
      gad_seek_block(snap,"HEAD");
    my_fread(&block1,sizeof(int),1,snap);
    my_fread(&head,sizeof(gad_header),1,snap);
    my_fread(&block2,sizeof(int),1,snap);
    gad_check_block(block1,block2);

    if(head.num_files!=nfils) {
      report_error(1,
		   "Header and existing files do not match %d != %d.\n"
		   "      There may be some files missing\n",
		   nfils,head.num_files);
    }

    *head_out=head;
    for(ii=0;ii<6;ii++) {
      (*head_out).npart[ii]=head.npartTotal[ii];
      (*head_out).num_files=1;
    }

    norm_vel=(float)(1./sqrt(head.time));

    fclose(snap);

    ulint np_read=0;
    ulint np_saved=0;
    for(ii=0;ii<nfils;ii++) {
      gad_header head_dum;
      ulint np_new,np_got_here=0,np_got_here_b;
      lint jj;

      if(nfils>1)
	sprintf(fname,"%s.%d",prefix,(int)ii);
      else
	sprintf(fname,"%s",prefix);
      snap=my_fopen(fname,"r");
#ifdef _DEBUG
      printf("  Node %d, reading file  %s \n",NodeThis,fname);
#endif //_DEBUG

      //Read header
      if(input==2)
	gad_seek_block(snap,"HEAD");
      my_fread(&block1,sizeof(int),1,snap);
      my_fread(&head_dum,sizeof(gad_header),1,snap);
      my_fread(&block2,sizeof(int),1,snap);
      gad_check_block(block1,block2);
      
      np_new=0;
      for(jj=0;jj<6;jj++)
	np_new+=head_dum.npart[jj];
#ifdef _DEBUG
      printf("  Node %d : %ld parts in file %ld \n",NodeThis,(long)np_new,(long)ii);
#endif //_DEBUG

      if(np_read+np_new>Npart_total) {
	report_error(1,
		     "files seem to contain too many particles\n"
		     "      file %s, %ld > %ld \n",
		     fname,(long)(np_read+np_new),(long)(Npart_total));
      }

      if(input==2)
	gad_seek_block(snap,"POS");
      my_fread(&block1,sizeof(int),1,snap);
      for(jj=0;jj<np_new;jj++) {
	int ax;
	int include_part=0;
	float *p=&(Pos[3*np_saved]);
	my_fread(p,sizeof(float),3,snap);
	for(ax=0;ax<3;ax++) {
	  if(p[ax]>=x_max)
	    p[ax]-=Lbox;
	  else if(p[ax]<x_min)
	    p[ax]+=Lbox;
	}
	if(p[0]>=x_left_proper && p[0]<x_right_proper)
	  include_part=1;
	else if(p[0]>=x_left_bleft && p[0]<x_right_bleft) {
	  include_part=1;
	  p[0]+=x_shift_bleft;
	}
	else if(p[0]>=x_left_bright && p[0]<x_right_bright) {
	  include_part=1;
	  p[0]+=x_shift_bright;
	}
	if(include_part) {
	  Ids[np_saved]=jj;
	  np_saved++;
	  np_got_here++;
	  if(np_saved>=Npart_alloc)
	    report_error(1,"Too many particles, increase buffer size %llu %llu\n",np_saved,Npart_alloc);
	}
      }
      my_fread(&block2,sizeof(int),1,snap);
      gad_check_block(block1,block2);

      np_read+=np_new;
      if(np_got_here==0) {
	fclose(snap);
	continue;
      }	

      np_got_here_b=0;
      if(input==2)
	gad_seek_block(snap,"VEL");
      my_fread(&block1,sizeof(int),1,snap);
      for(jj=0;jj<np_new;jj++) {
	float v[3];
	lint new_index=np_saved-np_got_here+np_got_here_b;
	my_fread(v,sizeof(float),3,snap);
	if(Ids[new_index]==jj) {
	  int ax;
	  np_got_here_b++;
	  for(ax=0;ax<3;ax++)
	    Vel[3*new_index+ax]=norm_vel*v[ax];
	}
      }
      if(np_got_here_b!=np_got_here)
	report_error(1,"This shouldn't have happened\n");
      my_fread(&block2,sizeof(int),1,snap);
      gad_check_block(block1,block2);

      np_got_here_b=0;
      if(input==2)
	gad_seek_block(snap,"ID");
      my_fread(&block1,sizeof(int),1,snap);
      for(jj=0;jj<np_new;jj++) {
	ulint id;
	lint new_index=np_saved-np_got_here+np_got_here_b;
	my_fread(&id,sizeof(ulint),1,snap);
	if(Ids[new_index]==jj) {
	  np_got_here_b++;
	  Ids[new_index]=id;
	}
      }
      if(np_got_here_b!=np_got_here)
	report_error(1,"This shouldn't have happened\n");
      my_fread(&block2,sizeof(int),1,snap);
      gad_check_block(block1,block2);

      fclose(snap);
    }
    if(np_read!=Npart_total) {
      report_error(1,
		   "#particles read disagrees with header: %ld != %ld\n",
		   (long)np_read,(long)(Npart_total));
    }
    *npart=np_saved;
  }
}

typedef struct {
  float eval;
  float evec[3];
} Esys;

static int compare_evals(const void *p1,const void *p2)
{
  Esys *e1=(Esys *)p1;
  Esys *e2=(Esys *)p2;

  return (e1->eval>e2->eval) - (e1->eval<e2->eval);
}

void write_output(char *prefix)
{
  FILE *fo;
  char fname[256];
  int ix,num_grids;

  sprintf(fname,"%s_dens.%04d",prefix,NodeThis);
  num_grids=1;
  fo=my_fopen(fname,"wb");
  my_fwrite(&(num_grids),sizeof(int),1,fo);
  my_fwrite(&(Ngrid),sizeof(int),1,fo);
  my_fwrite(&(Nx_here),sizeof(int),1,fo);
  for(ix=0;ix<Nx_here;ix++) {
    int iy;
    int ix_id=Ix0_here+ix;
    my_fwrite(&(ix_id),sizeof(int),1,fo);
    for(iy=0;iy<Ngrid;iy++) {
      lint index0=2*(Ngrid/2+1)*((lint)(iy+ix*Ngrid));
      my_fwrite(&(Dens_local[index0]),sizeof(float),Ngrid,fo);
    }
  }
  fclose(fo);

  if(TaskVel) {
    sprintf(fname,"%s_vel.%04d",prefix,NodeThis);
    num_grids=3;
    fo=my_fopen(fname,"wb");
    my_fwrite(&(num_grids),sizeof(int),1,fo);
    my_fwrite(&(Ngrid),sizeof(int),1,fo);
    my_fwrite(&(Nx_here),sizeof(int),1,fo);
    for(ix=0;ix<Nx_here;ix++) {
      int iy;
      int ix_id=Ix0_here+ix;
      my_fwrite(&(ix_id),sizeof(int),1,fo);
      for(iy=0;iy<Ngrid;iy++) {
	int iz;
	for(iz=0;iz<Ngrid;iz++) {
	  lint index=iz+2*(Ngrid/2+1)*((lint)(iy+ix*Ngrid));
	  my_fwrite(&(Vel_local[0][index]),sizeof(float),1,fo);
	  my_fwrite(&(Vel_local[1][index]),sizeof(float),1,fo);
	  my_fwrite(&(Vel_local[2][index]),sizeof(float),1,fo);
	}
      }
    }
    fclose(fo);
  }

  if(TaskSmooth) {
    sprintf(fname,"%s_dens_sm.%04d",prefix,NodeThis);
    num_grids=1;
    fo=my_fopen(fname,"wb");
    my_fwrite(&(num_grids),sizeof(int),1,fo);
    my_fwrite(&(Ngrid),sizeof(int),1,fo);
    my_fwrite(&(Nx_here),sizeof(int),1,fo);
    for(ix=0;ix<Nx_here;ix++) {
      int iy;
      int ix_id=Ix0_here+ix;
      my_fwrite(&(ix_id),sizeof(int),1,fo);
      for(iy=0;iy<Ngrid;iy++) {
	lint index0=2*(Ngrid/2+1)*((lint)(iy+ix*Ngrid));
	my_fwrite(&(Dens_sm_local[index0]),sizeof(float),Ngrid,fo);
      }
    }
    fclose(fo);
  }

  if(TaskFractalD) {
    sprintf(fname,"%s_dens_smD.%04d",prefix,NodeThis);
    num_grids=1;
    fo=my_fopen(fname,"wb");
    my_fwrite(&(num_grids),sizeof(int),1,fo);
    my_fwrite(&(Ngrid),sizeof(int),1,fo);
    my_fwrite(&(Nx_here),sizeof(int),1,fo);
    for(ix=0;ix<Nx_here;ix++) {
      int iy;
      int ix_id=Ix0_here+ix;
      my_fwrite(&(ix_id),sizeof(int),1,fo);
      for(iy=0;iy<Ngrid;iy++) {
	lint index0=2*(Ngrid/2+1)*((lint)(iy+ix*Ngrid));
	my_fwrite(&(Dens_smD_local[index0]),sizeof(float),Ngrid,fo);
      }
    }
    fclose(fo);
  }

  if(TaskLinvel) {
    sprintf(fname,"%s_linvel.%04d",prefix,NodeThis);
    num_grids=3;
    fo=my_fopen(fname,"wb");
    my_fwrite(&(num_grids),sizeof(int),1,fo);
    my_fwrite(&(Ngrid),sizeof(int),1,fo);
    my_fwrite(&(Nx_here),sizeof(int),1,fo);
    for(ix=0;ix<Nx_here;ix++) {
      int iy;
      int ix_id=Ix0_here+ix;
      my_fwrite(&(ix_id),sizeof(int),1,fo);
      for(iy=0;iy<Ngrid;iy++) {
	int iz;
	for(iz=0;iz<Ngrid;iz++) {
	  lint index=iz+2*(Ngrid/2+1)*((lint)(iy+ix*Ngrid));
	  my_fwrite(&(Lvel_local[0][index]),sizeof(float),1,fo);
	  my_fwrite(&(Lvel_local[1][index]),sizeof(float),1,fo);
	  my_fwrite(&(Lvel_local[2][index]),sizeof(float),1,fo);
	}
      }
    }
    fclose(fo);
  }

  if(TaskNlvel) {
    sprintf(fname,"%s_nlvel.%04d",prefix,NodeThis);
    num_grids=3;
    fo=my_fopen(fname,"wb");
    my_fwrite(&(num_grids),sizeof(int),1,fo);
    my_fwrite(&(Ngrid),sizeof(int),1,fo);
    my_fwrite(&(Nx_here),sizeof(int),1,fo);
    for(ix=0;ix<Nx_here;ix++) {
      int iy;
      int ix_id=Ix0_here+ix;
      my_fwrite(&(ix_id),sizeof(int),1,fo);
      for(iy=0;iy<Ngrid;iy++) {
	int iz;
	for(iz=0;iz<Ngrid;iz++) {
	  lint index=iz+2*(Ngrid/2+1)*((lint)(iy+ix*Ngrid));
	  my_fwrite(&(Nlvel_local[0][index]),sizeof(float),1,fo);
	  my_fwrite(&(Nlvel_local[1][index]),sizeof(float),1,fo);
	  my_fwrite(&(Nlvel_local[2][index]),sizeof(float),1,fo);
	}
      }
    }
    fclose(fo);
  }

  if(TaskTidal && !TaskTidalDiag) {
    sprintf(fname,"%s_tidal.%04d",prefix,NodeThis);
    num_grids=6;
    fo=my_fopen(fname,"wb");
    my_fwrite(&(num_grids),sizeof(int),1,fo);
    my_fwrite(&(Ngrid),sizeof(int),1,fo);
    my_fwrite(&(Nx_here),sizeof(int),1,fo);
    for(ix=0;ix<Nx_here;ix++) {
      int iy;
      int ix_id=Ix0_here+ix;
      my_fwrite(&(ix_id),sizeof(int),1,fo);
      for(iy=0;iy<Ngrid;iy++) {
	int iz;
	for(iz=0;iz<Ngrid;iz++) {
	  lint index=iz+2*(Ngrid/2+1)*((lint)(iy+ix*Ngrid));
	  my_fwrite(&(Tid_local[0][index]),sizeof(float),1,fo);
	  my_fwrite(&(Tid_local[1][index]),sizeof(float),1,fo);
	  my_fwrite(&(Tid_local[2][index]),sizeof(float),1,fo);
	  my_fwrite(&(Tid_local[3][index]),sizeof(float),1,fo);
	  my_fwrite(&(Tid_local[4][index]),sizeof(float),1,fo);
	  my_fwrite(&(Tid_local[5][index]),sizeof(float),1,fo);
	}
      }
    }
    fclose(fo);
  }

  if(TaskTidalDiag) {
    Esys leig[3];
    gsl_vector *eval=gsl_vector_alloc(3);
    gsl_matrix *tij=gsl_matrix_alloc(3,3);
    gsl_matrix *evec=gsl_matrix_alloc(3,3);
    gsl_eigen_symmv_workspace *w=gsl_eigen_symmv_alloc(3);

    sprintf(fname,"%s_tidal_eigsys.%04d",prefix,NodeThis);
    num_grids=12;
    fo=my_fopen(fname,"wb");
    my_fwrite(&(num_grids),sizeof(int),1,fo);
    my_fwrite(&(Ngrid),sizeof(int),1,fo);
    my_fwrite(&(Nx_here),sizeof(int),1,fo);

    for(ix=0;ix<Nx_here;ix++) {
      int iy;
      int ix_id=Ix0_here+ix;
      my_fwrite(&(ix_id),sizeof(int),1,fo);
      for(iy=0;iy<Ngrid;iy++) {
	int iz;
	for(iz=0;iz<Ngrid;iz++) {
	  int im;
	  lint index=iz+2*(Ngrid/2+1)*((lint)(iy+ix*Ngrid));
	  gsl_matrix_set(tij,0,0,Tid_local[0][index]);
	  gsl_matrix_set(tij,1,1,Tid_local[1][index]);
	  gsl_matrix_set(tij,2,2,Tid_local[2][index]);
	  gsl_matrix_set(tij,0,1,Tid_local[3][index]);
	  gsl_matrix_set(tij,1,0,Tid_local[3][index]);
	  gsl_matrix_set(tij,1,2,Tid_local[4][index]);
	  gsl_matrix_set(tij,2,1,Tid_local[4][index]);
	  gsl_matrix_set(tij,0,2,Tid_local[5][index]);
	  gsl_matrix_set(tij,2,0,Tid_local[5][index]);
	  gsl_eigen_symmv(tij,eval,evec,w);

	  for(im=0;im<3;im++) {
	    leig[im].eval=(float)(gsl_vector_get(eval,im));
	    leig[im].evec[0]=(float)(gsl_matrix_get(evec,0,im));
	    leig[im].evec[1]=(float)(gsl_matrix_get(evec,1,im));
	    leig[im].evec[2]=(float)(gsl_matrix_get(evec,2,im));
	  }
	  qsort(leig,3,sizeof(Esys),compare_evals);

	  for(im=0;im<3;im++) {
	    my_fwrite(&(leig[im].eval),sizeof(float),1,fo);
	    my_fwrite(&(leig[im].evec[0]),sizeof(float),1,fo);
	    my_fwrite(&(leig[im].evec[1]),sizeof(float),1,fo);
	    my_fwrite(&(leig[im].evec[2]),sizeof(float),1,fo);
	  }
	}
      }
    }
    fclose(fo);
    gsl_eigen_symmv_free(w);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(tij);
  }
}
