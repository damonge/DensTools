#include "common.h"

static void pos_2_ngp_serial(ulint np,float *pos,float *delta)
{
  if(NodeThis==0)
    printf("  Calculating NGP\n");
  memset(delta,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  
  lint ii;
  float i_agrid=1./Ngrid;

  for(ii=0;ii<np;ii++) {
    lint index;
    int ax,i0[3];
    float *p=&(pos[3*ii]);

    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(p[ax]*i_agrid+0.5);
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      else if(i0[ax]<0) i0[ax]+=Ngrid;
    }

    index=i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]));
    delta[index]++;
  }
}

static void pos_2_ngp_parallel(ulint np,float *pos,float *delta)
{
  lint ii;
  float i_agrid=1./Ngrid;

  if(NodeThis==0)
    printf("  Calculating NGP\n");
  memset(delta,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  
  for(ii=0;ii<np;ii++) {
    lint index;
    int ax,i0[3];
    float *p=&(pos[3*ii]);

    for(ax=0;ax<3;ax++)
      i0[ax]=(int)(p[ax]*i_agrid+0.5);
    if(i0[1]>=Ngrid) i0[1]-=Ngrid;
    if(i0[1]<0) i0[1]+=Ngrid;
    if(i0[2]>=Ngrid) i0[2]-=Ngrid;
    if(i0[2]<0) i0[2]+=Ngrid;

    if(i0[0]>=0 && i0[0]<Nx_here) {
      index=i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]));
      delta[index]++;
    }
  }
}

static void pos_2_cic_serial(ulint np,float *pos,float *delta)
{
  if(NodeThis==0)
    printf("  Calculating CIC\n");
  memset(delta,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));

  lint ii;
  float i_agrid=Ngrid/Lbox;

  for(ii=0;ii<np;ii++) {
    int ax,i0[3],i1[3];
    float a0[3],a1[3];
    float *p=&(pos[3*ii]);

    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(p[ax]*i_agrid);
      a1[ax]=p[ax]*i_agrid-i0[ax];
      a0[ax]=1-a1[ax];
      i1[ax]=i0[ax]+1;
      if(i0[ax]<0) i0[ax]+=Ngrid;
      if(i1[ax]<0) i1[ax]+=Ngrid;
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      if(i1[ax]>=Ngrid) i1[ax]-=Ngrid;
    }
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]))]+=a0[2]*a0[1]*a0[0];
    delta[i1[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]))]+=a1[2]*a0[1]*a0[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(i1[1]+Ngrid*i0[0]))]+=a0[2]*a1[1]*a0[0];
    delta[i1[2]+2*(Ngrid/2+1)*((lint)(i1[1]+Ngrid*i0[0]))]+=a1[2]*a1[1]*a0[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i1[0]))]+=a0[2]*a0[1]*a1[0];
    delta[i1[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i1[0]))]+=a1[2]*a0[1]*a1[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(i1[1]+Ngrid*i1[0]))]+=a0[2]*a1[1]*a1[0];
    delta[i1[2]+2*(Ngrid/2+1)*((lint)(i1[1]+Ngrid*i1[0]))]+=a1[2]*a1[1]*a1[0];
  }
}

static void pos_2_cic_parallel(ulint np,float *pos,float *delta)
{
  if(NodeThis==0)
    printf("  Calculating CIC\n");
  memset(delta,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));

  lint ii;
  float i_agrid=Ngrid/Lbox;
  for(ii=0;ii<np;ii++) {
    int ax,i0[3],i1[3];
    float a0[3],a1[3];
    float *p=&(pos[3*ii]);

    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(floorf(p[ax]*i_agrid));
      a1[ax]=p[ax]*i_agrid-i0[ax];
      a0[ax]=1-a1[ax];
      i1[ax]=i0[ax]+1;
    }
    for(ax=1;ax<3;ax++) {
      if(i0[ax]<0) i0[ax]+=Ngrid;
      if(i1[ax]<0) i1[ax]+=Ngrid;
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      if(i1[ax]>=Ngrid) i1[ax]-=Ngrid;
    }
    i0[0]-=Ix0_here;
    i1[0]-=Ix0_here;

    if(i0[0]>=0 && i0[0]<Nx_here) {
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]))]+=a0[2]*a0[1]*a0[0];
      delta[i1[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]))]+=a1[2]*a0[1]*a0[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(i1[1]+Ngrid*i0[0]))]+=a0[2]*a1[1]*a0[0];
      delta[i1[2]+2*(Ngrid/2+1)*((lint)(i1[1]+Ngrid*i0[0]))]+=a1[2]*a1[1]*a0[0];
    }
    if(i1[0]>=0 && i1[0]<Nx_here) {
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i1[0]))]+=a0[2]*a0[1]*a1[0];
      delta[i1[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i1[0]))]+=a1[2]*a0[1]*a1[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(i1[1]+Ngrid*i1[0]))]+=a0[2]*a1[1]*a1[0];
      delta[i1[2]+2*(Ngrid/2+1)*((lint)(i1[1]+Ngrid*i1[0]))]+=a1[2]*a1[1]*a1[0];
    }
  }
}

static void pos_2_tsc_serial(ulint np,float *pos,float *delta)
{
  if(NodeThis==0)
    printf("  Calculating TSC\n");
  memset(delta,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));

  lint ii;
  float i_agrid=Ngrid/Lbox;
  for(ii=0;ii<np;ii++) {
    int ax,i0[3],ip[3],im[3];
    float a0[3],ap[3],am[3];
    float *p=&(pos[3*ii]);
    
    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(p[ax]*i_agrid+0.5);
      a0[ax]=p[ax]*i_agrid-i0[ax];
      am[ax]=0.5*(0.5-a0[ax])*(0.5-a0[ax]);
      ap[ax]=0.5*(0.5+a0[ax])*(0.5+a0[ax]);
      a0[ax]=0.75-a0[ax]*a0[ax];
      ip[ax]=i0[ax]+1;
      im[ax]=i0[ax]-1;
      if(im[ax]<0) im[ax]+=Ngrid;
      if(i0[ax]<0) i0[ax]+=Ngrid;
      if(ip[ax]<0) ip[ax]+=Ngrid;
      if(im[ax]>=Ngrid) im[ax]-=Ngrid;
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      if(ip[ax]>=Ngrid) ip[ax]-=Ngrid;
    }
    
    delta[im[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*im[0]))]+=am[2]*am[1]*am[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*im[0]))]+=a0[2]*am[1]*am[0];
    delta[ip[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*im[0]))]+=ap[2]*am[1]*am[0];
    delta[im[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*im[0]))]+=am[2]*a0[1]*am[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*im[0]))]+=a0[2]*a0[1]*am[0];
    delta[ip[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*im[0]))]+=ap[2]*a0[1]*am[0];
    delta[im[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*im[0]))]+=am[2]*ap[1]*am[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*im[0]))]+=a0[2]*ap[1]*am[0];
    delta[ip[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*im[0]))]+=ap[2]*ap[1]*am[0];
    delta[im[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*i0[0]))]+=am[2]*am[1]*a0[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*i0[0]))]+=a0[2]*am[1]*a0[0];
    delta[ip[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*i0[0]))]+=ap[2]*am[1]*a0[0];
    delta[im[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]))]+=am[2]*a0[1]*a0[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]))]+=a0[2]*a0[1]*a0[0];
    delta[ip[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]))]+=ap[2]*a0[1]*a0[0];
    delta[im[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*i0[0]))]+=am[2]*ap[1]*a0[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*i0[0]))]+=a0[2]*ap[1]*a0[0];
    delta[ip[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*i0[0]))]+=ap[2]*ap[1]*a0[0];
    delta[im[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*ip[0]))]+=am[2]*am[1]*ap[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*ip[0]))]+=a0[2]*am[1]*ap[0];
    delta[ip[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*ip[0]))]+=ap[2]*am[1]*ap[0];
    delta[im[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*ip[0]))]+=am[2]*a0[1]*ap[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*ip[0]))]+=a0[2]*a0[1]*ap[0];
    delta[ip[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*ip[0]))]+=ap[2]*a0[1]*ap[0];
    delta[im[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*ip[0]))]+=am[2]*ap[1]*ap[0];
    delta[i0[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*ip[0]))]+=a0[2]*ap[1]*ap[0];
    delta[ip[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*ip[0]))]+=ap[2]*ap[1]*ap[0];
  }
}

static void pos_2_tsc_parallel(ulint np,float *pos,float *delta)
{
  if(NodeThis==0)
    printf("  Calculating TSC\n");
  memset(delta,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));

  lint ii;
  float i_agrid=Ngrid/Lbox;
  for(ii=0;ii<np;ii++) {
    int ax,i0[3],ip[3],im[3];
    float a0[3],ap[3],am[3];
    float *p=&(pos[3*ii]);
    
    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(floorf(p[ax]*i_agrid+0.5));
      a0[ax]=p[ax]*i_agrid-i0[ax];
      am[ax]=0.5*(0.5-a0[ax])*(0.5-a0[ax]);
      ap[ax]=0.5*(0.5+a0[ax])*(0.5+a0[ax]);
      a0[ax]=0.75-a0[ax]*a0[ax];
      ip[ax]=i0[ax]+1;
      im[ax]=i0[ax]-1;
    }
    for(ax=1;ax<3;ax++) {
      if(im[ax]<0) im[ax]+=Ngrid;
      if(i0[ax]<0) i0[ax]+=Ngrid;
      if(ip[ax]<0) ip[ax]+=Ngrid;
      if(im[ax]>=Ngrid) im[ax]-=Ngrid;
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      if(ip[ax]>=Ngrid) ip[ax]-=Ngrid;
    }
    im[0]-=Ix0_here;
    i0[0]-=Ix0_here;
    ip[0]-=Ix0_here;

    if(im[0]>=0 && im[0]<Nx_here) {
      delta[im[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*im[0]))]+=am[2]*am[1]*am[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*im[0]))]+=a0[2]*am[1]*am[0];
      delta[ip[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*im[0]))]+=ap[2]*am[1]*am[0];
      delta[im[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*im[0]))]+=am[2]*a0[1]*am[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*im[0]))]+=a0[2]*a0[1]*am[0];
      delta[ip[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*im[0]))]+=ap[2]*a0[1]*am[0];
      delta[im[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*im[0]))]+=am[2]*ap[1]*am[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*im[0]))]+=a0[2]*ap[1]*am[0];
      delta[ip[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*im[0]))]+=ap[2]*ap[1]*am[0];
    }    
    if(i0[0]>=0 && i0[0]<Nx_here) {
      delta[im[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*i0[0]))]+=am[2]*am[1]*a0[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*i0[0]))]+=a0[2]*am[1]*a0[0];
      delta[ip[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*i0[0]))]+=ap[2]*am[1]*a0[0];
      delta[im[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]))]+=am[2]*a0[1]*a0[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]))]+=a0[2]*a0[1]*a0[0];
      delta[ip[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*i0[0]))]+=ap[2]*a0[1]*a0[0];
      delta[im[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*i0[0]))]+=am[2]*ap[1]*a0[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*i0[0]))]+=a0[2]*ap[1]*a0[0];
      delta[ip[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*i0[0]))]+=ap[2]*ap[1]*a0[0];
    }
    if(ip[0]>=0 && ip[0]<Nx_here) {
      delta[im[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*ip[0]))]+=am[2]*am[1]*ap[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*ip[0]))]+=a0[2]*am[1]*ap[0];
      delta[ip[2]+2*(Ngrid/2+1)*((lint)(im[1]+Ngrid*ip[0]))]+=ap[2]*am[1]*ap[0];
      delta[im[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*ip[0]))]+=am[2]*a0[1]*ap[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*ip[0]))]+=a0[2]*a0[1]*ap[0];
      delta[ip[2]+2*(Ngrid/2+1)*((lint)(i0[1]+Ngrid*ip[0]))]+=ap[2]*a0[1]*ap[0];
      delta[im[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*ip[0]))]+=am[2]*ap[1]*ap[0];
      delta[i0[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*ip[0]))]+=a0[2]*ap[1]*ap[0];
      delta[ip[2]+2*(Ngrid/2+1)*((lint)(ip[1]+Ngrid*ip[0]))]+=ap[2]*ap[1]*ap[0];
    }
  }
}

void compute_overdensity(ulint np,ulint np_total,float *pos,float *delta,int interp_order)
{
  if(interp_order==0) {
    if(NNodes==1)
      pos_2_ngp_serial(np,pos,delta);
    else
      pos_2_ngp_parallel(np,pos,delta);
  }
  else if(interp_order==1) {
    if(NNodes==1)
      pos_2_cic_serial(np,pos,delta);
    else
      pos_2_cic_parallel(np,pos,delta);
  }
  else if(interp_order==2) {
    if(NNodes==1)
      pos_2_tsc_serial(np,pos,delta);
    else
      pos_2_tsc_parallel(np,pos,delta);
  }
  else
    report_error(1,"Wrong interpolation order\n");

  //Compute overdensity from particle number density
  int ix;
  lint ncells=Ngrid*((lint)(Ngrid*Ngrid));
  float i_densmean=(float)ncells/np_total;
  for(ix=0;ix<Nx_here;ix++) {
    int iy;
    for(iy=0;iy<Ngrid;iy++) {
      int iz;
      for(iz=0;iz<Ngrid;iz++) {
	lint index=iz+2*(Ngrid/2+1)*((lint)(iy+Ngrid*ix));
	delta[index]=delta[index]*i_densmean-1;
      }
    }
  }
}

static inline void add_weight(int ix,int iy,int iz,float w,float *dg,float **vg,float *v)
{
  lint index=iz+2*(Ngrid/2+1)*((lint)(iy+Ngrid*ix));
  dg[index]+=w;
  vg[0][index]+=w*v[0];
  vg[1][index]+=w*v[1];
  vg[2][index]+=w*v[2];
}

static void vel_2_ngp_serial(ulint np,float *pos,float *vel,float *densgrid,float **velgrid)
{
  if(NodeThis==0)
    printf("  Calculating NGP\n");
  memset(densgrid,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[0],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[1],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[2],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  
  lint ii;
  float i_agrid=1./Ngrid;

  for(ii=0;ii<np;ii++) {
    int ax,i0[3];
    float *p=&(pos[3*ii]);
    float *v=&(vel[3*ii]);

    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(p[ax]*i_agrid+0.5);
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      else if(i0[ax]<0) i0[ax]+=Ngrid;
    }

    add_weight(i0[0],i0[1],i0[2],1.,densgrid,velgrid,v);
  }
}

static void vel_2_ngp_parallel(ulint np,float *pos,float *vel,float *densgrid,float **velgrid)
{
  lint ii;
  float i_agrid=1./Ngrid;

  if(NodeThis==0)
    printf("  Calculating NGP\n");
  memset(densgrid,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[0],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[1],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[2],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  
  for(ii=0;ii<np;ii++) {
    int ax,i0[3];
    float *p=&(pos[3*ii]);
    float *v=&(vel[3*ii]);

    for(ax=0;ax<3;ax++)
      i0[ax]=(int)(p[ax]*i_agrid+0.5);
    if(i0[1]>=Ngrid) i0[1]-=Ngrid;
    if(i0[1]<0) i0[1]+=Ngrid;
    if(i0[2]>=Ngrid) i0[2]-=Ngrid;
    if(i0[2]<0) i0[2]+=Ngrid;

    if(i0[0]>=0 && i0[0]<Nx_here) {
      add_weight(i0[0],i0[1],i0[2],1.,densgrid,velgrid,v);
    }
  }
}

static void vel_2_cic_serial(ulint np,float *pos,float *vel,float *densgrid,float **velgrid)
{
  if(NodeThis==0)
    printf("  Calculating CIC\n");
  memset(densgrid,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[0],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[1],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[2],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));

  lint ii;
  float i_agrid=Ngrid/Lbox;

  for(ii=0;ii<np;ii++) {
    int ax,i0[3],i1[3];
    float a0[3],a1[3];
    float *p=&(pos[3*ii]);
    float *v=&(vel[3*ii]);

    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(p[ax]*i_agrid);
      a1[ax]=p[ax]*i_agrid-i0[ax];
      a0[ax]=1-a1[ax];
      i1[ax]=i0[ax]+1;
      if(i0[ax]<0) i0[ax]+=Ngrid;
      if(i1[ax]<0) i1[ax]+=Ngrid;
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      if(i1[ax]>=Ngrid) i1[ax]-=Ngrid;
    }
    add_weight(i0[0],i0[1],i0[2],a0[2]*a0[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],i0[1],i1[2],a1[2]*a0[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],i1[1],i0[2],a0[2]*a1[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],i1[1],i1[2],a1[2]*a1[1]*a0[0],densgrid,velgrid,v);
    add_weight(i1[0],i0[1],i0[2],a0[2]*a0[1]*a1[0],densgrid,velgrid,v);
    add_weight(i1[0],i0[1],i1[2],a1[2]*a0[1]*a1[0],densgrid,velgrid,v);
    add_weight(i1[0],i1[1],i0[2],a0[2]*a1[1]*a1[0],densgrid,velgrid,v);
    add_weight(i1[0],i1[1],i1[2],a1[2]*a1[1]*a1[0],densgrid,velgrid,v);
  }
}

static void vel_2_cic_parallel(ulint np,float *pos,float *vel,float *densgrid,float **velgrid)
{
  if(NodeThis==0)
    printf("  Calculating CIC\n");
  memset(densgrid,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[0],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[1],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[2],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));

  lint ii;
  float i_agrid=Ngrid/Lbox;
  for(ii=0;ii<np;ii++) {
    int ax,i0[3],i1[3];
    float a0[3],a1[3];
    float *p=&(pos[3*ii]);
    float *v=&(vel[3*ii]);

    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(floorf(p[ax]*i_agrid));
      a1[ax]=p[ax]*i_agrid-i0[ax];
      a0[ax]=1-a1[ax];
      i1[ax]=i0[ax]+1;
    }
    for(ax=1;ax<3;ax++) {
      if(i0[ax]<0) i0[ax]+=Ngrid;
      if(i1[ax]<0) i1[ax]+=Ngrid;
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      if(i1[ax]>=Ngrid) i1[ax]-=Ngrid;
    }
    i0[0]-=Ix0_here;
    i1[0]-=Ix0_here;

    if(i0[0]>=0 && i0[0]<Nx_here) {
      add_weight(i0[0],i0[1],i0[2],a0[2]*a0[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],i0[1],i1[2],a1[2]*a0[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],i1[1],i0[2],a0[2]*a1[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],i1[1],i1[2],a1[2]*a1[1]*a0[0],densgrid,velgrid,v);
    }
    if(i1[0]>=0 && i1[0]<Nx_here) {
      add_weight(i1[0],i0[1],i0[2],a0[2]*a0[1]*a1[0],densgrid,velgrid,v);
      add_weight(i1[0],i0[1],i1[2],a1[2]*a0[1]*a1[0],densgrid,velgrid,v);
      add_weight(i1[0],i1[1],i0[2],a0[2]*a1[1]*a1[0],densgrid,velgrid,v);
      add_weight(i1[0],i1[1],i1[2],a1[2]*a1[1]*a1[0],densgrid,velgrid,v);
    }
  }
}

static void vel_2_tsc_serial(ulint np,float *pos,float *vel,float *densgrid,float **velgrid)
{
  if(NodeThis==0)
    printf("  Calculating TSC\n");
  memset(densgrid,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[0],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[1],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[2],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));

  lint ii;
  float i_agrid=Ngrid/Lbox;
  for(ii=0;ii<np;ii++) {
    int ax,i0[3],ip[3],im[3];
    float a0[3],ap[3],am[3];
    float *p=&(pos[3*ii]);
    float *v=&(vel[3*ii]);
    
    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(p[ax]*i_agrid+0.5);
      a0[ax]=p[ax]*i_agrid-i0[ax];
      am[ax]=0.5*(0.5-a0[ax])*(0.5-a0[ax]);
      ap[ax]=0.5*(0.5+a0[ax])*(0.5+a0[ax]);
      a0[ax]=0.75-a0[ax]*a0[ax];
      ip[ax]=i0[ax]+1;
      im[ax]=i0[ax]-1;
      if(im[ax]<0) im[ax]+=Ngrid;
      if(i0[ax]<0) i0[ax]+=Ngrid;
      if(ip[ax]<0) ip[ax]+=Ngrid;
      if(im[ax]>=Ngrid) im[ax]-=Ngrid;
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      if(ip[ax]>=Ngrid) ip[ax]-=Ngrid;
    }
    
    add_weight(im[0],im[1],im[2],am[2]*am[1]*am[0],densgrid,velgrid,v);
    add_weight(im[0],im[1],i0[2],a0[2]*am[1]*am[0],densgrid,velgrid,v);
    add_weight(im[0],im[1],ip[2],ap[2]*am[1]*am[0],densgrid,velgrid,v);
    add_weight(im[0],i0[1],im[2],am[2]*a0[1]*am[0],densgrid,velgrid,v);
    add_weight(im[0],i0[1],i0[2],a0[2]*a0[1]*am[0],densgrid,velgrid,v);
    add_weight(im[0],i0[1],ip[2],ap[2]*a0[1]*am[0],densgrid,velgrid,v);
    add_weight(im[0],ip[1],im[2],am[2]*ap[1]*am[0],densgrid,velgrid,v);
    add_weight(im[0],ip[1],i0[2],a0[2]*ap[1]*am[0],densgrid,velgrid,v);
    add_weight(im[0],ip[1],ip[2],ap[2]*ap[1]*am[0],densgrid,velgrid,v);
    add_weight(i0[0],im[1],im[2],am[2]*am[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],im[1],i0[2],a0[2]*am[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],im[1],ip[2],ap[2]*am[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],i0[1],im[2],am[2]*a0[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],i0[1],i0[2],a0[2]*a0[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],i0[1],ip[2],ap[2]*a0[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],ip[1],im[2],am[2]*ap[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],ip[1],i0[2],a0[2]*ap[1]*a0[0],densgrid,velgrid,v);
    add_weight(i0[0],ip[1],ip[2],ap[2]*ap[1]*a0[0],densgrid,velgrid,v);
    add_weight(ip[0],im[1],im[2],am[2]*am[1]*ap[0],densgrid,velgrid,v);
    add_weight(ip[0],im[1],i0[2],a0[2]*am[1]*ap[0],densgrid,velgrid,v);
    add_weight(ip[0],im[1],ip[2],ap[2]*am[1]*ap[0],densgrid,velgrid,v);
    add_weight(ip[0],i0[1],im[2],am[2]*a0[1]*ap[0],densgrid,velgrid,v);
    add_weight(ip[0],i0[1],i0[2],a0[2]*a0[1]*ap[0],densgrid,velgrid,v);
    add_weight(ip[0],i0[1],ip[2],ap[2]*a0[1]*ap[0],densgrid,velgrid,v);
    add_weight(ip[0],ip[1],im[2],am[2]*ap[1]*ap[0],densgrid,velgrid,v);
    add_weight(ip[0],ip[1],i0[2],a0[2]*ap[1]*ap[0],densgrid,velgrid,v);
    add_weight(ip[0],ip[1],ip[2],ap[2]*ap[1]*ap[0],densgrid,velgrid,v);
  }
}

static void vel_2_tsc_parallel(ulint np,float *pos,float *vel,float *densgrid,float **velgrid)
{
  if(NodeThis==0)
    printf("  Calculating TSC\n");
  memset(densgrid,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[0],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[1],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  memset(velgrid[2],0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));

  lint ii;
  float i_agrid=Ngrid/Lbox;
  for(ii=0;ii<np;ii++) {
    int ax,i0[3],ip[3],im[3];
    float a0[3],ap[3],am[3];
    float *p=&(pos[3*ii]);
    float *v=&(vel[3*ii]);
    
    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(floorf(p[ax]*i_agrid+0.5));
      a0[ax]=p[ax]*i_agrid-i0[ax];
      am[ax]=0.5*(0.5-a0[ax])*(0.5-a0[ax]);
      ap[ax]=0.5*(0.5+a0[ax])*(0.5+a0[ax]);
      a0[ax]=0.75-a0[ax]*a0[ax];
      ip[ax]=i0[ax]+1;
      im[ax]=i0[ax]-1;
    }
    for(ax=1;ax<3;ax++) {
      if(im[ax]<0) im[ax]+=Ngrid;
      if(i0[ax]<0) i0[ax]+=Ngrid;
      if(ip[ax]<0) ip[ax]+=Ngrid;
      if(im[ax]>=Ngrid) im[ax]-=Ngrid;
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      if(ip[ax]>=Ngrid) ip[ax]-=Ngrid;
    }
    im[0]-=Ix0_here;
    i0[0]-=Ix0_here;
    ip[0]-=Ix0_here;

    if(im[0]>=0 && im[0]<Nx_here) {
      add_weight(im[0],im[1],im[2],am[2]*am[1]*am[0],densgrid,velgrid,v);
      add_weight(im[0],im[1],i0[2],a0[2]*am[1]*am[0],densgrid,velgrid,v);
      add_weight(im[0],im[1],ip[2],ap[2]*am[1]*am[0],densgrid,velgrid,v);
      add_weight(im[0],i0[1],im[2],am[2]*a0[1]*am[0],densgrid,velgrid,v);
      add_weight(im[0],i0[1],i0[2],a0[2]*a0[1]*am[0],densgrid,velgrid,v);
      add_weight(im[0],i0[1],ip[2],ap[2]*a0[1]*am[0],densgrid,velgrid,v);
      add_weight(im[0],ip[1],im[2],am[2]*ap[1]*am[0],densgrid,velgrid,v);
      add_weight(im[0],ip[1],i0[2],a0[2]*ap[1]*am[0],densgrid,velgrid,v);
      add_weight(im[0],ip[1],ip[2],ap[2]*ap[1]*am[0],densgrid,velgrid,v);
    }    
    if(i0[0]>=0 && i0[0]<Nx_here) {
      add_weight(i0[0],im[1],im[2],am[2]*am[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],im[1],i0[2],a0[2]*am[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],im[1],ip[2],ap[2]*am[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],i0[1],im[2],am[2]*a0[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],i0[1],i0[2],a0[2]*a0[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],i0[1],ip[2],ap[2]*a0[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],ip[1],im[2],am[2]*ap[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],ip[1],i0[2],a0[2]*ap[1]*a0[0],densgrid,velgrid,v);
      add_weight(i0[0],ip[1],ip[2],ap[2]*ap[1]*a0[0],densgrid,velgrid,v);
    }
    if(ip[0]>=0 && ip[0]<Nx_here) {
      add_weight(ip[0],im[1],im[2],am[2]*am[1]*ap[0],densgrid,velgrid,v);
      add_weight(ip[0],im[1],i0[2],a0[2]*am[1]*ap[0],densgrid,velgrid,v);
      add_weight(ip[0],im[1],ip[2],ap[2]*am[1]*ap[0],densgrid,velgrid,v);
      add_weight(ip[0],i0[1],im[2],am[2]*a0[1]*ap[0],densgrid,velgrid,v);
      add_weight(ip[0],i0[1],i0[2],a0[2]*a0[1]*ap[0],densgrid,velgrid,v);
      add_weight(ip[0],i0[1],ip[2],ap[2]*a0[1]*ap[0],densgrid,velgrid,v);
      add_weight(ip[0],ip[1],im[2],am[2]*ap[1]*ap[0],densgrid,velgrid,v);
      add_weight(ip[0],ip[1],i0[2],a0[2]*ap[1]*ap[0],densgrid,velgrid,v);
      add_weight(ip[0],ip[1],ip[2],ap[2]*ap[1]*ap[0],densgrid,velgrid,v);
    }
  }
}

void compute_velocity_and_overdensity(ulint np,ulint np_total,float *pos,float *vel,
				      float *delta,float **velgrid,int interp_order)
{
  if(interp_order==0) {
    if(NNodes==1)
      vel_2_ngp_serial(np,pos,vel,delta,velgrid);
    else
      vel_2_ngp_parallel(np,pos,vel,delta,velgrid);
  }
  else if(interp_order==1) {
    if(NNodes==1)
      vel_2_cic_serial(np,pos,vel,delta,velgrid);
    else
      vel_2_cic_parallel(np,pos,vel,delta,velgrid);
  }
  else if(interp_order==2) {
    if(NNodes==1)
      vel_2_tsc_serial(np,pos,vel,delta,velgrid);
    else
      vel_2_tsc_parallel(np,pos,vel,delta,velgrid);
  }
  else
    report_error(1,"Wrong interpolation order\n");

  //Compute overdensity from particle number density
  int ix;
  lint ncells=Ngrid*((lint)(Ngrid*Ngrid));
  float i_densmean=(float)ncells/np_total;
  for(ix=0;ix<Nx_here;ix++) {
    int iy;
    for(iy=0;iy<Ngrid;iy++) {
      int iz;
      for(iz=0;iz<Ngrid;iz++) {
	lint index=iz+2*(Ngrid/2+1)*((lint)(iy+Ngrid*ix));
	if(delta[index]>0) {
	  velgrid[0][index]/=delta[index];
	  velgrid[1][index]/=delta[index];
	  velgrid[2][index]/=delta[index];
	}
	delta[index]=delta[index]*i_densmean-1;
      }
    }
  }
}

void smooth_density_fourier(float r_smooth)
{
  int iy;
  fftwf_plan plan_tof,plan_tor;
  float x_smooth2=(float)(2*M_PI*r_smooth/Lbox);
  float norm=(float)(1./Ngrid);
  x_smooth2=x_smooth2*x_smooth2;
  norm=norm*norm*norm;

#ifdef _DEBUG
  printf("Node %d Planning\n",NodeThis);
#endif //_DEBUG
  plan_tof=fftwf_mpi_plan_dft_r2c_3d(Ngrid,Ngrid,Ngrid,Dens_local,Cdens_local,
				     MPI_COMM_WORLD,FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT);
  plan_tor=fftwf_mpi_plan_dft_c2r_3d(Ngrid,Ngrid,Ngrid,Cdens_local,Dens_local,
  				     MPI_COMM_WORLD,FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
#ifdef _DEBUG
  printf("Node %d Transforming\n",NodeThis);
#endif //_DEBUG
  fftwf_execute(plan_tof);
  fftwf_destroy_plan(plan_tof);

#ifdef _DEBUG
  printf("Node %d Smoothing\n",NodeThis);
#endif //_DEBUG
  for(iy=0;iy<Ny_here;iy++) {
    int ix;
    int iy_true=iy+Iy0_here;
    int iy0= iy_true<=(Ngrid/2) ? iy_true : iy_true-Ngrid;
    for(ix=0;ix<Ngrid;ix++) {
      int iz;
      int ix0= ix<=(Ngrid/2) ? ix : ix-Ngrid;
      for(iz=0;iz<Ngrid/2+1;iz++) {
	float x2=x_smooth2*(ix0*ix0+iy0*iy0+iz*iz);
	float sm=(float)(exp(-0.5*x2));
	long index=iz+(Ngrid/2+1)*((long)(ix+Ngrid*iy));

	Cdens_local[index]*=norm;
	Cdens_sm_local[index]=Cdens_local[index]*sm;
      }
    }
  }

#ifdef _DEBUG
  printf("Node %d Transforming back\n",NodeThis);
#endif //_DEBUG
  fftwf_execute(plan_tor);
  fftwf_destroy_plan(plan_tor);
}

void get_smoothed_density_real(void)
{
  fftwf_plan plan_tor;

#ifdef _DEBUG
  printf("Node %d Planning\n",NodeThis);
#endif //_DEBUG
  plan_tor=fftwf_mpi_plan_dft_c2r_3d(Ngrid,Ngrid,Ngrid,Cdens_sm_local,Dens_sm_local,
  				     MPI_COMM_WORLD,FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
#ifdef _DEBUG
  printf("Node %d Transforming back\n",NodeThis);
#endif //_DEBUG
  fftwf_execute(plan_tor);
  fftwf_destroy_plan(plan_tor);
}

void get_tidal_field(void)
{
  int iy;
  fftwf_plan plan_t_tor[6];

#ifdef _DEBUG
  printf("Node %d Planning\n",NodeThis);
#endif //_DEBUG
  for(iy=0;iy<6;iy++) {
    plan_t_tor[iy]=fftwf_mpi_plan_dft_c2r_3d(Ngrid,Ngrid,Ngrid,Ctid_local[iy],Tid_local[iy],
					     MPI_COMM_WORLD,FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
  }

#ifdef _DEBUG
  printf("Node %d computing tidal field\n",NodeThis);
#endif //_DEBUG
  for(iy=0;iy<Ny_here;iy++) {
    int ix;
    int iy_true=iy+Iy0_here;
    int iy0= iy_true<=(Ngrid/2) ? iy_true : iy_true-Ngrid;
    for(ix=0;ix<Ngrid;ix++) {
      int iz;
      int ix0= ix<=(Ngrid/2) ? ix : ix-Ngrid;
      for(iz=0;iz<Ngrid/2+1;iz++) {
	float i_kmod2=0;
	float kmod2=(float)(ix0*ix0+iy0*iy0+iz*iz);
	long index=iz+(Ngrid/2+1)*((long)(ix+Ngrid*iy));

	if(kmod2>0)
	  i_kmod2=1./kmod2;

	Ctid_local[0][index]=Cdens_sm_local[index]*ix0*ix0*i_kmod2; //xx
	Ctid_local[1][index]=Cdens_sm_local[index]*iy0*iy0*i_kmod2; //yy
	Ctid_local[2][index]=Cdens_sm_local[index]*iz *iz *i_kmod2; //zz
	Ctid_local[3][index]=Cdens_sm_local[index]*ix0*iy0*i_kmod2; //xy
	Ctid_local[4][index]=Cdens_sm_local[index]*iy0*iz *i_kmod2; //yz
	Ctid_local[5][index]=Cdens_sm_local[index]*iz *ix0*i_kmod2; //zx
      }
    }
  }

#ifdef _DEBUG
  printf("Node %d Transforming back\n",NodeThis);
#endif //_DEBUG
  for(iy=0;iy<6;iy++) {
    fftwf_execute(plan_t_tor[iy]);
    fftwf_destroy_plan(plan_t_tor[iy]);
  }
}

void get_linearized_velocity(gad_header head)
{
  //TODO correct normalization

  int iy;
  fftwf_plan plan_v_tor[3];
  float dk=2*M_PI/Lbox;

  double a=head.time;
  double hub=sqrt(head.Omega0/(a*a*a)+head.OmegaLambda+(1-head.Omega0-head.OmegaLambda)/(a*a));
  double omega_m=head.Omega0/(head.Omega0+head.OmegaLambda*a*a*a+
			      (1-head.Omega0-head.OmegaLambda)*a);
  float prefac_vel=100*pow(omega_m,0.55)*hub;

#ifdef _DEBUG
  printf("Node %d Planning\n",NodeThis);
#endif //_DEBUG
  for(iy=0;iy<3;iy++) {
    plan_v_tor[iy]=fftwf_mpi_plan_dft_c2r_3d(Ngrid,Ngrid,Ngrid,Clvel_local[iy],Lvel_local[iy],
					    MPI_COMM_WORLD,FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
  }

#ifdef _DEBUG
  printf("Node %d computing linearized velocity\n",NodeThis);
#endif //_DEBUG
  for(iy=0;iy<Ny_here;iy++) {
    int ix;
    int iy_true=iy+Iy0_here;
    int iy0= iy_true<=(Ngrid/2) ? iy_true : iy_true-Ngrid;
    for(ix=0;ix<Ngrid;ix++) {
      int iz;
      int ix0= ix<=(Ngrid/2) ? ix : ix-Ngrid;
      for(iz=0;iz<Ngrid/2+1;iz++) {
	float i_kmod2=0;
	float kmod2=(float)(ix0*ix0+iy0*iy0+iz*iz);
	long index=iz+(Ngrid/2+1)*((long)(ix+Ngrid*iy));

	if(kmod2>0)
	  i_kmod2=1./(dk*kmod2);

	Clvel_local[0][index]=I*prefac_vel*Cdens_sm_local[index]*ix0*i_kmod2; //vx
	Clvel_local[1][index]=I*prefac_vel*Cdens_sm_local[index]*iy0*i_kmod2; //vy
	Clvel_local[2][index]=I*prefac_vel*Cdens_sm_local[index]*iz *i_kmod2; //vz
      }
    }
  }

#ifdef _DEBUG
  printf("Node %d Transforming back\n",NodeThis);
#endif //_DEBUG
  for(iy=0;iy<3;iy++) {
    fftwf_execute(plan_v_tor[iy]);
    fftwf_destroy_plan(plan_v_tor[iy]);
  }
}


/*
typedef struct {
  double eval;
  double evec[3];
} Esys;

int compare_evals(const void *p1,const void *p2)
{
  Esys *e1=(Esys *)p1;
  Esys *e2=(Esys *)p2;

  return (e1->eval>e2->eval) - (e1->eval<e2->eval);
}

void decompose_tidal(float *tidal,float *delta,float *delta_sm,float *tidal_info)
{
  lint ngtot=((lint)(Ngrid*Ngrid))*Ngrid;

#pragma omp parallel default(none)		\
  shared(tidal,Ngrid,tidal_info,delta_sm,delta,ngtot)
  {
    int im;
    lint index;
    gsl_vector *eval=gsl_vector_alloc(3);
    gsl_matrix *tij=gsl_matrix_alloc(3,3);
    gsl_matrix *evec=gsl_matrix_alloc(3,3);
    gsl_eigen_symmv_workspace *w=gsl_eigen_symmv_alloc(3);

#pragma omp for
    for(index=0;index<ngtot;index++) {
      Esys leig[3];

      gsl_matrix_set(tij,0,0,tidal[0+6*index]);
      gsl_matrix_set(tij,1,1,tidal[1+6*index]);
      gsl_matrix_set(tij,2,2,tidal[2+6*index]);
      gsl_matrix_set(tij,0,1,tidal[3+6*index]);
      gsl_matrix_set(tij,1,0,tidal[3+6*index]);
      gsl_matrix_set(tij,0,2,tidal[4+6*index]);
      gsl_matrix_set(tij,2,0,tidal[4+6*index]);
      gsl_matrix_set(tij,1,2,tidal[5+6*index]);
      gsl_matrix_set(tij,2,1,tidal[5+6*index]);
      
      gsl_eigen_symmv(tij,eval,evec,w);
      
      for(im=0;im<3;im++) {
	leig[im].eval=gsl_vector_get(eval,im);
	leig[im].evec[0]=gsl_matrix_get(evec,0,im);
	leig[im].evec[1]=gsl_matrix_get(evec,1,im);
	leig[im].evec[2]=gsl_matrix_get(evec,2,im);
      }

      qsort(leig,3,sizeof(Esys),compare_evals);

#ifdef _DEBUG
      gsl_matrix *evec_here=gsl_matrix_alloc(3,3);
      gsl_matrix *tij_here=gsl_matrix_alloc(3,3);
      gsl_matrix *aux_here1=gsl_matrix_alloc(3,3);
      gsl_matrix *aux_here2=gsl_matrix_alloc(3,3);

      gsl_matrix_set_zero(evec_here);
      gsl_matrix_set_zero(tij_here);
      gsl_matrix_set_zero(aux_here1);
      gsl_matrix_set_zero(aux_here2);

      gsl_matrix_set(tij_here,0,0,tidal[0+6*index]);
      gsl_matrix_set(tij_here,1,1,tidal[1+6*index]);
      gsl_matrix_set(tij_here,2,2,tidal[2+6*index]);
      gsl_matrix_set(tij_here,0,1,tidal[3+6*index]);
      gsl_matrix_set(tij_here,1,0,tidal[3+6*index]);
      gsl_matrix_set(tij_here,0,2,tidal[4+6*index]);
      gsl_matrix_set(tij_here,2,0,tidal[4+6*index]);
      gsl_matrix_set(tij_here,1,2,tidal[5+6*index]);
      gsl_matrix_set(tij_here,2,1,tidal[5+6*index]);
      
      for(im=0;im<3;im++) {
	int jm;
	for(jm=0;jm<3;jm++)
	  gsl_matrix_set(evec_here,im,jm,leig[jm].evec[im]);
      }

      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,tij_here,evec_here,0,aux_here1);
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,evec_here,aux_here1,0,aux_here2);

      for(im=0;im<3;im++)
	printf("%lE ",leig[im].eval);
      printf("\n-\n");
      for(im=0;im<3;im++) {
	int jm;
	for(jm=0;jm<3;jm++)
	  printf("%lE ",gsl_matrix_get(aux_here2,im,jm));
	printf("\n");
      }
      gsl_matrix_free(tij_here);
      gsl_matrix_free(evec_here);
      gsl_matrix_free(aux_here1);
      gsl_matrix_free(aux_here2);
      scanf("%d",&im);
#endif //_DEBUG
      
      tidal_info[14*index+0] =delta[index];
      tidal_info[14*index+1] =delta_sm[index];
      tidal_info[14*index+2] =leig[0].eval;
      tidal_info[14*index+3] =leig[0].evec[0];
      tidal_info[14*index+4] =leig[0].evec[1];
      tidal_info[14*index+5] =leig[0].evec[2];

      tidal_info[14*index+6] =leig[1].eval;
      tidal_info[14*index+7] =leig[1].evec[0];
      tidal_info[14*index+8] =leig[1].evec[1];
      tidal_info[14*index+9] =leig[1].evec[2];

      tidal_info[14*index+10]=leig[2].eval;
      tidal_info[14*index+11]=leig[2].evec[0];
      tidal_info[14*index+12]=leig[2].evec[1];
      tidal_info[14*index+13]=leig[2].evec[2];
    } //end omp for

    gsl_eigen_symmv_free(w);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(tij);
  } //end omp parallel
}
*/
