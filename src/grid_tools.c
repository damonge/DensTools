#include "common.h"

static void pos_2_ngp_serial(ulint np,float *pos,float *delta)
{
  if(NodeThis==0)
    printf("  Calculating NGP\n");
  memset(delta,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  
  lint ii;
  float i_agrid=Ngrid/Lbox;

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
    delta[index]+=1.;
  }
}

static void pos_2_ngp_parallel(ulint np,float *pos,float *delta)
{
  lint ii;
  float i_agrid=Ngrid/Lbox;

  if(NodeThis==0)
    printf("  Calculating NGP\n");
  memset(delta,0,sizeof(float)*Nx_here*((lint)(Ngrid*Ngrid)));
  
  for(ii=0;ii<np;ii++) {
    lint index;
    int ax,i0[3];
    float *p=&(pos[3*ii]);

    for(ax=0;ax<3;ax++)
      i0[ax]=(int)(p[ax]*i_agrid+0.5);
    for(ax=1;ax<3;ax++) {
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      if(i0[ax]<0) i0[ax]+=Ngrid;
    }
    i0[0]-=Ix0_here;

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
      i0[ax]=(int)(p[ax]*i_agrid+0.5f);
      a0[ax]=p[ax]*i_agrid-i0[ax];
      am[ax]=0.5f*(0.5f-a0[ax])*(0.5f-a0[ax]);
      ap[ax]=0.5f*(0.5f+a0[ax])*(0.5f+a0[ax]);
      a0[ax]=0.75f-a0[ax]*a0[ax];
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
      i0[ax]=(int)(floorf(p[ax]*i_agrid+0.5f));
      a0[ax]=p[ax]*i_agrid-i0[ax];
      am[ax]=0.5f*(0.5f-a0[ax])*(0.5f-a0[ax]);
      ap[ax]=0.5f*(0.5f+a0[ax])*(0.5f+a0[ax]);
      a0[ax]=0.75f-a0[ax]*a0[ax];
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
  float i_agrid=Ngrid/Lbox;

  for(ii=0;ii<np;ii++) {
    int ax,i0[3];
    float *p=&(pos[3*ii]);
    float *v=&(vel[3*ii]);

    for(ax=0;ax<3;ax++) {
      i0[ax]=(int)(p[ax]*i_agrid+0.5f);
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      else if(i0[ax]<0) i0[ax]+=Ngrid;
    }

    add_weight(i0[0],i0[1],i0[2],1.0f,densgrid,velgrid,v);
  }
}

static void vel_2_ngp_parallel(ulint np,float *pos,float *vel,float *densgrid,float **velgrid)
{
  lint ii;
  float i_agrid=Ngrid/Lbox;

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
      i0[ax]=(int)(p[ax]*i_agrid+0.5f);
    for(ax=1;ax<3;ax++) {
      if(i0[ax]>=Ngrid) i0[ax]-=Ngrid;
      if(i0[ax]<0) i0[ax]+=Ngrid;
    }
    i0[0]-=Ix0_here;

    if(i0[0]>=0 && i0[0]<Nx_here) {
      add_weight(i0[0],i0[1],i0[2],1.0f,densgrid,velgrid,v);
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
      i0[ax]=(int)(p[ax]*i_agrid+0.5f);
      a0[ax]=p[ax]*i_agrid-i0[ax];
      am[ax]=0.5f*(0.5f-a0[ax])*(0.5f-a0[ax]);
      ap[ax]=0.5f*(0.5f+a0[ax])*(0.5f+a0[ax]);
      a0[ax]=0.75f-a0[ax]*a0[ax];
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
      i0[ax]=(int)(floorf(p[ax]*i_agrid+0.5f));
      a0[ax]=p[ax]*i_agrid-i0[ax];
      am[ax]=0.5f*(0.5f-a0[ax])*(0.5f-a0[ax]);
      ap[ax]=0.5f*(0.5f+a0[ax])*(0.5f+a0[ax]);
      a0[ax]=0.75f-a0[ax]*a0[ax];
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

static void get_tidal_field_fd(void)
{
  int iy;
  fftwf_plan plan_t_tor[6];
  fcomplex w=(fcomplex)(cexp(I*2*M_PI/Ngrid));

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
    fcomplex wy=(fcomplex)(cexp(I*2*M_PI*iy_true/Ngrid));
    fcomplex wx=1.0f;
    for(ix=0;ix<Ngrid;ix++) {
      int iz;
      fcomplex wz=1.0f;
      for(iz=0;iz<Ngrid/2+1;iz++) {
	long index=iz+(Ngrid/2+1)*((long)(ix+Ngrid*iy));
	if(index==0) {
	  Ctid_local[0][index]=0;
	  Ctid_local[1][index]=0;
	  Ctid_local[2][index]=0;
	  Ctid_local[3][index]=0;
	  Ctid_local[4][index]=0;
	  Ctid_local[5][index]=0;
	}
	else {
	  fcomplex i_kmod2=1.0f/(wx+1.0f/wx+wy+1.0f/wy+wz+1.0f/wz-6);
	  Ctid_local[0][index]=Cdens_sm_local[index]*(wx+1.0f/wx-2)*i_kmod2; //xx
	  Ctid_local[1][index]=Cdens_sm_local[index]*(wy+1.0f/wy-2)*i_kmod2; //yy
	  Ctid_local[2][index]=Cdens_sm_local[index]*(wz+1.0f/wz-2)*i_kmod2; //zz
	  Ctid_local[3][index]=Cdens_sm_local[index]*0.25f*(1.0f/(wx*wy)+wx*wy-wx/wy-wy/wx)*i_kmod2; //xy
	  Ctid_local[4][index]=Cdens_sm_local[index]*0.25f*(1.0f/(wy*wz)+wy*wz-wy/wz-wz/wy)*i_kmod2; //yz
	  Ctid_local[5][index]=Cdens_sm_local[index]*0.25f*(1.0f/(wz*wx)+wz*wx-wz/wx-wx/wz)*i_kmod2; //zx
	}

	wz*=w;
      }
      wx*=w;
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

static void get_tidal_field_ks(void)
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
	  i_kmod2=1.0f/kmod2;

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

void get_tidal_field(void)
{
  if(UseFD)
    get_tidal_field_fd();
  else
    get_tidal_field_ks();
}

static void get_linearized_velocity_fd(gad_header head)
{
  int iy;
  fftwf_plan plan_v_tor[3];
  float agrid=Lbox/Ngrid;
  fcomplex w=(fcomplex)(cexp(I*2*M_PI/Ngrid));

  double a=head.time;
  double hub=sqrt(head.Omega0/(a*a*a)+head.OmegaLambda+(1-head.Omega0-head.OmegaLambda)/(a*a));
  double omega_m=head.Omega0/(head.Omega0+head.OmegaLambda*a*a*a+
			      (1-head.Omega0-head.OmegaLambda)*a);
  float prefac_vel=(float)(100*pow(omega_m,0.55)*hub);

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
    fcomplex wy=(fcomplex)(cexp(I*2*M_PI*iy_true/Ngrid));
    fcomplex wx=1.0f;
    for(ix=0;ix<Ngrid;ix++) {
      int iz;
      fcomplex wz=1.0f;
      for(iz=0;iz<Ngrid/2+1;iz++) {
	long index=iz+(Ngrid/2+1)*((long)(ix+Ngrid*iy));
	if(index==0) {
	  Clvel_local[0][index]=0;
	  Clvel_local[1][index]=0;
	  Clvel_local[2][index]=0;
	}
	else {
	  fcomplex i_kmod2=1.0f/(wx+1.0f/wx+wy+1.0f/wy+wz+1.0f/wz-6);
	  Clvel_local[0][index]=prefac_vel*Cdens_sm_local[index]*agrid*0.5f*(1.0f/wx-wx)*i_kmod2;
	  Clvel_local[1][index]=prefac_vel*Cdens_sm_local[index]*agrid*0.5f*(1.0f/wy-wy)*i_kmod2;
	  Clvel_local[2][index]=prefac_vel*Cdens_sm_local[index]*agrid*0.5f*(1.0f/wz-wz)*i_kmod2;
	}
	wz*=w;
      }
      wx*=w;
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

static void get_linearized_velocity_ks(gad_header head)
{
  int iy;
  fftwf_plan plan_v_tor[3];
  float dk=(float)(2*M_PI/Lbox);

  double a=head.time;
  double hub=sqrt(head.Omega0/(a*a*a)+head.OmegaLambda+(1-head.Omega0-head.OmegaLambda)/(a*a));
  double omega_m=head.Omega0/(head.Omega0+head.OmegaLambda*a*a*a+
			      (1-head.Omega0-head.OmegaLambda)*a);
  float prefac_vel=(float)(100*pow(omega_m,0.55)*hub);

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
	  i_kmod2=(float)(1.0f/(dk*kmod2));

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

void get_linearized_velocity(gad_header head)
{
  if(UseFD)
    get_linearized_velocity_fd(head);
  else
    get_linearized_velocity_ks(head);
}

static void get_linearized_vpot_ks(gad_header head)
{
  int iy;
  fftwf_plan plan_p_tor;
  float dk=(float)(2*M_PI/Lbox);

  double a=head.time;
  double hub=sqrt(head.Omega0/(a*a*a)+head.OmegaLambda+(1-head.Omega0-head.OmegaLambda)/(a*a));
  double omega_m=head.Omega0/(head.Omega0+head.OmegaLambda*a*a*a+
			      (1-head.Omega0-head.OmegaLambda)*a);
  float prefac_vel=(float)(100*pow(omega_m,0.55)*hub);

#ifdef _DEBUG
  printf("Node %d Planning\n",NodeThis);
#endif //_DEBUG
  plan_p_tor=fftwf_mpi_plan_dft_c2r_3d(Ngrid,Ngrid,Ngrid,Cnlpot_local,Nlpot_local,
				       MPI_COMM_WORLD,FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);

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
	  i_kmod2=(float)(1.0f/(dk*dk*kmod2));

	Cnlpot_local[index]=prefac_vel*Cdens_sm_local[index]*i_kmod2; //phi_v=f*H*delta/k^2
      }
    }
  }

#ifdef _DEBUG
  printf("Node %d Transforming back\n",NodeThis);
#endif //_DEBUG
  fftwf_execute(plan_p_tor);
  fftwf_destroy_plan(plan_p_tor);
}

static void get_linearized_vpot_fd(gad_header head)
{
  int iy;
  fftwf_plan plan_p_tor;
  float agrid=Lbox/Ngrid;
  fcomplex w=(fcomplex)(cexp(I*2*M_PI/Ngrid));

  double a=head.time;
  double hub=sqrt(head.Omega0/(a*a*a)+head.OmegaLambda+(1-head.Omega0-head.OmegaLambda)/(a*a));
  double omega_m=head.Omega0/(head.Omega0+head.OmegaLambda*a*a*a+
			      (1-head.Omega0-head.OmegaLambda)*a);
  float prefac_vel=(float)(100*pow(omega_m,0.55)*hub);

#ifdef _DEBUG
  printf("Node %d Planning\n",NodeThis);
#endif //_DEBUG
  for(iy=0;iy<3;iy++)
    plan_p_tor=fftwf_mpi_plan_dft_c2r_3d(Ngrid,Ngrid,Ngrid,Cnlpot_local,Nlpot_local,
					 MPI_COMM_WORLD,FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);

#ifdef _DEBUG
  printf("Node %d computing linearized velocity\n",NodeThis);
#endif //_DEBUG
  for(iy=0;iy<Ny_here;iy++) {
    int ix;
    int iy_true=iy+Iy0_here;
    fcomplex wy=(fcomplex)(cexp(I*2*M_PI*iy_true/Ngrid));
    fcomplex wx=1.0f;
    for(ix=0;ix<Ngrid;ix++) {
      int iz;
      fcomplex wz=1.0f;
      for(iz=0;iz<Ngrid/2+1;iz++) {
	long index=iz+(Ngrid/2+1)*((long)(ix+Ngrid*iy));
	if(index==0) {
	  Clvel_local[0][index]=0;
	  Clvel_local[1][index]=0;
	  Clvel_local[2][index]=0;
	}
	else {
	  fcomplex i_kmod2=1.0f/(wx+1.0f/wx+wy+1.0f/wy+wz+1.0f/wz-6);
	  Cnlpot_local[index]=prefac_vel*Cdens_sm_local[index]*agrid*agrid*i_kmod2;
	}
	wz*=w;
      }
      wx*=w;
    }
  }

#ifdef _DEBUG
  printf("Node %d Transforming back\n",NodeThis);
#endif //_DEBUG
  fftwf_execute(plan_p_tor);
  fftwf_destroy_plan(plan_p_tor);
}

static void get_linearized_vpot(gad_header head)
{
  if(UseFD)
    get_linearized_vpot_ks(head);
  else
    get_linearized_vpot_ks(head);
}

void get_nonlinear_velocity(gad_header head)
{
  int ix;
  MPI_Status stat;
  int mod2_here=0;
  long slice_size=2*(Ngrid/2+1)*Ngrid;
  long ncells_total=Ngrid*((long)(Ngrid*Ngrid));
  double a=head.time;
  double hub=sqrt(head.Omega0/(a*a*a)+head.OmegaLambda+(1-head.Omega0-head.OmegaLambda)/(a*a));
  double omega_m=head.Omega0/(head.Omega0+head.OmegaLambda*a*a*a+
			      (1-head.Omega0-head.OmegaLambda)*a);
  float prefac_vel=(float)(100*pow(omega_m,0.55)*hub);
  double res_total=1000;
  int ngz=2*(Ngrid/2+1);
  float dx=Lbox/Ngrid;
  float i_dx=Ngrid/Lbox;
  double std_here=0,std_total=0;

  //Obtain linear solution
  get_linearized_vpot(head);

  //Compute variance of linear solution
  for(ix=0;ix<Nx_here;ix++) {
    int iy;
    long ix_0=ix=ngz*Ngrid;
    for(iy=0;iy<Ngrid;iy++) {
      int iz;
      long iy_0=iy*Ngrid;
      for(iz=0;iz<Ngrid;iz++) {
	long iz_0=iz;
	float phi=Nlpot_local[ix_0+iy_0+iz_0];
	std_here+=phi*phi;
      }
    }
  }
  MPI_Allreduce(&std_here,&std_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  std_total=sqrt(std_total/ncells_total);
  printf("STD %lE\n",std_total);

  //Transform back the smoothed density field
  get_smoothed_density_real();

  //Get slices to left and right of the density field to compute gradient
  MPI_Sendrecv(&(Dens_sm_local[(Nx_here-1)*slice_size]),slice_size,MPI_FLOAT,NodeRight,3,
	       SliceLeft_Dens,slice_size,MPI_FLOAT,NodeLeft,3,MPI_COMM_WORLD,&stat);
  MPI_Sendrecv(Dens_sm_local,slice_size,MPI_FLOAT,NodeLeft,4,
	       SliceRight_Dens,slice_size,MPI_FLOAT,NodeRight,4,MPI_COMM_WORLD,&stat);

#define MINDENS 1E-20
  //Solve non-linear equation using relaxation method
  int iter=0;
  double res_here=0;
  while(res_total>0.01) {
    //Get slices to left and right of the new potential to compute gradient
    MPI_Sendrecv(&(Nlpot_local[(Nx_here-1)*slice_size]),slice_size,MPI_FLOAT,NodeRight,1,
		 SliceLeft_Nlpot,slice_size,MPI_FLOAT,NodeLeft,1,MPI_COMM_WORLD,&stat);
    MPI_Sendrecv(Nlpot_local,slice_size,MPI_FLOAT,NodeLeft,2,
		 SliceRight_Nlpot,slice_size,MPI_FLOAT,NodeRight,2,MPI_COMM_WORLD,&stat);

    for(ix=0;ix<Nx_here;ix++) {
      int iy;
      long ix_0=ix;
      long ix_hi=ix+1;
      long ix_lo=ix-1;
      int ix_true=ix+Ix0_here;
      if(ix==0) ix_lo=Ngrid-1;
      if(ix==Ngrid-1) ix_hi=0;
      ix_0*=ngz*Ngrid;
      ix_hi*=ngz*Ngrid;
      ix_lo*=ngz*Ngrid;
      for(iy=0;iy<Ngrid;iy++) {
	int iz;
	long iy_0=iy;
	long iy_hi=iy+1;
	long iy_lo=iy-1;
	if(iy==0) iy_lo=Ngrid-1;
	if(iy==Ngrid-1) iy_hi=0;
	iy_0*=ngz;
	iy_hi*=ngz;
	iy_lo*=ngz;
	for(iz=0;iz<Ngrid;iz++) {
	  if((iz+iz+ix_true)%2==mod2_here) {
	    int ax;
	    float source,dens,inv_dens;
	    float ddphi[3],dphi[3],ddens[3];
	    long iz_0=iz;
	    long iz_hi=iz+1;
	    long iz_lo=iz-1;
	    if(iz==0) iz_lo=Ngrid-1;
	    if(iz==Ngrid-1) iz_hi=0;
	    dens=Dens_sm_local[ix_0+iy_0+iz_0];
	    if(1+dens<=MINDENS)
	      inv_dens=1./MINDENS;
	    else 
	      inv_dens=1./(1+dens);
	    if(ix==0) {
	      dphi[0]=0.5*(Nlpot_local[ix_hi+iy_0+iz_0]-SliceLeft_Nlpot[iy_0+iz_0]);
	      ddens[0]=0.5*(Dens_sm_local[ix_hi+iy_0+iz_0]-SliceLeft_Dens[iy_0+iz_0]);
	      ddphi[0]=Nlpot_local[ix_hi+iy_0+iz_0]+SliceLeft_Nlpot[iy_0+iz_0];
	    }
	    else if(ix==Nx_here-1) {
	      dphi[0]=0.5*(SliceRight_Nlpot[iy_0+iz_0]-Nlpot_local[ix_lo+iy_0+iz_0]);
	      ddens[0]=0.5*(SliceRight_Dens[iy_0+iz_0]-Dens_sm_local[ix_lo+iy_0+iz_0]);
	      ddphi[0]=SliceRight_Nlpot[iy_0+iz_0]+Nlpot_local[ix_lo+iy_0+iz_0];
	    }
	    else {
	      dphi[0]=0.5*(Nlpot_local[ix_hi+iy_0+iz_0]-Nlpot_local[ix_lo+iy_0+iz_0]);
	      ddens[0]=0.5*(Dens_sm_local[ix_hi+iy_0+iz_0]-Dens_sm_local[ix_lo+iy_0+iz_0]);
	      ddphi[0]=Nlpot_local[ix_hi+iy_0+iz_0]+Nlpot_local[ix_lo+iy_0+iz_0];
	    }
	    ddens[1]=0.5*(Dens_sm_local[ix_0+iy_hi+iz_0]-Dens_sm_local[ix_0+iy_lo+iz_0]);
	    ddens[2]=0.5*(Dens_sm_local[ix_0+iy_0+iz_hi]-Dens_sm_local[ix_0+iy_0+iz_lo]);
	    dphi[1]=0.5*(Nlpot_local[ix_0+iy_hi+iz_0]-Nlpot_local[ix_0+iy_lo+iz_0]);
	    dphi[2]=0.5*(Nlpot_local[ix_0+iy_0+iz_hi]-Nlpot_local[ix_0+iy_0+iz_lo]);
	    ddphi[1]=Nlpot_local[ix_0+iy_hi+iz_0]+Nlpot_local[ix_0+iy_lo+iz_0];
	    ddphi[2]=Nlpot_local[ix_0+iy_0+iz_hi]+Nlpot_local[ix_0+iy_0+iz_lo];

	    source=prefac_vel*dens*dens*dx*dx; // H*f*delta*h^2
	    for(ax=0;ax<3;ax++)
	      source+=dphi[ax]*ddens[ax]; // nabla(phi)*nabla(delta)
	    source*=inv_dens; //(H*f*delta*h^2+nabla(phi)*nabla(delta))/(1+delta)
	    for(ax=0;ax<3;ax++)
	      source+=ddphi[ax]; //Add rest of the laplacian
	    source/=6.;
	    res_here+=fabs(Nlpot_local[ix_0+iy_0+iz_0]-source);
	    Nlpot_local[ix_0+iy_0+iz_0]=source;
	  }
	}
      }
    }

    mod2_here=(mod2_here+1)%2;
    if(mod2_here==0) {
      iter++;
      res_total=0;
      MPI_Allreduce(&res_here,&res_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      res_total/=(ncells_total*std_total);
      res_here=0;
      if(NodeThis==0)
	printf(" %d-th iteration, residual = %lE\n",iter,res_total);
    }
  }

  //Take gradient of velocity potential (using finite differences here)
  MPI_Sendrecv(&(Nlpot_local[(Nx_here-1)*slice_size]),slice_size,MPI_FLOAT,NodeRight,1,
	       SliceLeft_Nlpot,slice_size,MPI_FLOAT,NodeLeft,1,MPI_COMM_WORLD,&stat);
  MPI_Sendrecv(Nlpot_local,slice_size,MPI_FLOAT,NodeLeft,2,
	       SliceRight_Nlpot,slice_size,MPI_FLOAT,NodeRight,2,MPI_COMM_WORLD,&stat);
  for(ix=0;ix<Nx_here;ix++) {
    int iy;
    long ix_0=ix;
    long ix_hi=ix+1;
    long ix_lo=ix-1;
    int ix_true=ix+Ix0_here;
    if(ix==0) ix_lo=Ngrid-1;
    if(ix==Ngrid-1) ix_hi=0;
    ix_0*=ngz*Ngrid;
    ix_hi*=ngz*Ngrid;
    ix_lo*=ngz*Ngrid;
    for(iy=0;iy<Ngrid;iy++) {
      int iz;
      long iy_0=iy;
      long iy_hi=iy+1;
      long iy_lo=iy-1;
      if(iy==0) iy_lo=Ngrid-1;
      if(iy==Ngrid-1) iy_hi=0;
      iy_0*=ngz;
      iy_hi*=ngz;
      iy_lo*=ngz;
      for(iz=0;iz<Ngrid;iz++) {
	long iz_0=iz;
	long iz_hi=iz+1;
	long iz_lo=iz-1;
	if(iz==0) iz_lo=Ngrid-1;
	if(iz==Ngrid-1) iz_hi=0;
	if(ix==0)
	  Nlvel_local[0][ix_0+iy_0+iz_0]=0.5*i_dx*(Nlpot_local[ix_hi+iy_0+iz_0]-SliceLeft_Nlpot[iy_0+iz_0]);
	else if(ix==Nx_here-1)
	  Nlvel_local[0][ix_0+iy_0+iz_0]=0.5*i_dx*(SliceRight_Nlpot[iy_0+iz_0]-Nlpot_local[ix_lo+iy_0+iz_0]);
	else
	  Nlvel_local[0][ix_0+iy_0+iz_0]=0.5*i_dx*(Nlpot_local[ix_hi+iy_0+iz_0]-Nlpot_local[ix_lo+iy_0+iz_0]);
	Nlvel_local[1][ix_0+iy_0+iz_0]=0.5*i_dx*(Nlpot_local[ix_0+iy_hi+iz_0]-Nlpot_local[ix_0+iy_lo+iz_0]);
	Nlvel_local[2][ix_0+iy_0+iz_0]=0.5*i_dx*(Nlpot_local[ix_0+iy_0+iz_hi]-Nlpot_local[ix_0+iy_0+iz_lo]);
      }
    }
  }
}
