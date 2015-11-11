void   make_shadow_photon(phtypeptr,phtypeptr,long,double, double);
void   advance_shadow_photon(phtypeptr,long,ctypeptr,int, double);


/********* function make shadow photon ******************************************************************/

void make_shadow_photon(sp,p,i,theta,phi)
     
     struct photon_packet_type sp[];
     struct photon_packet_type p[];  
     long i; 
     double theta,phi;
{
  
  
  sp[i].kx=sin(theta)*cos(phi);
  sp[i].ky=sin(theta)*sin(phi);
  sp[i].kz=cos(theta);

  sp[i].tau=0;
  sp[i].lum=p[i].lum;
  sp[i].weight=0;
  sp[i].starlum=p[i].starlum;
  
  sp[i].abs=p[i].abs;
  sp[i].scat=p[i].scat; 
  if (sp[i].abs!=0) sp[i].id=2;
  else if (sp[i].scat !=0) sp[i].id=1;

  sp[i].f=p[i].f;
  sp[i].x=p[i].x;
  sp[i].y=p[i].y;
  sp[i].z=p[i].z;
  sp[i].k_abs=p[i].k_abs;
  sp[i].s_scat=p[i].s_scat;
  sp[i].albedo=p[i].albedo;
  sp[i].cell=p[i].cell;
  sp[i].lcell=p[i].lcell;
  sp[i].lr=p[i].lr;
  sp[i].lx=p[i].lx;
  sp[i].ly=p[i].ly;
  sp[i].lz=p[i].lz;
  sp[i].multiplicity=1;
  sp[i].b=p[i].b;
  sp[i].r=p[i].r;
  sp[i].theta=p[i].theta;
  sp[i].phi=p[i].phi;
  sp[i].flag=p[i].flag;

  return;
}




from main--------------------


#ifdef INC_SHADOW_PHOTONS 
	     
	     make_shadow_photon(sp,p,i,theta_observer,phi_observer);
	     /* printf("kx ky kz out %.4e %.4e %.4e\n",sp[i].kx,sp[i].ky,sp[i].kz); */
	     advance_shadow_photon(sp,i,c,ncells,mfp_free);
	     if (sp[i].tau>=SHADOW_TAU_MAX-1.) sp[i].weight=0.; /* if tau very small no shadow photon */
	     else
	       {
		 sp[i].weight=1*(domega/(4.*M_PI))*exp(-sp[i].tau); /* NEED TO ACCOUNT FOR ASYMMETRIC SCATTETING!!!!! */    
		 p[i].weight=1 p[i].weight*(1-exp(-sp[i].tau));/* only for SED,iso purporses ? */ 
		 
#ifdef SPHERICAL /* count all photons that escape for the isophotal maps */
		 
		 sp[i].lr=sqrt(sp[i].lx*sp[i].lx+sp[i].ly*sp[i].ly+sp[i].lz*sp[i].lz);
		 sp[i].b=sqrt(sp[i].lr*sp[i].lr- pow(sp[i].lx*sp[i].kx+sp[i].ly*sp[i].ky+sp[i].lz*sp[i].kz,2.));
		 place_into_fbins(sp,i,&fbin,&total_fbin,0,0,TOTAL_LUM); /* add lum AND packets */
		 place_into_bbins(sp,i,&bbin,&total_bbin[0],&fbin_obs[0],f,0,0,TOTAL_LUM);      
		 
#endif
		 
		 
#ifndef SPHERICAL /* just count photons that escape in the observers direction (all shadow photons escape 
		     at observers direction */
		 for (l=0;l<NLBINS;l++)
		   for (m=0;m<NMBINS;m++)
		     { 	
     
		       place_into_fbins(sp,i,&fbin,&total_fbin,l,m,TOTAL_LUM);    
#ifdef DO_ISOMAP
		       place_into_isogrid(sp,i,&isogrid,theta_observer,phi_observer,&fbin_obs[0],l,m,dx,dy);
		       printf("hello\n");exit(0);
#endif
		     }
#endif
	       }	     
#endif
