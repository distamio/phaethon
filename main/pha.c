// PHAETHON MAIN PROGRAMME

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "grid.functions.h"
#include "inout.h"
#include "opacity.h"
#include "phaethon.functions.h"
#include "photons.h"
#include "spec.functions.h"
#include "functions.h"
#include "interactions.h"
#include "diagnostics.h"
                                /* Barnes-Hut tree routines */
#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#define   global                /* don't default to extern  */
#include "treecode.h"
#include "treedefs.h"
#include "phaethon.h"
#include "string.h"

#ifdef SPH_TREE
double mass=0;
#endif
int    nstar=-1;
double r_boundary;

int main(int argc, string argv[])
  
{
  int tid,nthreads;
  unsigned short newstar=0;
  int i,k,id;                 /* counters                   */
  int l,m; 
  struct xtable_type xtable; /* has the number of photons emitted at a specific frequency */
  char infilename[30];
#ifdef TIMING
    double rt0,rt1,rt2,rt3,rt4,rt5,rt6,tlock1=0,tlock2=0,tlock3=0,tlock4=0,tlock5=0,tlock6=0;
    FILE *fp;
#endif

  ctypeptr  c=(ctypeptr)malloc(sizeof(ctype)*MAX_CELLS);
  phtypeptr p=(phtypeptr)malloc(sizeof(phtype)*MAX_PHOTON_PACKETS);
  startypeptr star=(startypeptr)malloc(sizeof(sttype)*10);

#ifdef INC_SHADOW_PHOTONS 
  phtypeptr sp=(phtypeptr)malloc(sizeof(phtype)*10);  
#endif


#ifndef NO_OPA_TABLE

  read_opacity(&opacity);

#endif

  /* read mean opacity (in cgs units) */

  read_mean_opacity(&mean_op);


  /* table that has the frequencies at which the photons will be reemited */

  read_ftable(&reemit);

  /* read the table that has the reemision frequency probabilities */

  read_ptable(&reemit);
  /*  read_rtable(&reemit);*/

  /* create a table that has the opacities for each fpoint */

  create_opacity_table(&opa,&reemit,opacity);


 
#ifdef SPH_TREE 
#ifdef DUST_MIX
  printf("DUST MIX FOR SPH NOT READY!!\n");
  exit(0);
#endif
#endif

#ifdef DUST_MIX

  printf("\n Doing DUST MIX ...\n");
printf("\n WARNING: PARAM.OPA.H must be the same for both opacities!\n\n");

#ifndef NO_OPA_TABLE

  read_opacity2(&opacity2);

#endif

  /* read mean opacity (in cgs units) */

  read_mean_opacity2(&mean_op2);

  /* table that has the frequencies at which the photons will be reemited */

  read_ftable2(&reemit2);

  /* read the table that has the reemision frequency probabilities */

  read_ptable2(&reemit2);

  /* create a table that has the opacities for each fpoint */

  create_opacity_table(&opa2,&reemit2,opacity2);

#endif

#ifndef SPH_TREE

 // Making grid cells (creates incells.dat file) with the info of the system

    read_grid_params (opacity);
    mgrid();
#endif 
  
  // read run parameters

    read_params ();
 
 for (i=0;i<FS;i++)
     {
         switch (i)
           {
           case 0:lambda_obs[0]=F0;break;
           case 1:lambda_obs[1]=F1;break;
           case 2:lambda_obs[2]=F2;break; 
           case 3:lambda_obs[3]=F3;break;
           case 4:lambda_obs[4]=F4;break;
           case 5:lambda_obs[5]=F5;break;
           case 6:lambda_obs[6]=F6;break;
           case 7:lambda_obs[7]=F7;break;
           case 8:lambda_obs[8]=F8;break;
           case 9:lambda_obs[9]=F9;break;
           }

       }

 
 
  if (DTHETA>0.1*M_PI) printf("\nWARNING :: PHOTON COLLECTING ANGLE RELATIVELY LARGE");
  if (DTHETA>0.5*M_PI) printf("\nWARNING :: PHOTON COLLECTING ANGLE TOO LARGE! MAY CAUSE PROBLEMS!");

  if (EXTRA_SYMMETRY==1) printf("NOTE :: EXTRA SYMMETRY!\n");


#ifdef SPH_TREE

rw_dragon_file(argv[1],star,&mass,&nstar,&dragon_time, &dragon_step); 
sprintf(infilename,"phaethon.in.%s.dat",argv[1]); 
argv[1]=infilename;
make_cells(argc, argv, c, &ncells); 
record_potential_cell_info(c,ncells,mass);


if (RESTART==1)
  {  
    read_restart_info(&s_start,&output_counter_start,&k_fpoints_start,&photon_i_start,&xnum,&interactions,&rtime0);
    read_active_cells(c,global_index);
    printf("\n---------- RESTARTING RUN  -------------------------------\n");
    printf("\n RESTART       :: reading previous outcell.dat (sph) ... Done\n");
    printf(" RESTART POINT :: lum source: %d  photon: %.0f \n\n",s_start,(double)(photon_i_start));
    record_active_cells(c,ncells);
    record_active_cell_info(c,ncells,mass); 
  }

if (MORE==1)
  { 
    read_active_cells(c,global_index);
    printf("\n MORE :: reading previous outcell.dat (sph) ... Done\n\n");
    record_active_cells(c,ncells);
    record_active_cell_info(c,ncells,mass); 
  }

#endif


#ifdef ORIGIN
 printf("\nWARNING: ORIGIN defined. find_photon_origin.c in grid not used\n");
#endif	


r_boundary=R_MAX;


#ifndef SPH_TREE

 nstar=-1;   // total number of radiation sources (including ISRF) -1

  printf("\n:: Including radiation sources...\n");
 if (ISRF==1) 
	{
        printf(":: Including ISRF..\n");
	nstar=nstar+1;
        star[nstar].type=1;
	}

 if (MANUAL_STAR==1) 
	{
	printf(":: Including manual star...\n");
	nstar=nstar+1;
 	star[nstar].type=0;
	}

 if (nstar<0) 
	{
	printf(":: No radiative sources defined. Bye.\n");
	exit(0);
	}

if (RESTART==1)
  {  
    read_restart_info(&s_start,&output_counter_start,&k_fpoints_start,&photon_i_start,&xnum,&interactions,&rtime0);
    read_outcells(c, &ncells);
    read_cells_xyz(c, &ncells);
    printf("\n RESTART       :: reading previous outcell.dat  ... Done\n");
    printf(" RESTART POINT :: lum source: %d  photon: %.4e \n\n",s_start,(double)(photon_i_start));
    record_cells(c,ncells);  
  }
 else
   {
     if (MORE==0) {printf("\n::Reading incells.dat\n");read_cells(c,&ncells);}      
     else 
       { 
	 printf("\n::Reading outcells.dat\n"); 
	 read_outcells(c, &ncells);
         read_cells_xyz(c, &ncells);		
	 printf("\nAdding another lum source: reading previous outcell.dat ... Done\n\n");
	 record_cells(c,ncells); 
       }
   }
#endif  
  
  /* calculate photon step in empty space */
  
  mfp_free=find_mfp_free(c,ncells); 
  
#ifndef SPHERICAL
  create_dbins(&dbin);
#endif 
      
  create_fbins(&fbin,opacity);
  
  /* find out at which fbins are the given observation wavelenths */
  
  find_obs_fbins(&fbin_obs[0],fbin,&dlambda_obs[0]);
  
  write_lambda_opacity(fbin,opacity);  
  
#ifdef SPHERICAL
  create_bbins(&bbin,r_boundary);
  for (k=0;k<FS;k++) total_bbin[k]=0.; 
  printf("\n::SPHERICAL CASE\n");
#endif
  
  
#ifdef DO_ISOMAP
  make_isogrid(&isogrid);
#endif
 
 /*create file with time steps recorded :: BEGINNING OF TIME */
  
  make_step_info(&time0,&time1);

/*********** BEGIN STAR LOOP **************************************************/

// for msx-3D
// nstar=0; 

// star[0].type=1;
//


  for (s=s_start;s<=nstar;s=s+1)
    {            
      printf("\n--------------------   star: %d  star.type: %d -------------------------\n",s,star[s].type);
      
      k_fpoints=k_fpoints_start;
      output_counter=output_counter_start;

#ifdef DO_ISOMAP
      null_isogrid_npackets(&isogrid);  
#endif
      
 /* Initialize interaction counter table  */
      
      for (i=0;i<=50;i++) intercount[i]=0.;    
       
      if (star[s].type==0 || star[s].type==3) 
	{
#ifndef AUTO
        printf("Setting star parameters\n");
	set_star_parameters(star,s);
#else
	set_star_parameters_auto(star,s,dragon_step,dragon_time,c,p);
#endif
      
	}
      else 
	{	
	  if (star[s].type==1) set_isrf_parameters(star,s); 

	  else

	    if (star[s].type==2) set_star_isrf_parameters(star,s);
	}

           
#ifdef SPH_TREE
      
      if (star[s].type==0 && star[s].sg==1)
	{
	  printf("Making star grid\n");
	  make_star_grid(star,s,c,p,0,&ncells,&global_index);
	}  
      else
	star[s].r_star_grid=0.;
      
#endif
 
     
      /* distribute frequencies (independent of the star in x-space [dimensionless]) */
      
      if (star[s].type==1)    /** InterStellar Radiation Field **/ 
	{
	  read_ISRF(&bisrf); 
	  find_xtable_isrf(&xtable,&bisrf,star,s);
	  star[s].lum=(4*M_PI*pow(star[s].radius,2.)*M_PI*bisrf.int_total);
	  printf("\n:: ISRF total luminosity   :%.5f Lsun\n",star[s].lum/SOLAR_LUM);
	  background[0]=0;  
	  
	  for (i=1;i<FS;i++)
	    
	    background[i]= isrf(lambda_obs[i],&bisrf)*1.e17; /* in MJy/sterad */
	  /* background[i]= black_isrf(lambda_obs[i],bisrf)*
	     (c_CONST/(1e-4*1e-4*lambda_obs[i]*lambda_obs[i]));*/
	}
      
      else
	{
#ifdef STAR_SED_FILE
	  if (star[s].sedfile==1)    /** star with given SED **/ 
	    {

	      read_star_sed(star,s); 

	      find_xtable_star_sed(&xtable,star,s);
	      
	      star[s].lum=(4*M_PI*pow(star[s].radius,2.)*M_PI*star[s].int_total);
	      
	      star[s].temp=pow(star[s].lum/(star[s].dilution*4*M_PI*pow(star[s].radius,2)*sigma_CONST),0.25); 
	      star[s].r_dust=((star[s].radius/2.)*pow(star[s].temp/DUST_TEMP,3));
	      if (star[s].r_dust<R_DUST) star[s].r_dust=R_DUST;
	      
	      
	      printf(":: StarSedFile:%d\n\n",star[s].sedfile);
	      printf(":: Star lum   :%.5f Lsun\n",star[s].lum/SOLAR_LUM);
	      printf(":: Star.Rdust :%.4e (AU)\n",star[s].r_dust/AU);
	      printf(":: Star.Rdust :%.4e (cm)\n",star[s].r_dust);
	      printf(":: Star.Temp  :%.4e (K)\n",star[s].temp);
	      printf(":: Star.Radius:%.4e (Rsun)\n",star[s].radius/SOLAR_RADIUS);
	    }
	  else
#endif	    
	    find_xtable(&xtable,star,s); /* stellar blackbody radiation  */
	  
	}
    

      
    if (RESTART==1) 
	{ 
	  xtable.num[k_fpoints_start]=xnum;
	  printf("\n RESTART :: Reading interactions file...\n");
	  if (photon_i_start!=0) read_interactions(&intercount[0],s_start,star[s].npackets);
	  if (photon_i_start==0) {for (i=0;i<=50;i++) intercount[i]=0.; interactions=0;}  
	  record_interactions(&intercount[0],(double)(star[s].npackets), interactions,s,star[s].dint,photon_i);
	}
    
 
      
      /***************************************************************************************/
      /************ start injecting packets into the envelope ********************************/
      /***************************************************************************************/
      
      printf("\n*** Doing Monte Carlo simulation... Star %d : Number of photons: %.2e\n\n",s,star[s].npackets);   
 

 
#ifdef _OPENMP
    
#pragma omp parallel private(i,photon_i)   
      {
#pragma omp for schedule(dynamic, 600)     
      for (photon_i=photon_i_start+1;photon_i<=(long)(star[s].npackets); photon_i++) 
#else

	{
	for (photon_i=photon_i_start+1;photon_i<=(star[s].npackets); photon_i++)	    
#endif
         {	
#ifdef _OPENMP 
	  i=omp_get_thread_num();
#else
	  i=0;
#endif
	  /* if (i>1) printf("p_i %d\n",i); */
	 	  
#ifdef INFO_OUTPUT
	  printf("doing photon %.0f out of %0.f\n",(double) (photon_i), (1.*star[s].npackets) ); 
#endif

   

	  /* set photons positions & directions */
	  if (star[s].type==0 || star[s].type==3)  
	    {
	    init_star_photons(i,p,star,s);  /* stellar radiation  */
	  	
	    }
	  else
	    {
	      if (star[s].type==1 || star[s].type==2) /* 3D case ISRF*/
		
		init_isrf_photons_3D(i,p,star,s);  /* init_photon_params_box(i,p); */
	      
	      else
		
		if (star[s].type==4)      /* 1D case ISRF */
		  
		  init_isrf_photons_1D(i,p,star,s);  
	    }
	  	
#pragma omp critical (photon_freq)
	  {  /* set photons frequency & luminosity */
#ifdef TIMING 
	    rt1=time(NULL);
#endif
	  
	    assign_photon_freq_lum(i,&k_fpoints,p,&xtable,star,s,&pk_fpoint,&pk_num); /***/
#ifdef TIMING
	    tlock1=tlock1+time(NULL)-rt1;
#endif
	  }

	  find_photon_fpoint(p,i,&reemit);

#ifdef INFO_OUTPUT
 printf("new lambda (micron): %7.1f fpoint:%5d\n",1e4*c_CONST/p[i].f,p[i].fpoint);
#endif


	   find_photon_angles(p,i,c);	  
  

	  if (star[s].type!=3)
	    {

	      assign_photon_opacity(p,i,&opa, p[i].fpoint);

#ifdef DUST_MIX
	      assign_photon_opacity2(p,i,&opa2, p[i].fpoint);	 
#endif 
	      p[i].cell=find_photon_cell(p,i,c,star,s); 
	    }

	  else

	    {
	      p[i].cell=find_photon_cell(p,i,c,star,s); 
	      
	      if (p[i].cell>=0 || p[i].cell==-11)
		
		{ 
		  assign_photon_opacity(p,i,&opa, p[i].fpoint);
		  
#ifdef DUST_MIX
		  
		  assign_photon_opacity2(p,i,&opa2, p[i].fpoint);	 
#endif 
		}
	    }
	  
#ifdef NO_INTERACTION  /*   if i don't want the photons to interact at all  */
	  p[i].cell=-1; 
#endif
	  
	  
	  p[i].lcell=p[i].cell;
	  p[i].lx= p[i].x;
	  p[i].ly= p[i].y;
	  p[i].lz= p[i].z;
	  p[i].lr=p[i].r;

	  p[i].origin=0;
	  

	  while (p[i].cell>=0 || p[i].cell==-11)       
	    {
	     
	      /* last interaction position info */
	      
	    
	      p[i].lcell=p[i].cell;
	      p[i].lx= p[i].x;
	      p[i].ly= p[i].y;
	      p[i].lz= p[i].z;
	      p[i].lr=p[i].r;
	      /*
		advances photon till it's outside the grid or till tau=0 , i.e. 
		reaches an interaction point. By the end of advance the photon
		is either outside the grid or at an interaction cell. 
		(dens not zero) 
	      */
    
	     
	      advance_photon(p,i,c,ncells,star,s);
	      	      
	      /* DESTROY PHOTONS  if (interactions>1e9) p[i].cell=-1; */
	      
	      if (p[i].cell>=0) 
		{ 

#pragma omp atomic		  		  
		  interactions++; /***/
		  

#ifdef INFO_OUTPUT
		  printf(" %.0f %d %.0f\n",(double)(photon_i),interactions, (double)(p[i].cell));        
#endif  			 

#ifdef DUST_MIX
	     
		  if (p[i].dust_kind==1) /* second kind of dust */
		    {
		      if (p[i].albedo2<(double) (rand())/RAND_MAX) 
			{		
#pragma omp critical (absorb2)
			  {
			    p[i].temp=absorb_photon2(p[i].lum,p[i].cell,c,&mean_op2); 
			  }
			  reemit_photon2(p,i,&reemit2);
			  /* re- compute photon opacities */
			  assign_photon_opacity(p,i,&opa, p[i].fpoint);
			  assign_photon_opacity2(p,i,&opa2, p[i].fpoint);		     
			}
		      else
			{
			  scatter_photon(p,i); 
#pragma omp atomic
			  c[p[i].cell].nscat++;
			}
		    }
		  else
#endif
		    {
		      if (p[i].albedo<(double) (rand())/RAND_MAX) 
			{		
#pragma omp critical (absorb)
			  {
#ifdef TIMING
			    rt2=time(NULL);
#endif
			    p[i].temp=absorb_photon(p[i].lum,p[i].cell,c,&mean_op);
#ifdef TIMING
			    tlock2=tlock2+time(NULL)-rt1;
#endif
			  } /***/
#ifdef TIMING
			  rt3=time(NULL);
#endif
 			  reemit_photon(p,i,&reemit); 
#ifdef TIMING
			  tlock3=tlock3+time(NULL)-rt3;	
#endif
			  /* re- compute photon opacities */
#ifdef TIMING	
			  rt6=time(NULL);
#endif
			  assign_photon_opacity(p, i , &opa , p[i].fpoint); 
#ifdef DUST_MIX
			  assign_photon_opacity2(p,i,&opa2, p[i].fpoint);
#endif		     
#ifdef TIMING
			  tlock6=tlock6+time(NULL)-rt6;
#endif
			}
		      else
			{ 
#ifdef TIMING
			  rt4=time(NULL);
#endif
			  scatter_photon(p,i);	
#pragma omp atomic
			  c[p[i].cell].nscat++;	
#ifdef TIMING
			  tlock4=tlock4+time(NULL)-rt4;
#endif
			}	
		    }

		  /**** SHADOW PHOTONS: BEGIN *****/
		  
		  /* Here goes */	     
		  
		  /***** SHADOW PHOTONS: END *****/
		  
		}	       
	      
	    }	  	  
	  
	  /* find particle id */
     
	  if (p[i].abs>0) p[i].id=4;	  

	  else 
	    
	    if (p[i].scat>0) p[i].id=3;

	
	  
	  find_photon_origin(p,i,star,s);
	  
	  /* update interaction table */
#pragma omp critical (output)
	  {	
#ifdef TIMING
	    rt5=time(NULL);   
#endif
	    update_interaction_table(p,i,&intercount[0],star[s].tint,star[s].dint);/***/

	    if (p[i].origin==1) source_lum=source_lum+p[i].lum;/***/	    
	    
#ifdef SPHERICAL /* count all photons that escape for the isophotal maps */
	    
	    /* calculate impact factor */
	    p[i].b=sqrt(p[i].lr*p[i].lr- pow(p[i].lx*p[i].kx+p[i].ly*p[i].ky+p[i].lz*p[i].kz,2)); 	    
	    place_into_fbins(p,i,&fbin,&total_fbin,0.,0.,star,s); /***/
	    place_into_bbins(p,i,&bbin,&total_bbin[0],&fbin_obs[0],r_boundary);      /***/	    
#else
	    /* just count photons that escape in the observers direction */
	    /* no symmetry at all */	    	    

	    if (DIM==3 && ISODIM==3)  /***/ 
	      
	      find_fbin_isogrid_3D(p,i,dbin,&isogrid,&fbin,&fbin_obs[0], &total_fbin,star,s);
	    
	    /* some symmetry: phi symmetry and theta, -theta symmetry */
	    else
	      {
		if (DIM==2 || ISODIM==2 ) 
		  
		  find_fbin_isogrid_2D_I (p,i,dbin,&isogrid,&fbin,&fbin_obs[0], &total_fbin,star,s);
		
		/* just phi symmetry  */
		
		else
		  {
		    if (DIM==4 || ISODIM==4) 
		      
		      find_fbin_isogrid_2D_II(p,i,dbin,&isogrid,&fbin,&fbin_obs[0], &total_fbin,star,s);
		  }
	      }
	    
#endif
	    /* printf ("%7.0f %7d %7d %d \n", (double)(photon_i), k_fpoints, p[i].fpoint, i);*/ /***/ 
	    

	    /* record_photons_info(p,photon_i);      */
	    
	    if (photon_i==star[s].output+output_counter && photon_i!=star[s].npackets)/* conditional output */ /***/
	      {  
		
		output_counter=star[s].output+output_counter;
		printf("star: %d --> DONE photon %.2e out of %.2e\n", s, (double)(photon_i),(double)(star[s].npackets));
		update_step_info((double)(photon_i),(double)(star[s].npackets), time0,&time1,s,rtime0);
                
		write_restart_info(s,output_counter,pk_fpoint,(double) (photon_i),pk_num,interactions,time0,rtime0);
            
#ifdef SPH_TREE
		record_active_cells(c,ncells);
		record_active_cell_info(c,ncells,mass); 
#else	      
		record_cells(c,ncells); 	 
#endif	
             
		record_interactions(&intercount[0],(double)(star[s].npackets), interactions,s,star[s].dint,photon_i);

                write_sed(&fbin,&bisrf,dbin,r_boundary,&sed_flag); 

		init_fbins(&fbin,opacity); /* needs to be zero  when input from previous file */		

#ifdef SPHERICAL                
                
#ifdef DO_ISOMAP                
                find_bbin_intensity(&bbin,&lambda_obs[0],&dlambda_obs[0],star,s,&iso_flag);      
                write_isophotal_map2(&bbin,"isodata"); 
                make_isogrid_spherical_case(&bbin,&isogrid,&background[0]);  /* ... transform to pixels */
                isomap_idl(&isogrid);
                null_isogrid_npackets(&isogrid);  	       
                null_bbins_npackets(&bbin);
#endif 

#else 
		
#ifdef DO_ISOMAP  
	       find_isogrid_intensity(&isogrid,&lambda_obs[0],&background[0],dbin,&dlambda_obs[0],star,s,r_boundary,&iso_flag);
	       isomap_idl(&isogrid);	       
#ifdef STATISTICS 
	       isomap_statistics_idl(&isogrid,s);
#endif 
            null_isogrid_npackets(&isogrid);  	       

#endif 
#endif  
	      }
#ifdef TIMING
  tlock5=tlock5+time(NULL)-rt5;
#endif
	  }
	}
      }

      
      photon_i=star[s].npackets;
      
//***************************************************************************************/
//*********** END OF PROPAGATING PHOTONS LOOP  **********************************/
//***************************************************************************************/
      

      
      /* 2D correct temperature improving statistics for some symmetry */
      
      if (DIM==2)    
	{
	  fix_temp(c); 
	}
      
      /* write final outcell.dat file (cells+info) */ 
      
     
      printf("star: %d --> DONE photon %.2e out of %.2e\n", s,(double)(photon_i) ,(double)(star[s].npackets));
      
      write_restart_info(s+1,0,0,0,0,0,time0,rtime0);

#ifdef SPH_TREE
      record_active_cells(c,ncells);
      record_active_cell_info(c,ncells,mass);
#else
      record_cells(c,ncells); 
#endif
      
      /* step info file */
      final_update_step_info((double)(star[s].npackets),(double)(star[s].npackets), time0,&time1,s,rtime0);
      
      /* record interaction file */
      record_interactions(&intercount[0],star[s].npackets, interactions,s,star[s].dint,photon_i);
      
      printf("\n- Percent of %d stellar photons f-binned (within specified lamba-range): %.2f %s\n",s, 100.*total_fbin/star[s].npackets,"%");
      printf("\n- Total %d star photons interactions       : %.0f \n",s, interactions);
    

    
#ifdef SPHERICAL
 #ifdef DO_ISOMAP
 find_bbin_intensity(&bbin,&lambda_obs[0],&dlambda_obs[0],star,s,&iso_flag);      
 write_isophotal_map2(&bbin,"isodata"); 
 make_isogrid_spherical_case(&bbin,&isogrid,&background[0]);  /* ... transform to pixels */
 isomap_idl(&isogrid);
#endif 

#else
       
#ifdef DO_ISOMAP 
      find_isogrid_intensity(&isogrid,&lambda_obs[0],&background[0],dbin,&dlambda_obs[0],star,s,r_boundary,&iso_flag);
      isomap_idl(&isogrid);
#endif  
      
#ifdef STATISTICS 
      isomap_statistics_idl(&isogrid,s);
#endif

#endif

      interactions=0;
      total_fbin=0; 
      
       write_sed(&fbin,&bisrf,dbin,r_boundary,&sed_flag);
 
    }
    
  /******* END STAR LOOP ********************************************************************************/
  





 time2=time(NULL); 
 
 printf("\n- TOTAL CPU TIME: %.2f min \n",(time2-time0)/60);


#ifdef TIMING
 fp=fopen("timing","w");

 fprintf(fp, "\n-- TIME--pc-\nfreq   :%5.2f \nabsorb :%5.2f \nreemit :%5.2f \nscatter:%5.2f \nupdate :%5.2f\nassign :%5.2f\n",
	100*tlock1/(time2-time0),100*tlock2/(time2-time0),100*tlock3/(time2-time0),
	 100*tlock4/(time2-time0),100*tlock5/(time2-time0),100*tlock6/(time2-time0) );

 fprintf(fp, "\n-- TIME-----\nfreq   :%5.2f \nabsorb :%5.2f \nreemit :%5.2f \nscatter:%5.2f \nupdate :%5.2f\nassign :%5.2f\n",
	tlock1 ,tlock2 ,tlock3 ,tlock4 ,tlock5 ,tlock6);
 fclose(fp);
#endif


 return (0);

}

