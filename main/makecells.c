

#include "structs.h"

                         /* Barnes-Hut tree routines */

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#define global                                  /* don't default to extern  */
#include "treecode.h"
#include "treedefs.h" 
#include "inout.h" 
#include <stdio.h>
#include <stdlib.h> 
#include "getparam.h"

local void startrun(void);

void make_cells(int argc, string argv[],struct cell_type c[], int *ncells_ptr)
     
{
  long counter=1;
  double mass_scale,size_scale,dens_scale; 
  nodeptr q,w; 

string defv[] = {       ";Hierarchical N-body code "
#if defined(QUICKSCAN)
			"(quick scan)",
#else
			"(theta scan)",
#endif
			"in=",                      ";Input file with initial conditions",
			"out=",                     ";Output file of N-body frames",
#if defined(USEFREQ)
			"freq=32.0",                ";Leapfrog integration frequency",
#else
			"dtime=1/32",               ";Leapfrog integration timestep",
#endif
			"eps=0.025",                ";Density smoothing length",
#if !defined(QUICKSCAN)
			"theta=1.0",                ";Force accuracy parameter",
#endif
			"usequad=false",            ";If true, use quad moments",
			"options=",                 ";Various control options",
			"tstop=2.0",                ";Time to stop integration",
#if defined(USEFREQ)
			"freqout=4.0",              ";Data output frequency",
#else
			"dtout=1/4",                ";Data output timestep",
#endif
			"nbody=4096",               ";Number of bodies for test run",
			"seed=123",                 ";Random number seed for test run",
			"save=",                    ";Write state file as code runs",
				"restore=",                 ";Continue run from state file",
			"VERSION=1.4",              ";Joshua Barnes  February 21 2001",
			NULL,
};




 
#ifdef TREE_INFO_OUT 
  printf("count body=1     x           y           z       c.index  c.size     nsph   'index'  mass    dens\n");
#endif
 

  global_index=0;

  printf("\n\n-------------------- TREE INFO -----------------------------------\n\n");
  initparam(argv, defv);                      /* initialize param access  */
  headline = defv[0] + 1;                     /* skip ";" in headline     */
  startrun();                                 /* get params & input data  */                 
  /* startoutput(); */                             /* activate output code     */

  printf("Number of SPH particles: %8d%", nbody);
  maketree(bodytab, nbody);                   /* construct tree structure */ 

  printf("\nroot:  x %.4e  y  %.4e  z %.4e\n", Pos((nodeptr)(root))[0], Pos(root)[1], Pos(root)[2]);

  printf("ncell parameter : %d\n", ncell);
 


#ifdef TREE_INFO_OUT 
  printf("count body=1     x           y           z       c.index  c.size     nsph   'index'  mass    dens\n");
#endif
 
/* loop over all cells (i.e grouped sph particles) !!!! */
  
  counter =1;
  
  for (q=More(root);q!=Next(root);q=Next(q)){    
    
    {size_scale=csize(q);mass_scale=Mass(q);dens_scale=cdens(q);}
    
    while (q!=NULL){                       
      
      if (cnum(q)>1) 
	{    
	  w=q;
	  
	  /* printf("###%4d (%d) %7.4e %6.4e %6.4e %6d   %6.4e   %.0f   %4d %6.4e  %6.4e \n",
	     counter,Type(q),Pos(q)[0],Pos(q)[1],Pos(q)[2],cid(q),csize(q),cnum(q),cindex(q),
	     Mass(q),Mass(q)/csize(q)); */	  
	  
	  c[cid(q)-1].x=Pos(q)[0];
	  c[cid(q)-1].y=Pos(q)[1];
	  c[cid(q)-1].z=Pos(q)[2];
	  c[cid(q)-1].mass=dens_scale*pow(csize(q),3);
	  c[cid(q)-1].size=csize(q);
	  c[cid(q)-1].dens=dens_scale;
	  c[cid(q)-1].temp=0.;
	  c[cid(q)-1].nsph=cnum(q);
	  c[cid(q)-1].nabs=0.;
	  c[cid(q)-1].nscat=0.;
	  c[cid(q)-1].abslum=0.;
	  c[cid(q)-1].newcell=0;
	  c[cid(q)-1].advsize=size_scale;
#ifdef TREE_INFO_OUT
	  
	  printf("###%4d (%d) %7.4e %6.4e %6.4e %6d   %6.4e   %7.0f   %4d %6.4e  %6.4e \n",
		 counter,Type(q),Pos(q)[0],Pos(q)[1],Pos(q)[2],
		 cid(q),c[cid(q)-1].size,cnum(q),cindex(q),c[cid(q)-1].mass,c[cid(q)-1].dens);
#endif
	  counter=counter+1;
	  
	  
	  if (More(q)==NULL) 
	    break; 
	  else 
	    {
	      if (cnum(q)>25) 
		{size_scale=csize(q);mass_scale=Mass(q);dens_scale=cdens(q);}
	      q=More(w); 
	      if (cnum(q)>25) 
		{size_scale=csize(q);mass_scale=Mass(q);dens_scale=cdens(q);}
	    }
	}
      
      else
	{    
	  	  
	  w=q;
	  /*
	    printf("###%4d (%d) %7.4e %6.4e %6.4e %6d   %6.4e   %.0f   %4d %6.4e   %6.4e \n",
	    counter,Type(q),Pos(q)[0],Pos(q)[1],Pos(q)[2],cid(q),csize(q),
	    cnum(q),cindex(q),Mass(q),Mass(q)/csize(q));
	  */ 
	  c[cid(q)-1].x=Pos(q)[0];
	  c[cid(q)-1].y=Pos(q)[1];
	  c[cid(q)-1].z=Pos(q)[2];
	  c[cid(q)-1].mass=dens_scale*pow(csize(q),3);
	  c[cid(q)-1].size=csize(q);
	  c[cid(q)-1].dens=dens_scale;
	  c[cid(q)-1].temp=0.;
	  c[cid(q)-1].nsph=cnum(q);
	  c[cid(q)-1].nabs=0.;
	  c[cid(q)-1].nscat=0.;
	  c[cid(q)-1].abslum=0.;
	  c[cid(q)-1].newcell=0;
	  c[cid(q)-1].advsize=size_scale;

#ifdef TREE_INFO_OUT
	  

	  printf("###%4d (%d) %7.4e %6.4e %6.4e %6d   %6.4e   %7.0f   %4d %6.4e  %6.4e \n",
		 counter,Type(q),Pos(q)[0],Pos(q)[1],Pos(q)[2],
		 cid(q),c[cid(q)-1].size,cnum(q),cindex(q),c[cid(q)-1].mass,c[cid(q)-1].dens);
	
#endif
	  counter=counter+1;
	  break; 
	}      
    }    
  }
  
  /* printf("\ncellsum size %.0f  rsize %0.f\n",pow(cellsum,1/3.), rsize);   */
  
  *ncells_ptr=global_index-1;
  printf("ncells (counting starts from 0) : %d\n",*ncells_ptr); 
  printf("root size %.4e\n",rsize);
  printf("global index : %d\n", global_index);

  boundary=0.5*rsize;
  write_cells(c,*ncells_ptr);   

return;

}

/*
 * STARTRUN: startup hierarchical N-body code.
 */

local void startrun(void)
{
#if !defined(USEFREQ)
    double dt1, dt2;
#endif

    infile = getparam("in");                    /* set I/O file names       */
    outfile = getparam("out");
    savefile = getparam("save");
    if (strnull(getparam("restore"))) {         /* if starting a new run    */
        eps = getdparam("eps");                 /* get input parameters     */
#if defined(USEFREQ)
        freq = getdparam("freq");
#else
        dtime = (sscanf(getparam("dtime"), "%lf/%lf", &dt1, &dt2) == 2 ?
                 dt1 / dt2 : getdparam("dtime"));
#endif
#if !defined(QUICKSCAN)
        theta = getdparam("theta");
#endif
        usequad = getbparam("usequad");
        tstop = getdparam("tstop");
#if defined(USEFREQ)
        freqout = getdparam("freqout");
#else
        dtout = (sscanf(getparam("dtout"), "%lf/%lf", &dt1, &dt2) == 2 ?
                 dt1 / dt2 : getdparam("dtout"));
#endif
        options = getparam("options");
        if (! strnull(infile))                  /* if data file was given   */
            inputdata();                        /* then read inital data    */
        else {                                  /* else make initial data   */	  
	  error("you need to input the sph file mate! \n");
	}
        rsize = 1.0;                            /* start root w/ unit cube  */
        nstep = 0;                              /* begin counting steps     */
        tout = tnow;                            /* schedule first output    */
    } else {                                    /* else restart old run     */
        restorestate(getparam("restore"));      /* read in state file       */
        if (getparamstat("eps") & ARGPARAM)     /* if given, set new params */
            eps = getdparam("eps");
#if !defined(QUICKSCAN)
        if (getparamstat("theta") & ARGPARAM)
            theta = getdparam("theta");
#endif
        if (getparamstat("usequad") & ARGPARAM)
            usequad = getbparam("usequad");
        if (getparamstat("options") & ARGPARAM)
            options = getparam("options");
        if (getparamstat("tstop") & ARGPARAM)
            tstop = getdparam("tstop");
#if defined(USEFREQ)
        if (getparamstat("freqout") & ARGPARAM)
            freqout = getdparam("freqout");
        if (scanopt(options, "new-tout"))       /* if output time reset     */
            tout = tnow + 1 / freqout;          /* then offset from now     */
#else
            dtout = (sscanf(getparam("dtout"), "%lf/%lf", &dt1, &dt2) == 2 ?
                      dt1 / dt2 : getdparam("dtout"));
        if (scanopt(options, "new-tout"))       /* if output time reset     */
            tout = tnow + dtout;                /* then offset from now     */
#endif
    }
}

