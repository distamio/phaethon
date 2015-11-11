local void mynewtree(void);                     /* flush existing tree      */
local void mynewtree(void)
{
    static bool firstcall = TRUE;
    nodeptr p;
    int i;
    if (! firstcall) {                          /* if cells to reclaim      */
        p = (nodeptr) root;                     /* start with the root      */
        while (p != NULL)                       /* loop scanning tree       */
            if (Type(p) == CELL) {              /* if we found a cell to    */
                Next(p) = freecell;             /* then save existing list  */
                freecell = p;                   /* and add it to the front  */
                p = More(p);                    /* then scan down tree      */
            } else                              /* else, skip over bodies   */
                p = Next(p);                    /* by going on to the next  */
    } else                                      /* else nothing to reclaim  */
        firstcall = FALSE;                      /* so just note it          */
    root = NULL;                                /* flush existing tree      */
    ncell = 0;                                  /* reset cell count         */
    real_ncell=1;                               /* reset real cell count         */
    for (i = 0; i < NSUB-1; i++) nextcell[i]=i+1; /* set next list for root children */
}

/*
 * MY-MAKETREE: initialize tree structure for hierarchical force calculation.
 */

void mymaketree(bodyptr btab, int nbody)
{
    double cpustart;
    bodyptr p;
    int i;

    cpustart = cputime();                       /* record time at start     */
    mynewtree();                                /* flush existing tree, etc */
    root = mymakecell();                        /* allocate the root cell   */
    CLRV(Pos(root));                            /* initialize the midpoint  */
    expandbox(btab, nbody);                     /* and expand cell to fit   */
    for (p = btab; p < btab+nbody; p++)         /* loop over all bodies     */
        myloadbody(p);                          /* insert each into tree    */
    bh86 = scanopt(options, "bh86");            /* set flags for alternate  */
    sw94 = scanopt(options, "sw94");            /* ...cell opening criteria */
    if (bh86 && sw94)                           /* can't have both at once  */
        error("maketree: incompatible options bh86 and sw94\n");
    tdepth = 0;                                 /* init count of levels     */
    for (i = 0; i < MAXLEVEL; i++)              /* and init tree histograms */
      cellhist[i] = subnhist[i] = 0;
    hackcofm(root, rsize, 0);                   /* find c-of-m coords, etc  */
    threadtree((nodeptr) root, NULL);           /* add next and more links  */
    if (usequad)                                /* if including quad terms  */
      hackquad(root);                           /* find quadrupole moments  */
    cputree = cputime() - cpustart;             /* store elapsed CPU time   */

}

/*
 * MY-MAKETREE: initialize tree structure for hierarchical force calculation.
 */

void mymaketree2(bodyptr btab, int nbody)
{
    double cpustart;
    bodyptr p;
    int i;

    cpustart = cputime();                       /* record time at start     */
    mynewtree();                                /* flush existing tree, etc */
    root = mymakecell();                        /* allocate the root cell   */
    CLRV(Pos(root));                            /* initialize the midpoint  */
    expandbox(btab, nbody);                     /* and expand cell to fit   */
    for (p = btab; p < btab+nbody; p++)         /* loop over all bodies     */
        myloadbody(p);                          /* insert each into tree    */
    bh86 = scanopt(options, "bh86");            /* set flags for alternate  */
    sw94 = scanopt(options, "sw94");            /* ...cell opening criteria */
    if (bh86 && sw94)                           /* can't have both at once  */
        error("maketree: incompatible options bh86 and sw94\n");
    tdepth = 0;                                 /* init count of levels     */
    for (i = 0; i < MAXLEVEL; i++)              /* and init tree histograms */
        cellhist[i] = subnhist[i] = 0;
    /* hackcofm(root, rsize, 0); */                 /* find c-of-m coords, etc  */
    /* threadtree((nodeptr) root, NULL);  */        /* add next and more links  */
    /* if (usequad)   */                             /* if including quad terms  */
    /*  hackquad(root);   */                      /* find quadrupole moments  */
    cputree = cputime() - cpustart;             /* store elapsed CPU time   */

}






void mymaketree(bodyptr, int);          /* construct tree structure       */
void mymaketree2(bodyptr, int);         /* construct tree structure   no threading *
