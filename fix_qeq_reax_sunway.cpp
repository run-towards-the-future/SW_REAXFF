/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

     Hybrid and sub-group capabilities: Ray Shan (Sandia)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_qeq_reax_sunway.h"
#include "pair_reaxc_sunway.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "pair.h"
#include "respa.h"
#include "memory.h"
#include "citeme.h"
#include "error.h"
#include "reaxc_defs_sunway.h"
#include "sunway.h"
#include "gptl.h"

extern"C"
{
  void sparse_matvec_C(void *param);
  void sparse_matvec_C_spawn(void *param);
  void sparse_matvec_C_join();
  void compute_H_Full_C(void *param);
}

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace REAXC_SUNWAY_NS;
#define EV_TO_KCAL_PER_MOL 14.4
//#define DANGER_ZONE     0.95
//#define LOOSE_ZONE      0.7
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MIN_NBRS 100

static const char cite_fix_qeq_reax[] =
  "fix qeq/reax command:\n\n"
  "@Article{Aktulga12,\n"
  " author = {H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama},\n"
  " title = {Parallel reactive molecular dynamics: Numerical methods and algorithmic techniques},\n"
  " journal = {Parallel Computing},\n"
  " year =    2012,\n"
  " volume =  38,\n"
  " pages =   {245--259}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixQEqReaxSunway::FixQEqReaxSunway(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), pertype_option(NULL)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_qeq_reax);

  if (narg<8 || narg>9) error->all(FLERR,"Illegal fix qeq/reax command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix qeq/reax command");

  swa = force->numeric(FLERR,arg[4]);
  swb = force->numeric(FLERR,arg[5]);
  tolerance = force->numeric(FLERR,arg[6]);
  int len = strlen(arg[7]) + 1;
  pertype_option = new char[len];
  strcpy(pertype_option,arg[7]);

  // dual CG support only available for USER-OMP variant
  // check for compatibility is in Fix::post_constructor()
  dual_enabled = 0;
  if (narg == 9) {
    if (strcmp(arg[8],"dual") == 0) dual_enabled = 1;
    else error->all(FLERR,"Illegal fix qeq/reax command");
  }
  shld = NULL;

  n = n_cap = 0;
  N = nmax = 0;
  m_fill = m_cap = 0;
  pack_flag = 0;
  s = NULL;
  t = NULL;
  nprev = 5;

  Hdia_inv = NULL;
  b_s = NULL;
  b_t = NULL;
  b_prc = NULL;
  b_prm = NULL;

  // CG
  p = NULL;
  q = NULL;
  r = NULL;
  d = NULL;

  // H matrix
  H.firstnbr = NULL;
  H.numnbrs = NULL;
  H.jlist = NULL;
  H.val = NULL;

  // dual CG support
  // Update comm sizes for this fix
  if (dual_enabled) comm_forward = comm_reverse = 2;
  else comm_forward = comm_reverse = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  reaxc = NULL;
  reaxc = (PairReaxCSunway *) force->pair_match("reax/c",0);

  if (reaxc) {
    s_hist = t_hist = NULL;
    grow_arrays(atom->nmax);
    atom->add_callback(0);
    for (int i = 0; i < atom->nmax; i++)
      for (int j = 0; j < nprev; ++j)
        s_hist[i][j] = t_hist[i][j] = 0;
  }
}

/* ---------------------------------------------------------------------- */

FixQEqReaxSunway::~FixQEqReaxSunway()
{
  if (copymode) return;

  delete[] pertype_option;

  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  memory->destroy(s_hist);
  memory->destroy(t_hist);

  deallocate_storage();
  deallocate_matrix();

  memory->destroy(shld);

  if (!reaxflag) {
    memory->destroy(chi);
    memory->destroy(eta);
    memory->destroy(gamma);
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::post_constructor()
{
  pertype_parameters(pertype_option);
  if (dual_enabled)
    error->all(FLERR,"Dual keyword only supported with fix qeq/reax/omp");
}

/* ---------------------------------------------------------------------- */

int FixQEqReaxSunway::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::pertype_parameters(char *arg)
{
  if (strcmp(arg,"reax/c") == 0) {
    reaxflag = 1;
    Pair *pair = force->pair_match("reax/c",0);
    if (pair == NULL) error->all(FLERR,"No pair reax/c for fix qeq/reax");

    int tmp;
    chi = (double *) pair->extract("chi",tmp);
    eta = (double *) pair->extract("eta",tmp);
    gamma = (double *) pair->extract("gamma",tmp);
    if (chi == NULL || eta == NULL || gamma == NULL)
      error->all(FLERR,
                 "Fix qeq/reax could not extract params from pair reax/c");
    return;
  }

  int i,itype,ntypes;
  double v1,v2,v3;
  FILE *pf;

  reaxflag = 0;
  ntypes = atom->ntypes;

  memory->create(chi,ntypes+1,"qeq/reax:chi");
  memory->create(eta,ntypes+1,"qeq/reax:eta");
  memory->create(gamma,ntypes+1,"qeq/reax:gamma");

  if (comm->me == 0) {
    if ((pf = fopen(arg,"r")) == NULL)
      error->one(FLERR,"Fix qeq/reax parameter file could not be found");

    for (i = 1; i <= ntypes && !feof(pf); i++) {
      fscanf(pf,"%d %lg %lg %lg",&itype,&v1,&v2,&v3);
      if (itype < 1 || itype > ntypes)
        error->one(FLERR,"Fix qeq/reax invalid atom type in param file");
      chi[itype] = v1;
      eta[itype] = v2;
      gamma[itype] = v3;
    }
    if (i <= ntypes) error->one(FLERR,"Invalid param file for fix qeq/reax");
    fclose(pf);
  }

  MPI_Bcast(&chi[1],ntypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&eta[1],ntypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamma[1],ntypes,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::allocate_storage()
{
  nmax = atom->nmax;

  memory->create(s,nmax,"qeq:s");
  memory->create(t,nmax,"qeq:t");
  //printf("nmax = %d, sizeof(float) = %d, sizeof(double) = %d\n", nmax, sizeof(float), sizeof(double));

  memory->create(Hdia_inv,nmax,"qeq:Hdia_inv");
  memory->create(b_s,nmax,"qeq:b_s");
  memory->create(b_t,nmax,"qeq:b_t");
  memory->create(b_prc,nmax,"qeq:b_prc");
  memory->create(b_prm,nmax,"qeq:b_prm");

  // dual CG support
  int size = nmax;
  if (dual_enabled) size*= 2;

  memory->create(p,size,"qeq:p");
  memory->create(q,size,"qeq:q");
  memory->create(r,size,"qeq:r");
  memory->create(d,size,"qeq:d");

  //memory->create(m,size,"qeq:m");//Pipeline PCG
  //memory->create(v,size,"qeq:v");
  //memory->create(u,size,"qeq:u");
  //memory->create(w,size,"qeq:w");
  //memory->create(z,size,"qeq:z");
  
  
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::deallocate_storage()
{
  memory->destroy(s);
  memory->destroy(t);

  memory->destroy( Hdia_inv );
  memory->destroy( b_s );
  memory->destroy( b_t );
  memory->destroy( b_prc );
  memory->destroy( b_prm );

  memory->destroy( p );
  memory->destroy( q );
  memory->destroy( r );
  memory->destroy( d );

  //memory->destroy( m );
  //memory->destroy( v );
  //memory->destroy( u );
  //memory->destroy( w );
  //memory->destroy( z );

  
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::reallocate_storage()
{
  deallocate_storage();
  allocate_storage();
  init_storage();
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::allocate_matrix()
{
  int i,ii,inum,m;
  int *ilist, *numneigh;

  int mincap;
  double safezone;

  if (reaxflag) {
    mincap = reaxc->system->mincap;
    safezone = reaxc->system->safezone;
  } else {
    mincap = MIN_CAP;
    safezone = SAFE_ZONE;
  }

  
  // determine the total space for the H matrix
  if (reaxc) {
    inum = reaxc->listfull->inum;
    ilist = reaxc->listfull->ilist;
    numneigh = reaxc->listfull->numneigh;
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
  }

  /* original method for memory allocation */
  //n = atom->nlocal + atom->nghost;
  //n_cap = MAX( (int)(n * safezone), mincap);
  //m = 0;
  //for (ii = 0; ii < inum; ii++) {
  //  i = ilist[ii];
  //  m += numneigh[i];
  //}
  //m_cap = MAX( (int)(m * safezone), mincap * MIN_NBRS);
  //H.n = n_cap;
  //H.m = m_cap;
  //memory->create(H.firstnbr,n_cap,"qeq:H.firstnbr");
  //memory->create(H.numnbrs,n_cap,"qeq:H.numnbrs");
  //memory->create(H.jlist,m_cap,"qeq:H.jlist");
  //memory->create(H.val,m_cap,"qeq:H.val");
 
  /* modified method for memory allocation */
  //inum = nlocal = 20736;
  maxHlist = 600;
  int natoms = maxHlist * (atom->nlocal+64);
  memory->create(H.firstnbr,natoms,"qeq:H.firstnbr");
  memory->create(H.numnbrs, natoms,"qeq:H.numnbrs");
  memory->create(H.jlist,   natoms,"qeq:H.jlist");
  memory->create(H.val,     natoms,"qeq:H.val");

  //H_newtoff.n = n_cap;
  //H_newtoff.m = m_cap;
  //memory->create(H_newtoff.firstnbr,n_cap,"qeq:H_newtoff.firstnbr");
  //memory->create(H_newtoff.numnbrs,n_cap,"qeq:H_newtoff.numnbrs");
  //memory->create(H_newtoff.jlist,m_cap,"qeq:H_newtoff.jlist");
  //memory->create(H_newtoff.val,m_cap,"qeq:H_newtoff.val");

}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::deallocate_matrix()
{
  memory->destroy( H.firstnbr );
  memory->destroy( H.numnbrs );
  memory->destroy( H.jlist );
  memory->destroy( H.val );

  //memory->destroy( H_newtoff.firstnbr );
  //memory->destroy( H_newtoff.numnbrs );
  //memory->destroy( H_newtoff.jlist );
  //memory->destroy( H_newtoff.val );

}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::reallocate_matrix()
{
  deallocate_matrix();
  allocate_matrix();
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/reax requires atom attribute q");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/reax group has no atoms");

  // need a half neighbor list w/ Newton off and ghost neighbors
  // built whenever re-neighboring occurs

  // int irequest = neighbor->request(this,instance_me);
  // neighbor->requests[irequest]->pair = 0;
  // neighbor->requests[irequest]->fix = 1;
  // neighbor->requests[irequest]->newton = 2;
  // neighbor->requests[irequest]->ghost = 1;

  init_shielding();
  init_taper();

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::init_list(int id, NeighList *ptr)
{
    list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::init_shielding()
{
  int i,j;
  int ntypes;

  ntypes = atom->ntypes;
  if (shld == NULL)
  {
    memory->create(shld,ntypes+1,ntypes+1,"qeq:shielding");
  }

  for (i = 1; i <= ntypes; ++i)
    for (j = 1; j <= ntypes; ++j)
      shld[i][j] = pow( gamma[i] * gamma[j], -1.5);
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::init_taper()
{
  double d7, swa2, swa3, swb2, swb3;

  if (fabs(swa) > 0.01 && comm->me == 0)
    error->warning(FLERR,"Fix qeq/reax has non-zero lower Taper radius cutoff");
  if (swb < 0)
    error->all(FLERR, "Fix qeq/reax has negative upper Taper radius cutoff");
  else if (swb < 5 && comm->me == 0)
    error->warning(FLERR,"Fix qeq/reax has very low Taper radius cutoff");

  d7 = pow( swb - swa, 7);
  swa2 = SQR( swa);
  swa3 = CUBE( swa);
  swb2 = SQR( swb);
  swb3 = CUBE( swb);

  Tap[7] =  20.0 / d7;
  Tap[6] = -70.0 * (swa + swb) / d7;
  Tap[5] =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
  Tap[4] = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3) / d7;
  Tap[3] = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3) / d7;
  Tap[2] =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
  Tap[1] = 140.0 * swa3 * swb3 / d7;
  Tap[0] = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 +
            7.0*swa*swb3*swb3 + swb3*swb3*swb) / d7;
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::setup_pre_force(int vflag)
{
  deallocate_storage();
  allocate_storage();

  init_storage();

  deallocate_matrix();
  allocate_matrix();

  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::setup_pre_force_respa(int vflag, int ilevel)
{
  if (ilevel < nlevels_respa-1) return;
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::min_setup_pre_force(int vflag)
{
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::init_storage()
{
  int NN;

  if (reaxc)
    NN = reaxc->listfull->inum + reaxc->listfull->gnum;
  else
    NN = list->inum + list->gnum;

  for (int i = 0; i < NN; i++) {
    Hdia_inv[i] = 1. / eta[atom->type[i]];
    b_s[i] = -chi[atom->type[i]];
    b_t[i] = -1.0;
    b_prc[i] = 0;
    b_prm[i] = 0;
    s[i] = t[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::pre_force(int vflag)
{
  GPTLstart("fix qeq reax pre force");
  double t_start, t_end;

  if (update->ntimestep % nevery) return;
  if (comm->me == 0) t_start = MPI_Wtime();

  n = atom->nlocal;
  N = atom->nlocal + atom->nghost;

  // grow arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) reallocate_storage();
  if (n > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();
  GPTLstart("fix qeq reax pre force init matvec");
  init_matvec();
  GPTLstop("fix qeq reax pre force init matvec");
  
  GPTLstart("fix qeq reax pre force CG");
  
  int sz = reaxc->listfull->inum + reaxc->listfull->gnum;
  if (dual_enabled) sz*= 2;
  memory->create(m,sz,"qeq:m");//Pipeline PCG
  memory->create(v,sz,"qeq:v");
  memory->create(u,sz,"qeq:u");
  memory->create(w,sz,"qeq:w");
  memory->create(z,sz,"qeq:z");
  memory->create(ss,sz,"qeq:ss");
  memory->create(n_test,sz,"qeq:p_test");

  //matvecs_s = CG_new(b_s, s);    	// CG on s - parallel
  //matvecs_t = CG_new(b_t, t);       // CG on t - parallel
  
  matvecs_s = CG_v2(b_s, s);    	// CG on s - parallel
  matvecs_t = CG_v2(b_t, t);       // CG on t - parallel

  memory->destroy(m);
  memory->destroy(v);
  memory->destroy(u);
  memory->destroy(w);
  memory->destroy(z);
  memory->destroy(ss);
  memory->destroy(n_test);

  GPTLstop("fix qeq reax pre force CG");
  
  
  matvecs = matvecs_s + matvecs_t;

  GPTLstart("fix qeq reax pre force calculate_Q");
  calculate_Q();
  GPTLstop("fix qeq reax pre force calculate_Q");
  if (comm->me == 0) 
  {
    t_end = MPI_Wtime();
    qeq_time = t_end - t_start;
  }
  GPTLstop("fix qeq reax pre force");
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::pre_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::init_matvec()
{
  /* fill-in H matrix */
  GPTLstart("init matvec compute_H");
  int nn;//, ii;
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj, flag;
  double dx, dy, dz, r_sqr;
  const double SMALL = 0.0001;

  int *type = atom->type;
  tagint *tag = atom->tag;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (reaxc) 
  {
    nn = reaxc->listfull->inum;
    inum = reaxc->listfull->inum + reaxc->listfull->gnum;
    ilist = reaxc->listfull->ilist;
    numneigh = reaxc->listfull->numneigh;
    firstneigh = reaxc->listfull->firstneigh;
  } 
  else 
  {
    nn = list->inum;
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  }
  
  m_fill = 0;
  r_sqr = 0;
  double swbsq = swb * swb;

  fix_qeq_pack_t param;
  param.swbsq           = swbsq;
  param.shld            = shld;
  param.groupbit        = groupbit;
  param.m_fill          = m_fill;
  param.list.ilist      = ilist;
  param.list.numneigh   = numneigh;
  param.list.firstneigh = firstneigh;
  param.list.inum       = inum;
  param.atom.ntypes     = atom->ntypes;
  param.atom.type       = atom->type;
  param.atom.tag        = atom->tag;
  param.atom.mask       = atom->mask;
  param.atom.nlocal     = atom->nlocal;
  param.atom.x          = static_cast<double(*)[3]>((void*)atom->x[0]);
  param.H.firstnbr      = H.firstnbr;
  param.H.numnbrs       = H.numnbrs;
  param.H.val           = H.val;
  param.H.jlist         = H.jlist;
  param.Tap             = Tap;
  param.maxHlist        = maxHlist;

  //GPTLstart("compute H Full c");
  compute_H_Full_C(&param);
  //GPTLstop("compute H Full c");
  //m_fill = param.m_fill;


  ////compute_H_Full_Newtoff();
  //compute_H_Full();
  //int i;
  //// for (i = 0; i < 10; i ++){
  ////   printf("%6d %f\t", H.jlist[i], H.val[i]);
  //// }
  //// puts("===");
  //// for (i = 0; i < 10; i ++){
  ////   printf("%6d %f\t", H_newtoff.jlist[i], H_newtoff.val[i]);
  //// }
  //// puts("");

  GPTLstop("init matvec compute_H");
  
  for (ii = 0; ii < nn; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {

      /* init pre-conditioner for H and init solution vectors */
      Hdia_inv[i] = 1. / eta[ atom->type[i] ];
      b_s[i]      = -chi[ atom->type[i] ];
      b_t[i]      = -1.0;

      /* linear extrapolation for s & t from previous solutions */
      //s[i] = 2 * s_hist[i][0] - s_hist[i][1];
      //t[i] = 2 * t_hist[i][0] - t_hist[i][1];

      /* quadratic extrapolation for s & t from previous solutions */
      //s[i] = s_hist[i][2] + 3 * ( s_hist[i][0] - s_hist[i][1] );
      t[i] = t_hist[i][2] + 3 * ( t_hist[i][0] - t_hist[i][1]);

      /* cubic extrapolation for s & t from previous solutions */
      s[i] = 4*(s_hist[i][0]+s_hist[i][2])-(6*s_hist[i][1]+s_hist[i][3]);
      //t[i] = 4*(t_hist[i][0]+t_hist[i][2])-(6*t_hist[i][1]+t_hist[i][3]);
    }
  }

  pack_flag = 2;
  comm->forward_comm_fix(this); //Dist_vector( s );
  pack_flag = 3;
  comm->forward_comm_fix(this); //Dist_vector( t );
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::compute_H()
{
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj, flag;
  double dx, dy, dz, r_sqr;
  const double SMALL = 0.0001;

  int *type = atom->type;
  tagint *tag = atom->tag;
  double **x = atom->x;
  int *mask = atom->mask;

  if (reaxc) {
    inum = reaxc->listfull->inum;
    ilist = reaxc->listfull->ilist;
    numneigh = reaxc->listfull->numneigh;
    firstneigh = reaxc->listfull->firstneigh;
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  }

  // fill in the H matrix
  m_fill = 0;
  r_sqr = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      jlist = firstneigh[i];
      jnum = numneigh[i];
      H.firstnbr[i] = m_fill;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];

        dx = x[j][0] - x[i][0];
        dy = x[j][1] - x[i][1];
        dz = x[j][2] - x[i][2];
        r_sqr = SQR(dx) + SQR(dy) + SQR(dz);

        flag = 0;
        if (r_sqr <= SQR(swb)) {
          if (j < n) flag = 1;
          else if (tag[i] < tag[j]) flag = 1;
          else if (tag[i] == tag[j]) {
            if (dz > SMALL) flag = 1;
            else if (fabs(dz) < SMALL) {
              if (dy > SMALL) flag = 1;
              else if (fabs(dy) < SMALL && dx > SMALL)
                flag = 1;
	    }
	  }
	}

        if (flag) {
          H.jlist[m_fill] = j;
          H.val[m_fill] = calculate_H( sqrt(r_sqr), shld[type[i]][type[j]]);
          m_fill++;
        }
      }
      H.numnbrs[i] = m_fill - H.firstnbr[i];
    }
  }

  if (m_fill >= H.m) {
    char str[128];
    sprintf(str,"H matrix size has been exceeded: m_fill=%d H.m=%d\n",
             m_fill, H.m);
    error->warning(FLERR,str);
    error->all(FLERR,"Fix qeq/reax has insufficient QEq matrix size");
  }
}

static inline int llf(double dx, double dy, double dz){
  double SMALL = 0.0001;
  if (dz < -SMALL) return 1;
  if (dz < SMALL){
    if (dy < -SMALL) return 1;
    if (dy < SMALL && dx < -SMALL) return 1;
  }
  return 0;
}
static inline int urb(double dx, double dy, double dz){
  double SMALL = 0.0001;
  if (dz > SMALL) return 1;
  if (dz > -SMALL){
    if (dy > SMALL) return 1;
    if (dy > -SMALL && dx > SMALL) return 1;
  }
  return 0;
}

void FixQEqReaxSunway::compute_H_Full()
{
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj, flag;
  double dx, dy, dz, r_sqr;
  const double SMALL = 0.0001;

  int *type = atom->type;
  tagint *tag = atom->tag;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (reaxc) 
  {
    inum = reaxc->listfull->inum + reaxc->listfull->gnum;
    ilist = reaxc->listfull->ilist;
    numneigh = reaxc->listfull->numneigh;
    firstneigh = reaxc->listfull->firstneigh;
  } 
  else 
  {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  }
  
  m_fill = 0;
  r_sqr = 0;
  double swbsq = swb * swb;

  fix_qeq_pack_t param;
  param.swbsq           = swbsq;
  param.shld            = shld;
  param.groupbit        = groupbit;
  param.m_fill          = m_fill;
  param.list.ilist      = ilist;
  param.list.numneigh   = numneigh;
  param.list.firstneigh = firstneigh;
  param.list.inum       = inum;
  param.atom.ntypes     = atom->ntypes;
  param.atom.type       = atom->type;
  param.atom.tag        = atom->tag;
  param.atom.mask       = atom->mask;
  param.atom.nlocal     = atom->nlocal;
  param.atom.x          = static_cast<double(*)[3]>((void*)atom->x[0]);
  param.H.firstnbr      = H.firstnbr;
  param.H.numnbrs       = H.numnbrs;
  param.H.val           = H.val;
  param.H.jlist         = H.jlist;
  param.Tap             = Tap;
  GPTLstart("compute H Full c");
  compute_H_Full_C(&param);
  GPTLstop("compute H Full c");
  m_fill = param.m_fill;
}

void FixQEqReaxSunway::compute_H_Full_Newtoff()
{
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj, flag;
  double dx, dy, dz, r_sqr;
  const double SMALL = 0.0001;

  int *type = atom->type;
  tagint *tag = atom->tag;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (reaxc) {
    inum = reaxc->listfull->inum;
    ilist = reaxc->listfull->ilist;
    numneigh = reaxc->listfull->numneigh;
    firstneigh = reaxc->listfull->firstneigh;
  } 
  else  {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  }
  
  m_fill = 0;
  r_sqr = 0;
  double swbsq = swb * swb;
  for (ii = 0; ii < inum; ii++) {  
    //i = ilist[ii];
    i = ii;
    if (mask[i] & groupbit) {
      jlist = firstneigh[i];
      jnum = numneigh[i];
      //H_newtoff.firstnbr[i] = m_fill;
  
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
  
        dx = x[j][0] - x[i][0];
        dy = x[j][1] - x[i][1];
        dz = x[j][2] - x[i][2];
        r_sqr = SQR(dx) + SQR(dy) + SQR(dz);
  
        if (r_sqr <= swbsq) {
          //H_newtoff.jlist[m_fill] = j;
          //H_newtoff.val[m_fill] = calculate_H( sqrt(r_sqr), shld[type[i]][type[j]]);
          m_fill++;
          //printf("%d\n", m_fill);
        }//if
      }//for-jj
      //H_newtoff.numnbrs[i] = m_fill - H_newtoff.firstnbr[i];
    }//if
  }//for-ii
  //for (i = 0; i < 10; i ++){
  //  printf("%6d %f\t", H_newtoff.jlist[i], H_newtoff.val[i]);
  //}
  //puts("");

  // fix_qeq_pack_t param;
  // param.swbsq           = swbsq;
  // param.shld            = shld;
  // param.groupbit        = groupbit;
  // param.m_fill          = m_fill;
  // param.list.ilist      = ilist;
  // param.list.numneigh   = numneigh;
  // param.list.firstneigh = firstneigh;
  // param.list.inum       = inum;
  // param.atom.ntypes     = atom->ntypes;
  // param.atom.type       = atom->type;
  // param.atom.tag        = atom->tag;
  // param.atom.mask       = atom->mask;
  // param.atom.nlocal     = atom->nlocal;
  // param.atom.x          = static_cast<double(*)[3]>((void*)atom->x[0]);
  // param.H.firstnbr      = H.firstnbr;
  // param.H.numnbrs       = H.numnbrs;
  // param.H.val           = H.val;
  // param.H.jlist         = H.jlist;
  // param.Tap             = Tap;
  // GPTLstart("compute H Full c");
  // compute_H_Full_C(&param);
  // GPTLstop("compute H Full c");
  //m_fill = param.m_fill;
}

/* ---------------------------------------------------------------------- */

double FixQEqReaxSunway::calculate_H( double r, double gamma)
{
  double Taper, denom;

  Taper = Tap[7] * r + Tap[6];
  Taper = Taper * r + Tap[5];
  Taper = Taper * r + Tap[4];
  Taper = Taper * r + Tap[3];
  Taper = Taper * r + Tap[2];
  Taper = Taper * r + Tap[1];
  Taper = Taper * r + Tap[0];

  denom = r * r * r + gamma;
  denom = pow(denom,0.3333333333333);

  return Taper * EV_TO_KCAL_PER_MOL / denom;
}
/* ------------------------new CG ---------------------------------------- */
int FixQEqReaxSunway::CG_v2( double *b, double *x)
{
  int  i, j, imax;
  double tmp, alpha, beta, b_norm;
  double sig_old, sig_new;
  double deta, heta;
  int nn, jj;
  int *ilist;
  if (reaxc) {
    nn = reaxc->listfull->inum;
    ilist = reaxc->listfull->ilist;
  } else {
    nn = list->inum;
    ilist = list->ilist;
  }
  imax =200;
  pack_flag = 1;
  double dot_local[2], dot_global[2];

  /**********Pipeline PCG*****************/
  fix_qeq_pack_t param;
  param.x          = d;
  param.b          = q;
  param.eta        = eta;
  param.groupbit   = groupbit;
  param.list.inum  = nn;
  param.list.nums  = atom->nlocal + atom->nghost;
  param.list.ilist = ilist;
  param.H.firstnbr = H.firstnbr;
  param.H.numnbrs  = H.numnbrs;
  param.H.jlist    = H.jlist;
  param.H.val      = H.val;
  param.atom.ntypes= atom->ntypes;
  param.atom.mask  = atom->mask;
  param.atom.type  = atom->type;

  //int myrank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //double time0, time1, time2, time3;
  //time0 = MPI_Wtime();

  GPTLstart("CG before for:");
  GPTLstart("CG matrix 1");
  sparse_matvec( &H, x, q);             //q = H * x;
  GPTLstop("CG matrix 1");
  
  GPTLstart("CG cmp 1");
  for (jj = 0; jj < nn; ++jj) 
  {
    j = ilist[jj];
    if (atom->mask[j] & groupbit)
    {
      r[j] = b[j] - 1 * q[j];
      u[j] = r[j] * Hdia_inv[j];        //u = H' * r;
    }
    d[j] = u[j];                        //d = u;
  }
  GPTLstop("CG cmp 1");

  GPTLstart("CG comm 1");
  comm->forward_comm_fix(this); 
  GPTLstop("CG comm 1");

  GPTLstart("CG matrix 2");
  sparse_matvec( &H, d, q );            //w = q  = H' * u;
  GPTLstop("CG matrix 2");

  GPTLstart("CG cmp 2");
  for (jj = 0; jj < nn; ++jj) 
  {
    j = ilist[jj];
    w[j] = q[j];                        //w = q;
    if (atom->mask[j] & groupbit)
    {
      m[j] = q[j] * Hdia_inv[j];        //m = H' * q;
    }
    d[j] = m[j];                        //d = m;
  }                                     
  GPTLstop("CG cmp 2");

  GPTLstart("CG comm 2");
  comm->forward_comm_fix(this); 
  GPTLstop("CG comm 2");

  GPTLstart("CG matrix 3");
  sparse_matvec_C_spawn(&param);
  
  double my_sum[3], my_res[3];
  my_sum[0] = my_sum[1] = my_sum[2] = 0.0;
  my_res[0] = my_res[1] = my_res[2] = 0.0;

  for (jj = 0; jj < nn; ++jj) 
  {
    j = ilist[jj];
    p[j] = u[j];                        //p = u
    ss[j] = w[j];                       //s = w
    v[j] = m[j];                        //v(q) = m;
    if (atom->mask[j] & groupbit)
    {
      my_sum[0] += SQR( b[j]);
      my_sum[1] += u[j] * r[j];
      my_sum[2] += u[j] * w[j];
    }
  }                                     
  MPI_Allreduce(my_sum, my_res, 3, MPI_DOUBLE, MPI_SUM, world);
  b_norm  = sqrt(my_res[0]);
  sig_old = my_res[1];
  deta    = my_res[2];
  heta    = deta;                       //heta = deta;
  //alpha   = sig_old / heta;             //alpha = sig_old / heta;
  alpha   = sig_old / deta;             //alpha = sig_old / heta;
  //sig_new = sig_old;
  dot_global[0] = sig_old;

  sparse_matvec_C_join();
  GPTLstop("CG matrix 3");

  for (jj = 0; jj < nn; ++jj) 
  {
    z[jj] = q[jj];                        //z = q(n) 
  }                                     

  GPTLstop("CG before for:");
  
  //time1 = MPI_Wtime();
  GPTLstart("CG all for:");
  //dot_global[0] = sig_new;
  for (i = 1; i < imax && sqrt(dot_global[0]) / b_norm > tolerance; ++i) 
  {
    GPTLstart("cmp1 in for");
    dot_local[0] = dot_local[1] = 0.0;
    for (jj = 0; jj < nn; jj ++)
    {
      r[jj] -= alpha * ss[jj];
      u[jj] -= alpha * v[jj];
      w[jj] -= alpha * z[jj];
      dot_local[0] += u[jj] * r[jj];
      dot_local[1] += u[jj] * w[jj];      
      if (atom->mask[jj] & groupbit)
        d[jj] = w[jj] * Hdia_inv[jj];     //m = H' * w;
    }
    GPTLstop("cmp1 in for");

    GPTLstart("comm in for");
    comm->forward_comm_fix(this); 
    GPTLstop("comm in for");

    GPTLstart("sparse matvec in for");
    sparse_matvec_C_spawn(&param);
    MPI_Allreduce( dot_local, dot_global, 2, MPI_DOUBLE, MPI_SUM, world);
    beta = dot_global[0] / sig_old;         //beta = sig_new / sig_old;
    heta = dot_global[1] - beta * beta * heta;
    for (jj = 0; jj < nn; jj ++)
    {
      x[jj]   += alpha * p[jj];
      p[jj]   = u[jj] + p[jj] * beta;
      ss[jj]  = w[jj] + ss[jj] * beta;
      v[jj]   = d[jj] + v[jj] * beta;
    }
    alpha = dot_global[0] / heta;           //alpha = sig_new / heta;
    sparse_matvec_C_join();
    GPTLstop("sparse matvec in for");
    
    GPTLstart("cmp2 in for");
    for (jj = 0; jj < nn; jj ++)
    {
      z[jj] = q[jj] + z[jj] * beta;
    }
    
    sig_old = dot_global[0];
    GPTLstop("cmp2 in for");
  }
  GPTLstop("CG all for:");
  //time2 = MPI_Wtime();
  //printf("tot:%lf, %lf, %lf\n", time2-time0, time1-time0, time2-time1);

  if (i >= imax && comm->me == 0) 
  {
    char str[128];
    sprintf(str,"Fix qeq/reax CG convergence failed after %d iterations "
            "at " BIGINT_FORMAT " step",i,update->ntimestep);
    error->warning(FLERR,str);
  }
  return i;
}

int FixQEqReaxSunway::CG_new( double *b, double *x)
{
  int  i, j, imax;
  double tmp, alpha, beta, b_norm;
  double sig_old, sig_new;
  double deta, heta;
  int nn, jj;
  int *ilist;
  if (reaxc) {
    nn = reaxc->listfull->inum;
    ilist = reaxc->listfull->ilist;
  } else {
    nn = list->inum;
    ilist = list->ilist;
  }
  imax =200;
  pack_flag = 1;

  ///**********Pipeline PCG*****************/
  ////double *m, *u, *v, *w, *z, *s, *n_test;
  //double *m, *u, *v, *w, *z, *ss, *n_test;
  //int sz = reaxc->listfull->inum + reaxc->listfull->gnum;
  //if (dual_enabled) sz*= 2;
  //memory->create(m,sz,"qeq:m");//Pipeline PCG
  //memory->create(v,sz,"qeq:v");
  //memory->create(u,sz,"qeq:u");
  //memory->create(w,sz,"qeq:w");
  //memory->create(z,sz,"qeq:z");
  //memory->create(ss,sz,"qeq:ss");
  //memory->create(n_test,sz,"qeq:p_test");

  fix_qeq_pack_t param;
  param.x          = d;
  param.b          = q;
  param.eta        = eta;
  param.groupbit   = groupbit;
  param.list.inum  = nn;
  param.list.nums  = atom->nlocal + atom->nghost;
  param.list.ilist = ilist;
  param.H.firstnbr = H.firstnbr;
  param.H.numnbrs  = H.numnbrs;
  param.H.jlist    = H.jlist;
  param.H.val      = H.val;
  param.atom.ntypes= atom->ntypes;
  param.atom.mask  = atom->mask;
  param.atom.type  = atom->type;

  //int myrank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //double time0, time1, time2, time3;
  //time0 = MPI_Wtime();

  //printf("pack_flag = %d\n", pack_flag);

  sparse_matvec( &H, x, q);             //q = H * x;

  //vector_sum( r , 1.,  b, -1., q, nn);  //r = b -1 * q;
  for (jj = 0; jj < nn; ++jj) 
  {
    j = ilist[jj];
    if (atom->mask[j] & groupbit)
    {
      
      r[j] = b[j] - 1 * q[j];
      u[j] = r[j] * Hdia_inv[j];        //u = H' * r;
    }
    d[j] = u[j];                        //d = u;
  }

  comm->forward_comm_fix(this); 
  sparse_matvec( &H, d, q );            //w = q  = H' * u;
  //comm->reverse_comm_fix(this); 
  for (jj = 0; jj < nn; ++jj) 
  {
    j = ilist[jj];
    w[j] = q[j];                        //w = q;
    if (atom->mask[j] & groupbit)
    {
      m[j] = q[j] * Hdia_inv[j];        //m = H' * q;
    }
    d[j] = m[j];                        //d = m;
  }                                     

  comm->forward_comm_fix(this); 
  sparse_matvec( &H, d, q);            //v(n) = q = H * m;
  //comm->reverse_comm_fix(this);       
  
  double my_sum[3], my_res[3];
  my_sum[0] = my_sum[1] = my_sum[2] = 0.0;
  my_res[0] = my_res[1] = my_res[2] = 0.0;

  for (jj = 0; jj < nn; ++jj) 
  {
    j = ilist[jj];
    p[j] = u[j];                      //p = u
    ss[j] = w[j];                      //s = w
    v[j] = m[j];                      //v(q) = m;
    z[j] = q[j];                      //z = q(n) 
    if (atom->mask[j] & groupbit)
    {
      my_sum[0] += SQR( b[j]);
      my_sum[1] += u[j] * r[j];
      my_sum[2] += u[j] * w[j];
    }
  }                                     
  MPI_Allreduce(my_sum, my_res, 3, MPI_DOUBLE, MPI_SUM, world);
  b_norm  = sqrt(my_res[0]);
  sig_old = my_res[1];
  deta    = my_res[2];

  //b_norm = parallel_norm(b, nn);
  //sig_old = parallel_dot(u, r, nn);     //sig_old = dot(u, r);
  //deta = parallel_dot(u, w, nn);        //deta = dot(u, w);
  heta = deta;                           //heta = deta;
  alpha = sig_old / heta;                //alpha = sig_old / heta;
  sig_new = sig_old;

  //time1 = MPI_Wtime();

  for (i = 1; i < imax && sqrt(sig_new) / b_norm > tolerance; ++i) 
  //for (i = 1; i < imax; ++i) 
  {
    double dot_local[2], dot_global[2];
    dot_local[0] = dot_local[1] = 0.0;
    for (jj = 0; jj < nn; jj ++)
    {
      x[jj] += alpha * p[jj];
      r[jj] -= alpha * ss[jj];
      u[jj] -= alpha * v[jj];
      w[jj] -= alpha * z[jj];
      dot_local[0] += u[jj] * r[jj];
      dot_local[1] += u[jj] * w[jj];      
      if (atom->mask[jj] & groupbit)
        d[jj] = w[jj] * Hdia_inv[jj];        //m = H' * w;

    }
    // // vector_add(x,  alpha, p, nn );        //x = x + alpha * d;
    // // vector_add(r, -alpha, s, nn );        //r = r - alpha * u;
    // // vector_add(u, -alpha, v, nn );        //p = p - alpha * v;
    // // vector_add(w, -alpha, z, nn );        //w = w - alpha * z;
    // sig_new=parallel_dot(u, r, nn);       //sig_new=dot(p, r);
    // deta = parallel_dot(u, w, nn);        //deta = dot(p, w);

    //for (jj = 0; jj < nn; ++jj) 
    //{
    //  if (atom->mask[jj] & groupbit)
    //    d[jj] = w[jj] * Hdia_inv[jj];        //m = H' * w;
    //  //d[j] = m[j];
    //}                                     
    
    //for (jj = 0; jj < nn; jj ++)
    //{
    //  dot_local[0] += u[jj] * r[jj];
    //  dot_local[1] += u[jj] * w[jj];
    //}

    comm->forward_comm_fix(this); 

    GPTLstart("sparse matvec in for");
    sparse_matvec_C_spawn(&param);

    //for (jj = 0; jj < nn; jj ++)
    //{
    //  dot_local[0] += u[jj] * r[jj];
    //  dot_local[1] += u[jj] * w[jj];
    //}
    MPI_Allreduce( dot_local, dot_global, 2, MPI_DOUBLE, MPI_SUM, world);
    sig_new = dot_global[0];
    deta = dot_global[1];
    beta = sig_new / sig_old;             //beta = sig_new / sig_old;
    heta = deta - beta * beta * heta;
    alpha = sig_new / heta;               //alpha = sig_new / heta;
    for (jj = 0; jj < nn; jj ++)
    {
      p[jj] = u[jj] + p[jj] * beta;
      ss[jj] = w[jj] + ss[jj] * beta;
      v[jj] = d[jj] + v[jj] * beta;
      //z[jj] = q[jj] + z[jj] * beta;
    }

    sparse_matvec_C_join();
    GPTLstop("sparse matvec in for");

    //sparse_matvec( &H, d, q);            //q = H * m;
    //comm->reverse_comm_fix(this); 
    
    for (jj = 0; jj < nn; jj ++)
    {
      // p[jj] = u[jj] + p[jj] * beta;
      // s[jj] = w[jj] + s[jj] * beta;
      // v[jj] = d[jj] + v[jj] * beta;
      z[jj] = q[jj] + z[jj] * beta;
    }
    // vector_sum( p, 1., u, beta, p, nn );  //d = p+ beta * d;
    // vector_sum( s, 1., w, beta, s, nn );  //u = w+ beta * u;
    // vector_sum( v, 1., d, beta, v, nn );  //v = m+ beta * v;
    // vector_sum( z, 1., q, beta, z, nn );  //z = q+ beta * z;
    sig_old = sig_new;
  }
  //time2 = MPI_Wtime();
  //printf("tot:%lf, %lf, %lf\n", time2-time0, time1-time0, time2-time1);


  //memory->destroy(m);
  //memory->destroy(v);
  //memory->destroy(u);
  //memory->destroy(w);
  //memory->destroy(z);
  //memory->destroy(ss);
  //memory->destroy(n_test);

  if (i >= imax && comm->me == 0) 
  {
    char str[128];
    sprintf(str,"Fix qeq/reax CG convergence failed after %d iterations "
            "at " BIGINT_FORMAT " step",i,update->ntimestep);
    error->warning(FLERR,str);
  }
  return i;
}

/* ---------------------------------------------------------------------- */
int FixQEqReaxSunway::CG( double *b, double *x)
{
  int  i, j, imax;
  double tmp, alpha, beta, b_norm;
  double sig_old, sig_new;
  double deta, heta;
  int nn, jj;
  int *ilist;
  if (reaxc) {
    nn = reaxc->listfull->inum;
    ilist = reaxc->listfull->ilist;
  } else {
    nn = list->inum;
    ilist = list->ilist;
  }
  imax =200;
  pack_flag = 1;

  ///**********Pipeline PCG*****************/
  double *m, *u, *v, *w, *z, *s, *n_test;
  int sz = reaxc->listfull->inum + reaxc->listfull->gnum;
  if (dual_enabled) sz*= 2;
  memory->create(m,sz,"qeq:m");//Pipeline PCG
  memory->create(v,sz,"qeq:v");
  memory->create(u,sz,"qeq:u");
  memory->create(w,sz,"qeq:w");
  memory->create(z,sz,"qeq:z");
  memory->create(s,sz,"qeq:s");
  memory->create(n_test,sz,"qeq:p_test");

  fix_qeq_pack_t param;
  param.x          = d;
  param.b          = q;
  param.eta        = eta;
  param.groupbit   = groupbit;
  param.list.inum  = nn;
  param.list.nums  = atom->nlocal + atom->nghost;
  param.list.ilist = ilist;
  param.H.firstnbr = H.firstnbr;
  param.H.numnbrs  = H.numnbrs;
  param.H.jlist    = H.jlist;
  param.H.val      = H.val;
  param.atom.ntypes= atom->ntypes;
  param.atom.mask  = atom->mask;
  param.atom.type  = atom->type;

  sparse_matvec( &H, x, q);             //q = H * x;

  vector_sum( r , 1.,  b, -1., q, nn);  //r = b -1 * q;
  for (jj = 0; jj < nn; ++jj) 
  {
    j = ilist[jj];
    if (atom->mask[j] & groupbit)
    {
      u[j] = r[j] * Hdia_inv[j];        //u = H' * r;
    }
    d[j] = u[j];                        //d = u;
  }

  comm->forward_comm_fix(this); 
  sparse_matvec( &H, d, q );            //w = q  = H' * u;
  //comm->reverse_comm_fix(this); 
  for (jj = 0; jj < nn; ++jj) 
  {
    j = ilist[jj];
    w[j] = q[j];                        //w = q;
    if (atom->mask[j] & groupbit)
    {
      m[j] = q[j] * Hdia_inv[j];        //m = H' * q;
    }
    d[j] = m[j];                        //d = m;
  }                                     

  comm->forward_comm_fix(this); 
  sparse_matvec( &H, d, q);            //v(n) = q = H * m;
  //comm->reverse_comm_fix(this);       
  
  for (jj = 0; jj < nn; ++jj) 
  {
    j = ilist[jj];
    p[j] = u[j];                      //p = u
    s[j] = w[j];                      //s = w
    v[j] = m[j];                      //v(q) = m;
    z[j] = q[j];                      //z = q(n) 
  }                                     

  b_norm = parallel_norm(b, nn);
  sig_old = parallel_dot(u, r, nn);     //sig_old = dot(u, r);
  deta = parallel_dot(u, w, nn);        //deta = dot(u, w);
  heta = deta;                           //heta = deta;
  alpha = sig_old / heta;                //alpha = sig_old / heta;
  sig_new = sig_old;

  for (i = 1; i < imax && sqrt(sig_new) / b_norm > tolerance; ++i) 
  //for (i = 1; i < imax; ++i) 
  {
    double dot_local[2], dot_global[2];
    dot_local[0] = dot_local[1] = 0.0;
    for (jj = 0; jj < nn; jj ++){
      x[jj] += alpha * p[jj];
      r[jj] -= alpha * s[jj];
      u[jj] -= alpha * v[jj];
      w[jj] -= alpha * z[jj];
    }
    // // vector_add(x,  alpha, p, nn );        //x = x + alpha * d;
    // // vector_add(r, -alpha, s, nn );        //r = r - alpha * u;
    // // vector_add(u, -alpha, v, nn );        //p = p - alpha * v;
    // // vector_add(w, -alpha, z, nn );        //w = w - alpha * z;
    // sig_new=parallel_dot(u, r, nn);       //sig_new=dot(p, r);
    // deta = parallel_dot(u, w, nn);        //deta = dot(p, w);

    for (jj = 0; jj < nn; ++jj) 
    {
      if (atom->mask[jj] & groupbit)
        d[jj] = w[jj] * Hdia_inv[jj];        //m = H' * w;
      //d[j] = m[j];
    }                                     
    
    comm->forward_comm_fix(this); 
    sparse_matvec_C_spawn(&param);

    for (jj = 0; jj < nn; jj ++){
      dot_local[0] += u[jj] * r[jj];
      dot_local[1] += u[jj] * w[jj];
    }
    MPI_Allreduce( dot_local, dot_global, 2, MPI_DOUBLE, MPI_SUM, world);
    sig_new = dot_global[0];
    deta = dot_global[1];
    beta = sig_new / sig_old;             //beta = sig_new / sig_old;
    heta = deta - beta * beta * heta;
    alpha = sig_new / heta;               //alpha = sig_new / heta;
    for (jj = 0; jj < nn; jj ++){
      p[jj] = u[jj] + p[jj] * beta;
      s[jj] = w[jj] + s[jj] * beta;
      v[jj] = d[jj] + v[jj] * beta;
    }

    sparse_matvec_C_join();
    //sparse_matvec( &H, d, q);            //q = H * m;
    //comm->reverse_comm_fix(this); 
    
    for (jj = 0; jj < nn; jj ++){
      // p[jj] = u[jj] + p[jj] * beta;
      // s[jj] = w[jj] + s[jj] * beta;
      // v[jj] = d[jj] + v[jj] * beta;
      z[jj] = q[jj] + z[jj] * beta;
    }
    // vector_sum( p, 1., u, beta, p, nn );  //d = p+ beta * d;
    // vector_sum( s, 1., w, beta, s, nn );  //u = w+ beta * u;
    // vector_sum( v, 1., d, beta, v, nn );  //v = m+ beta * v;
    // vector_sum( z, 1., q, beta, z, nn );  //z = q+ beta * z;
    sig_old = sig_new;
  }
  memory->destroy(m);
  memory->destroy(v);
  memory->destroy(u);
  memory->destroy(w);
  memory->destroy(z);
  memory->destroy(s);
  memory->destroy(n_test);

  if (i >= imax && comm->me == 0) {
    char str[128];
    sprintf(str,"Fix qeq/reax CG convergence failed after %d iterations "
            "at " BIGINT_FORMAT " step",i,update->ntimestep);
    error->warning(FLERR,str);
  }
  return i;
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::sparse_matvec( sparse_matrix *A, double *x, double *b)
{
  int i, j, itr_j;
  int nn, NN, ii;
  int *ilist;
  if (reaxc) 
  {
    nn = reaxc->listfull->inum;
    NN = reaxc->listfull->inum + reaxc->listfull->gnum;
    ilist = reaxc->listfull->ilist;
  } 
  else 
  {
    nn = list->inum;
    NN = list->inum + list->gnum;
    ilist = list->ilist;
  }

  fix_qeq_pack_t param;
  param.x          = x;
  param.b          = b;
  param.eta        = eta;
  param.groupbit   = groupbit;
  param.list.inum  = nn;
  param.list.nums  = NN;
  param.list.ilist = ilist;
  param.H.firstnbr = A->firstnbr;
  param.H.numnbrs  = A->numnbrs;
  param.H.jlist    = A->jlist;
  param.H.val      = A->val;
  param.atom.ntypes=atom->ntypes;
  param.atom.mask  = atom->mask;
  param.atom.type  = atom->type;
  GPTLstart("slave sparse matvec c");
  sparse_matvec_C(&param);
  GPTLstop("slave sparse matvec c");
  return;

  GPTLstart("slave sparse matvec c");
  //for (ii = 0; ii < NN; ++ii) 
  for (ii = 0; ii < nn; ++ii) 
  {
    i = ilist[ii];
    //i = ii;
    if(atom->mask[i] & groupbit)
    {
      if(ii < nn) b[i] = eta[atom->type[i]] * x[i];
      else  b[i] = 0;
    }
    
    //if(A->numnbrs[i] > 512)
    //  printf("numnbrs[%d] = %d\n", i, A->numnbrs[i]);

    for(itr_j=A->firstnbr[i]; itr_j < A->firstnbr[i]+A->numnbrs[i]; itr_j++) 
    {
      j = A->jlist[itr_j];
      b[i] += A->val[itr_j] * x[j];
    }//for-itr_j
  }//for-ii
  GPTLstop("slave sparse matvec c");
}

void FixQEqReaxSunway::sparse_matvec_noff( sparse_matrix *A, double *x, double *b)
{
  int i, j, itr_j;
  int nn, NN, ii;
  int *ilist;
  if (reaxc) 
  {
    nn = reaxc->listfull->inum;
    NN = reaxc->listfull->inum + reaxc->listfull->gnum;
    ilist = reaxc->listfull->ilist;
  } 
  else 
  {
    nn = list->inum;
    NN = list->inum + list->gnum;
    ilist = list->ilist;
  }


  for (ii = 0; ii < nn; ++ii) 
  {
    i = ilist[ii];
    //i = ii;
    if(atom->mask[i] & groupbit)
    {
      if(ii < nn) b[i] = eta[atom->type[i]] * x[i];
      else  b[i] = 0;
    }
    for(itr_j=A->firstnbr[i]; itr_j < A->firstnbr[i]+A->numnbrs[i]; itr_j++) 
    {
      j = A->jlist[itr_j];
      b[i] += A->val[itr_j] * x[j];
    }//for-itr_j
  }//for-ii

  //for (ii = 0; ii < nn; ++ii) 
  //{
  //  i = ilist[ii];
  //  if (atom->mask[i] & groupbit) b[i] = eta[atom->type[i]] * x[i];
  //}
  //for (ii = nn; ii < NN; ++ii) 
  //{
  //  i = ilist[ii];
  //  if (atom->mask[i] & groupbit) b[i] = 0;
  //}
  //for (ii = 0; ii < NN; ++ii) 
  //{
  //  i = ilist[ii];
  //  if (atom->mask[i] & groupbit) 
  //  {
  //    for (itr_j=A->firstnbr[i]; itr_j < A->firstnbr[i]+A->numnbrs[i]; itr_j++) 
  //    {
  //      j = A->jlist[itr_j];
  //      b[i] += A->val[itr_j] * x[j];
  //    }//for-itr_j
  //  }//if
  //}//for-ii

}

void FixQEqReaxSunway::sparse_matvec_comm(sparse_matrix *H, double *x, double *b){
  int pack_flag_orig = pack_flag;
  pack_flag = PACK_SPARSE_MATVEC;
  sparse_x = x;
  sparse_b = b;
  comm->forward_comm_fix(this);
  sparse_matvec(H, x, b);
  comm->reverse_comm_fix(this);
  pack_flag = pack_flag_orig;
}
/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::calculate_Q()
{
  int i, k;
  double u, s_sum, t_sum;
  double *q = atom->q;

  int nn, ii;
  int *ilist;

  if (reaxc) {
    nn = reaxc->listfull->inum;
    ilist = reaxc->listfull->ilist;
  } else {
    nn = list->inum;
    ilist = list->ilist;
  }

  //s_sum = parallel_vector_acc( s, nn);
  //t_sum = parallel_vector_acc( t, nn);
  //u = s_sum / t_sum;
  
  double my_acc[2], res[2];
  my_acc[0] = my_acc[1] = 0.0;
  res[0] = res[1] = 0.0;
  for (ii = 0; ii < nn; ++ii) 
  {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
    {
      my_acc[0] += s[i];
      my_acc[1] += t[i];
    }
  }

  MPI_Allreduce(my_acc, res, 2, MPI_DOUBLE, MPI_SUM, world);
  s_sum = res[0];
  t_sum = res[1];
  u = s_sum / t_sum;

  for (ii = 0; ii < nn; ++ii) 
  {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) 
    {
      q[i] = s[i] - u * t[i];

      /* backup s & t */
      for (k = 4; k > 0; --k) {
        s_hist[i][k] = s_hist[i][k-1];
        t_hist[i][k] = t_hist[i][k-1];
      }
      s_hist[i][0] = s[i];
      t_hist[i][0] = t[i];
    }
  }

  pack_flag = 4;
  comm->forward_comm_fix(this); //Dist_vector( atom->q );
}

/* ---------------------------------------------------------------------- */

int FixQEqReaxSunway::pack_forward_comm(int n, int *list, double *buf,
                                  int pbc_flag, int *pbc)
{
  int m;

  if (pack_flag == 1)
    for(m = 0; m < n; m++) buf[m] = d[list[m]];
  else if (pack_flag == 2)
    for(m = 0; m < n; m++) buf[m] = s[list[m]];
  else if (pack_flag == 3)
    for(m = 0; m < n; m++) buf[m] = t[list[m]];
  else if (pack_flag == 4)
    for(m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  else if (pack_flag == PACK_SPARSE_MATVEC)
    for(m = 0; m < n; m++) buf[m] = sparse_x[list[m]];
  else if (pack_flag == 5) 
  {
    m = 0;
    for(int i = 0; i < n; i++) {
      int j = 2 * list[i];
      buf[m++] = d[j  ];
      buf[m++] = d[j+1];
    }
    return m;
  }
  return n;
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1)
    for(m = 0, i = first; m < n; m++, i++) d[i] = buf[m];
  else if (pack_flag == 2)
    for(m = 0, i = first; m < n; m++, i++) s[i] = buf[m];
  else if (pack_flag == 3)
    for(m = 0, i = first; m < n; m++, i++) t[i] = buf[m];
  else if (pack_flag == 4)
    for(m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
  else if (pack_flag == PACK_SPARSE_MATVEC)
    for (m = 0, i = first; m < n; m ++, i ++) sparse_x[i] = buf[m];
  else if (pack_flag == 5) {
    int last = first + n;
    m = 0;
    for(i = first; i < last; i++) {
      int j = 2 * i;
      d[j  ] = buf[m++];
      d[j+1] = buf[m++];
    }
  }

}

/* ---------------------------------------------------------------------- */

int FixQEqReaxSunway::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  if (pack_flag == PACK_SPARSE_MATVEC){
    for (m = 0, i = first; m < n; m ++, i ++) buf[m] = sparse_b[i];
  } else if (pack_flag == 5) {
    m = 0;
    int last = first + n;
    for(i = first; i < last; i++) {
      int indxI = 2 * i;
      buf[m++] = q[indxI  ];
      buf[m++] = q[indxI+1];
    }
    return m;
  } else {
    for (m = 0, i = first; m < n; m++, i++) buf[m] = q[i];
    return n;
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m,i;
  if (pack_flag == PACK_SPARSE_MATVEC){
    for (m = 0; m < n; m ++) sparse_b[list[m]] += buf[m];
  } else if (pack_flag == 5) {
    int m = 0;
    for(int i = 0; i < n; i++) {
      int indxI = 2 * list[i];
      q[indxI  ] += buf[m++];
      q[indxI+1] += buf[m++];
    }
  } else {
    for (int m = 0; m < n; m++) q[list[m]] += buf[m];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixQEqReaxSunway::memory_usage()
{
  double bytes;

  bytes = atom->nmax*nprev*2 * sizeof(double); // s_hist & t_hist
  bytes += atom->nmax*11 * sizeof(double); // storage
  bytes += n_cap*2 * sizeof(int); // matrix...
  bytes += m_cap * sizeof(int);
  bytes += m_cap * sizeof(double);

  if (dual_enabled)
    bytes += atom->nmax*4 * sizeof(double); // double size for q, d, r, and p

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

void FixQEqReaxSunway::grow_arrays(int nmax)
{
  memory->grow(s_hist,nmax,nprev,"qeq:s_hist");
  memory->grow(t_hist,nmax,nprev,"qeq:t_hist");
}

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

void FixQEqReaxSunway::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < nprev; m++) {
    s_hist[j][m] = s_hist[i][m];
    t_hist[j][m] = t_hist[i][m];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixQEqReaxSunway::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nprev; m++) buf[m] = s_hist[i][m];
  for (int m = 0; m < nprev; m++) buf[nprev+m] = t_hist[i][m];
  return nprev*2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixQEqReaxSunway::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nprev; m++) s_hist[nlocal][m] = buf[m];
  for (int m = 0; m < nprev; m++) t_hist[nlocal][m] = buf[nprev+m];
  return nprev*2;
}

/* ---------------------------------------------------------------------- */

double FixQEqReaxSunway::parallel_norm( double *v, int n)
{
  int  i;
  double my_sum, norm_sqr;

  int ii;
  int *ilist;

  if (reaxc)
    ilist = reaxc->listfull->ilist;
  else
    ilist = list->ilist;

  my_sum = 0.0;
  norm_sqr = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_sum += SQR( v[i]);
  }

  MPI_Allreduce( &my_sum, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world);

  return sqrt( norm_sqr);
}

/* ---------------------------------------------------------------------- */

double FixQEqReaxSunway::parallel_dot( double *v1, double *v2, int n)
{
  int  i;
  double my_dot, res;

  int ii;
  int *ilist;

  if (reaxc)
    ilist = reaxc->listfull->ilist;
  else
    ilist = list->ilist;

  my_dot = 0.0;
  res = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
    {
      my_dot += v1[i] * v2[i];
    }
  }
  MPI_Allreduce( &my_dot, &res, 1, MPI_DOUBLE, MPI_SUM, world);
  return res;
}

/* ---------------------------------------------------------------------- */

double FixQEqReaxSunway::parallel_vector_acc( double *v, int n)
{
  int  i;
  double my_acc, res;

  int ii;
  int *ilist;

  if (reaxc)
    ilist = reaxc->listfull->ilist;
  else
    ilist = list->ilist;

  my_acc = 0.0;
  res = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_acc += v[i];
  }

  MPI_Allreduce( &my_acc, &res, 1, MPI_DOUBLE, MPI_SUM, world);

  return res;
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::vector_sum( double* dest, double c, double* v,
                                double d, double* y, int k)
{
  int kk;
  int *ilist;
  int nn = k;

  if (reaxc)
    ilist = reaxc->listfull->ilist;
  else
    ilist = list->ilist;

  //for (--k; k>=0; --k) 
  for (k=0; k<nn; k++) 
  {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] = c * v[kk] + d * y[kk];
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxSunway::vector_add( double* dest, double c, double* v, int k)
{
  int kk;
  int *ilist;

  if (reaxc)
    ilist = reaxc->listfull->ilist;
  else
    ilist = list->ilist;

  for (--k; k>=0; --k) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] += c * v[kk];
  }
}
