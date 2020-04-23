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

#include <stdio.h>
#include <string.h>
#include "fix_nve_sunway.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "sunway.h"
#include "gptl.h"
#include "fix_nve_sw64.h"
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVESunway::FixNVESunway(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"nve/sphere") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix nve command");

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVESunway::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVESunway::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVESunway::initial_integrate(int vflag)
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  GPTLstart("fix nve initial integratesunway");
  fix_nve_param_t pm;
  pm.x = (double(*)[3])(void*)atom->x[0];
  pm.v = (double(*)[3])(void*)atom->v[0];
  pm.f = (double(*)[3])(void*)atom->f[0];
  pm.rmass = atom->rmass;
  pm.mass = atom->mass;
  pm.dtv = dtv;
  pm.dtf = dtf;
  pm.type = atom->type;
  pm.mask = atom->mask;
  pm.nlocal = nlocal;
  pm.groupbit = groupbit;
  pm.ntypes = atom->ntypes;
  fix_nve_initial_integrate(&pm);
  GPTLstop("fix nve initial integratesunway");
}

/* ---------------------------------------------------------------------- */

void FixNVESunway::final_integrate()
{
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  GPTLstart("fix nve final integrate sunway");
  fix_nve_param_t pm;
  pm.x = (double(*)[3])(void*)atom->x[0];
  pm.v = (double(*)[3])(void*)atom->v[0];
  pm.f = (double(*)[3])(void*)atom->f[0];
  pm.rmass = atom->rmass;
  pm.mass = atom->mass;
  pm.dtv = dtv;
  pm.dtf = dtf;
  pm.type = atom->type;
  pm.mask = atom->mask;
  pm.nlocal = nlocal;
  pm.groupbit = groupbit;
  pm.ntypes = atom->ntypes;
  fix_nve_final_integrate(&pm);
  GPTLstop("fix nve final integrate sunway");
}

/* ---------------------------------------------------------------------- */

void FixNVESunway::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVESunway::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVESunway::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
