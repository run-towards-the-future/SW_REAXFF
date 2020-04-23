#include "sunway.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "reaxc_defs_sunway.h"
#include "reaxc_ctypes_sunway.h"
#include "reaxc_inlines_sw64.h"
#include "simd.h"

static inline int llf(double *dvec)
{
  double SMALL = 0.0001;
  if(dvec[2] < -SMALL || (dvec[2]<SMALL && (dvec[1] < -SMALL || 
    (dvec[1] > -SMALL && dvec[0]>SMALL))))
    return 1;
  return 0;
}
static inline int urb(double *dvec)
{
  double SMALL = 0.0001;
  if(dvec[2] > SMALL || (dvec[2] > -SMALL && (dvec[1] > SMALL || 
    (dvec[1] > -SMALL && dvec[0] > SMALL))))
    return 1;
  return 0;
}


#ifdef MPE
#include <mpi.h>
#include <athread.h>
//#define LWPF_UNITS U(NONBOND)
//#include "lwpf2.h"
extern SLAVE_FUN(vdW_Coulomb_Energy_Full_C_para_general)(vdw_coulomb_pack_t *);
extern SLAVE_FUN(vdW_Coulomb_Energy_Full_C_para_special)(vdw_coulomb_pack_t *);
extern SLAVE_FUN(vdW_Coulomb_Energy_cpe_evflag1)(vdw_coulomb_pack_t *);
extern SLAVE_FUN(vdW_Coulomb_Energy_cpe_evflag0)(vdw_coulomb_pack_t *);
extern void Validate_Lists_C(vdw_coulomb_pack_t * param);
extern void sort_bonds_C(vdw_coulomb_pack_t * param); 
int r = 0;
void vdW_Coulomb_Energy_Full_C_test_err(vdw_coulomb_pack_t * param)
{
  reax_system_c *system       = param->system;
  control_params *control     = param->control;
  simulation_data *data       = param->data;
  storage *workspace          = param->workspace;
  reax_list **lists           = param->lists;

  int i, j, pj, natoms, nlocal;
  int start_i, end_i, flag;
  rc_tagint orig_i, orig_j;
  double p_vdW1, p_vdW1i;
  double powr_vdW1, powgi_vdW1;
  double tmp, r_ij, fn13, exp1, exp2;
  double Tap, dTap, dfn13, CEvd, CEclmb, de_core;
  double dr3gamij_1, dr3gamij_3;
  double e_ele, e_vdW, e_core, SMALL = 0.0001;
  double e_lg, de_lg, r_ij5, r_ij6, re6;
  rvec temp, ext_press;
  two_body_parameters *twbp;
  far_neighbor_data_full *nbr_pj;
  reax_list *far_nbrs;

  // Tallying variables:
  double pe_vdw, f_tmp, delij[3];

  natoms = system->N;
  nlocal = system->n;
  far_nbrs = (*lists) + FAR_NBRS_FULL;
  p_vdW1 = system->reax_param.gp.l[28];
  p_vdW1i = 1.0 / p_vdW1;
  e_core = 0;
  e_vdW = 0;
  e_lg = de_lg = 0.0;

  //for( i = 0; i < natoms; ++i ) 
  for( i = 0; i < nlocal; ++i ) 
  {
    //if (system->my_atoms[i].type < 0) continue;
    if (system->packed_atoms[i].type < 0) continue;
    start_i = Start_Index(i, far_nbrs);
    end_i   = End_Index(i, far_nbrs);
    //orig_i  = system->my_atoms[i].orig_id;
    orig_i  = system->packed_atoms[i].orig_id;

    for( pj = start_i; pj < end_i; ++pj ) 
    {
      nbr_pj = &(far_nbrs->select.far_nbr_list_full[pj]);
      j = nbr_pj->nbr;
      //if (system->my_atoms[j].type < 0) continue;
      //orig_j  = system->my_atoms[j].orig_id;
      if (system->packed_atoms[j].type < 0) continue;
      orig_j  = system->packed_atoms[j].orig_id;

      flag = 1;

      //if (i >= nlocal && j >= nlocal) continue;
      //if (i >= nlocal)
      //{
      //  if (orig_j > orig_i) continue;
      //  if (orig_j == orig_i && urb(nbr_pj->dvec)) continue;
      //}
      //if (j >= nlocal)
      //{
      //  if (orig_i > orig_j) continue;
      //  if (orig_i == orig_j && llf(nbr_pj->dvec)) continue;
      //}
      

      r_ij = nbr_pj->d;
      //twbp = &(system->reax_param.tbp[ system->my_atoms[i].type ]
      //         [ system->my_atoms[j].type ]);
      twbp = &(system->reax_param.tbp[ system->packed_atoms[i].type ]
               [ system->packed_atoms[j].type ]);


      Tap = workspace->Tap[7] * r_ij + workspace->Tap[6];
      Tap = Tap * r_ij + workspace->Tap[5];
      Tap = Tap * r_ij + workspace->Tap[4];
      Tap = Tap * r_ij + workspace->Tap[3];
      Tap = Tap * r_ij + workspace->Tap[2];
      Tap = Tap * r_ij + workspace->Tap[1];
      Tap = Tap * r_ij + workspace->Tap[0];

      dTap = 7*workspace->Tap[7] * r_ij + 6*workspace->Tap[6];
      dTap = dTap * r_ij + 5*workspace->Tap[5];
      dTap = dTap * r_ij + 4*workspace->Tap[4];
      dTap = dTap * r_ij + 3*workspace->Tap[3];
      dTap = dTap * r_ij + 2*workspace->Tap[2];
      dTap += workspace->Tap[1]/r_ij;

      /*vdWaals Calculations*/
      if(system->reax_param.gp.vdw_type==1 || system->reax_param.gp.vdw_type==3)
      { // shielding
        powr_vdW1 = pow(r_ij, p_vdW1);
        powgi_vdW1 = pow( 1.0 / twbp->gamma_w, p_vdW1);

        fn13 = pow( powr_vdW1 + powgi_vdW1, p_vdW1i );
        exp1 = exp( twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );
        exp2 = exp( 0.5 * twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );

        e_vdW = twbp->D * (exp1 - 2.0 * exp2);
        data->my_en.e_vdW += Tap * e_vdW;

        dfn13 = pow( powr_vdW1 + powgi_vdW1, p_vdW1i - 1.0) *
          pow(r_ij, p_vdW1 - 2.0);

        CEvd = dTap * e_vdW -
          Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2) * dfn13;
      }
      else
      { // no shielding
        exp1 = exp( twbp->alpha * (1.0 - r_ij / twbp->r_vdW) );
        exp2 = exp( 0.5 * twbp->alpha * (1.0 - r_ij / twbp->r_vdW) );

        e_vdW = twbp->D * (exp1 - 2.0 * exp2);
        data->my_en.e_vdW += Tap * e_vdW;

        CEvd = dTap * e_vdW -
          Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2) / r_ij;
      }

      if(system->reax_param.gp.vdw_type==2 || system->reax_param.gp.vdw_type==3)
      { // inner wall
        e_core = twbp->ecore * exp(twbp->acore * (1.0-(r_ij/twbp->rcore)));
        data->my_en.e_vdW += Tap * e_core;

        de_core = -(twbp->acore/twbp->rcore) * e_core;
        CEvd += dTap * e_core + Tap * de_core / r_ij;

        //  lg correction, only if lgvdw is yes
        if (control->lgflag) {
          r_ij5 = pow( r_ij, 5.0 );
          r_ij6 = pow( r_ij, 6.0 );
          re6 = pow( twbp->lgre, 6.0 );
          e_lg = -(twbp->lgcij/( r_ij6 + re6 ));
          data->my_en.e_vdW += Tap * e_lg;

          de_lg = -6.0 * e_lg *  r_ij5 / ( r_ij6 + re6 ) ;
          CEvd += dTap * e_lg + Tap * de_lg / r_ij;
        }

      }

      /*Coulomb Calculations*/
      dr3gamij_1 = ( r_ij * r_ij * r_ij + twbp->gamma );
      dr3gamij_3 = pow( dr3gamij_1 , 0.33333333333333 );

      tmp = Tap / dr3gamij_3;
      data->my_en.e_ele += e_ele =
        C_ele * system->packed_atoms[i].q * system->packed_atoms[j].q * tmp;
        //C_ele * system->my_atoms[i].q * system->my_atoms[j].q * tmp;

      //CEclmb = C_ele * system->my_atoms[i].q * system->my_atoms[j].q *
      CEclmb = C_ele * system->packed_atoms[i].q * system->packed_atoms[j].q *
        ( dTap -  Tap * r_ij / dr3gamij_1 ) / dr3gamij_3;

      /* tally into per-atom energy */
      if( system->evflag || system->vflag_atom) {
        pe_vdw = Tap * (e_vdW + e_core + e_lg);
        //rvec_ScaledSum( delij, 1., system->my_atoms[i].x,
        //                -1., system->my_atoms[j].x );
        rvec_ScaledSum( delij, 1., system->packed_atoms[i].x,
                        -1., system->packed_atoms[j].x );

        f_tmp = -(CEvd + CEclmb);
        //system->ev_tally_full(i,pe_vdw,e_ele,f_tmp,delij[0],
        //                                delij[1],delij[2]);
        ev_tally_full(i, pe_vdw, e_ele, f_tmp, delij[0],delij[1],delij[2],
                      system->eflag_global, system->vflag_global, 
                      system->eflag_atom, system->vflag_atom,
                      &(system->eng_vdwl), &(system->eng_coul),
                      system->virial, system->eatom, system->vatom);
      }

      if( control->virial == 0 ) {
        //rvec_ScaledAdd( workspace->f[i], -(CEvd + CEclmb), nbr_pj->dvec );
        rvec_ScaledAdd( workspace->fCdDelta[i], -(CEvd + CEclmb), nbr_pj->dvec );
        // if (i == 0){
        //   double *f = workspace->f[i];
        //   double *d = nbr_pj->dvec;;
        //   double t[3];
        //   rvec_Copy(t, system->x[j]);
        //   rvec_ScaledAdd(t, -1, system->x[i]);
        //   printf("%d %f %f %f %f\n", j, -(CEvd + CEclmb), t[0], t[1], t[2]);
        //   printf("%d %f %f %f %f\n", j, -(CEvd + CEclmb), d[0], d[1], d[2]);
        //   rvec_Scale(t, -(CEvd +CEclmb), nbr_pj->dvec);
        //   printf("%d %f %f %f %f\n", j, -(CEvd + CEclmb), t[0], t[1], t[2]);
        //   printf("%d %f %f %f %f\n", j, -(CEvd + CEclmb), f[0], f[1], f[2]);
        // }
        // if (i == 0 || j == 0)
        //   printf("%5d %5d %f\n", i, j, -(CEvd + CEclmb));
        //rvec_ScaledAdd( workspace->f[j], +(CEvd + CEclmb), nbr_pj->dvec );
      }
      else { /* NPT, iNPT or sNPT */
        rvec_Scale( temp, CEvd + CEclmb, nbr_pj->dvec );

        //rvec_ScaledAdd( workspace->f[i], -1., temp );
        rvec_ScaledAdd( workspace->fCdDelta[i], -1., temp );
        //rvec_Add( workspace->f[j], temp );

        rvec_iMultiply( ext_press, nbr_pj->rel_box, temp );
        //rvec_Add( data->my_ext_press, ext_press );
        // if (i == 0 || j == 0)
        //   printf("%d %d %f %f %f\n", i, j, ext_press[0], ext_press[1], ext_press[2]);
        rvec_ScaledAdd( data->my_ext_press, 0.5, ext_press );
      }
    }//for-pj
    system->virial[0] += system->packed_atoms[i].x[0] * workspace->fCdDelta[i][0];
    system->virial[1] += system->packed_atoms[i].x[1] * workspace->fCdDelta[i][1];
    system->virial[2] += system->packed_atoms[i].x[2] * workspace->fCdDelta[i][2];
    system->virial[3] += system->packed_atoms[i].x[0] * workspace->fCdDelta[i][1];
    system->virial[4] += system->packed_atoms[i].x[0] * workspace->fCdDelta[i][2];
    system->virial[5] += system->packed_atoms[i].x[1] * workspace->fCdDelta[i][2];

  }    

    //Compute_Polarization_Energy( system, data );
}

void vdW_Coulomb_Energy_Full_C(vdw_coulomb_pack_t * param)
{
  if(athread_idle() == 0)
    athread_init();
  
  //perf_config_t conf;
  //conf.pcrc = PCRC_ALL;
  //conf.pcr0 = PC0_CYCLE;
  //conf.pcr1 = PC1_CYCLE;
  ////conf.pcr2 = PC2_N_GLD;
  //conf.pcr2 = PC2_N_DMA_REQ;
  //lwpf_init(&conf);
  //int myrank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int tmp = param->system->vflag_global;
  param->system->vflag_global = 1;


  //2019/11/10
  reax_system_c *system       = param->system;
  control_params *control     = param->control;
  simulation_data *data       = param->data;
  storage *workspace          = param->workspace;
  reax_list **lists           = param->lists;

  reax_list *far_nbrs         = (*lists) + FAR_NBRS_FULL;
  param->index                = far_nbrs->index;
  param->end_index            = far_nbrs->end_index;
  param->far_nbr_list_full    = far_nbrs->select.far_nbr_list_full;


  param->n                    = system->n;
  param->N                    = system->N;
  param->evflag               = system->evflag;
  param->eflag_global         = system->eflag_global;
  param->eflag_atom           = system->eflag_atom;
  param->vflag_global         = system->vflag_global;
  param->vflag_atom           = system->vflag_atom;
  param->vdw_type             = system->reax_param.gp.vdw_type; 
  param->num_atom_types       = system->reax_param.num_atom_types;
  param->l28                  = system->reax_param.gp.l[28];
  
  param->wksp_tap0 = workspace->Tap[0];
  param->wksp_tap1 = workspace->Tap[1];
  param->wksp_tap2 = workspace->Tap[2];
  param->wksp_tap3 = workspace->Tap[3];
  param->wksp_tap4 = workspace->Tap[4];
  param->wksp_tap5 = workspace->Tap[5];
  param->wksp_tap6 = workspace->Tap[6];
  param->wksp_tap7 = workspace->Tap[7];

  int ntypes = system->reax_param.num_atom_types;
  tbp_boost_t tbpb[ntypes][ntypes];
  two_body_parameters **tbpr = system->reax_param.tbp;
  vdw_coulomb_tbp pack_tbpr[ntypes][ntypes];
  double packed_eng[16]={0};
  
  double p_vdW1 = system->reax_param.gp.l[28];
  //printf("p_vdW1 = %lf, %lf\n", p_vdW1, 1/p_vdW1);

  
  int i, j;
    
  //for (i = 0; i < ntypes; i++)//ntypes=4;
  //{
  //  for (j = 0; j < ntypes; j++)
  //  {
  //    tbpb[i][j].powgi_vdW1     = pow(tbpr[i][j].gamma_w, -p_vdW1);
  //    tbpb[i][j].r_vdWinv       = 1.0 / tbpr[i][j].r_vdW;
  //    tbpb[i][j].alpha_div_rvdW = tbpr[i][j].alpha / tbpr[i][j].r_vdW;
  //    if(system->reax_param.gp.vdw_type==2 || system->reax_param.gp.vdw_type==3)
  //    { 
  //      tbpb[i][j].a_div_rcore = tbpr[i][j].acore / tbpr[i][j].rcore;
  //      tbpb[i][j].rcoreinv = 1.0 / tbpr[i][j].rcore;
  //    }
  //    double lgresq = tbpr[i][j].lgre * tbpr[i][j].lgre;
  //    tbpb[i][j].re6 = lgresq * lgresq * lgresq;

  //    pack_tbpr[i][j].alpha = tbpr[i][j].alpha;
  //    pack_tbpr[i][j].ecore = tbpr[i][j].ecore;
  //    pack_tbpr[i][j].acore = tbpr[i][j].acore;
  //    pack_tbpr[i][j].gamma = tbpr[i][j].gamma;
  //    pack_tbpr[i][j].lgcij = tbpr[i][j].lgcij;
  //    pack_tbpr[i][j].D     = tbpr[i][j].D;
  //  }
  //}

  

  //if(control->virial == 0 && (param->vdw_type == 1) && system->vflag_atom==0 && system->eflag_atom==0 && system->vflag_global==1 && system->eflag_global==1 && system->evflag==1)
  if(control->virial == 0 && (param->vdw_type == 1) && system->vflag_atom==0 && system->eflag_atom==0 && system->vflag_global==1 )
  {

    //double powgi_vdW1[ntypes][ntypes];
    //double r_vdWinv[ntypes][ntypes];
    //double alpha_div_rvdW[ntypes][ntypes];
    //double a_div_rcore[ntypes][ntypes];
    //double rcoreinv[ntypes][ntypes];
    //double re6[ntypes][ntypes];

    //double alpha[ntypes][ntypes];
    //double ecore[ntypes][ntypes];
    //double acore[ntypes][ntypes];
    //double gamma[ntypes][ntypes];
    //double lgcij[ntypes][ntypes];
    //double D[ntypes][ntypes];
    double params[ntypes*12][ntypes];

    param->powgi_vdW1     = &params[0][0];
    param->r_vdWinv       = &params[ntypes][0];
    param->alpha_div_rvdW = &params[2 * ntypes][0];
    param->alpha          = &params[3 * ntypes][0];
    param->gamma          = &params[4 * ntypes][0];
    param->D              = &params[5 * ntypes][0];
    param->a_div_rcore    = &params[6 * ntypes][0];
    param->rcoreinv       = &params[7 * ntypes][0];
    param->re6            = &params[8 * ntypes][0];
    param->ecore          = &params[9 * ntypes][0];
    param->acore          = &params[10 * ntypes][0];
    param->lgcij          = &params[11 * ntypes][0];


    for (i = 0; i < ntypes; i++)//ntypes=4;
    {
      for (j = 0; j < ntypes; j++)
      {
        param->powgi_vdW1[i*ntypes+j] = pow(tbpr[i][j].gamma_w, -p_vdW1);
        param->r_vdWinv[i*ntypes+j] = 1.0 / tbpr[i][j].r_vdW;
        param->alpha_div_rvdW[i*ntypes+j] = tbpr[i][j].alpha / tbpr[i][j].r_vdW;
        if(system->reax_param.gp.vdw_type==2 || system->reax_param.gp.vdw_type==3)
        { 
          param->a_div_rcore[i*ntypes+j] = tbpr[i][j].acore / tbpr[i][j].rcore;
          param->rcoreinv[i*ntypes+j] = 1.0 / tbpr[i][j].rcore;
        }
        double lgresq = tbpr[i][j].lgre * tbpr[i][j].lgre;
        param->re6[i*ntypes+j] = lgresq * lgresq * lgresq;

        param->alpha[i*ntypes+j] = tbpr[i][j].alpha;
        param->ecore[i*ntypes+j] = tbpr[i][j].ecore;
        param->acore[i*ntypes+j] = tbpr[i][j].acore;
        param->gamma[i*ntypes+j] = tbpr[i][j].gamma;
        param->lgcij[i*ntypes+j] = tbpr[i][j].lgcij;
        param->D[i*ntypes+j] = tbpr[i][j].D;
      }
    }


    //param->powgi_vdW1     = powgi_vdW1    ;
    //param->r_vdWinv       = r_vdWinv      ;
    //param->alpha_div_rvdW = alpha_div_rvdW;
    //param->a_div_rcore    = a_div_rcore   ;
    //param->rcoreinv       = rcoreinv      ;
    //param->re6            = re6           ;
    //param->alpha  = alpha;
    //param->ecore  = ecore;
    //param->acore  = acore;
    //param->gamma  = gamma;
    //param->lgcij  = lgcij;
    //param->D      = D    ;


    param->tbpb                 = tbpb;
    param->pack_tbpr            = pack_tbpr;
    param->packed_eng           = packed_eng;
    param->packed_atoms         = system->packed_atoms;
    param->f                    = workspace->fCdDelta;
    
    //printf("special\n");
    //athread_spawn(vdW_Coulomb_Energy_Full_C_para_special, param);
    if(system->evflag)
      athread_spawn(vdW_Coulomb_Energy_cpe_evflag1, param);
    else
      athread_spawn(vdW_Coulomb_Energy_cpe_evflag0, param);

    
    Validate_Lists_C(param);
    sort_bonds_C(param);
    athread_join();

    system->eng_vdwl      += packed_eng[0];
    system->eng_coul      += packed_eng[1];
    system->virial[0]     += packed_eng[2];
    system->virial[1]     += packed_eng[3];
    system->virial[2]     += packed_eng[4];
    system->virial[3]     += packed_eng[5];
    system->virial[4]     += packed_eng[6];
    system->virial[5]     += packed_eng[7];
    data->my_ext_press[0] += packed_eng[8];
    data->my_ext_press[1] += packed_eng[9];
    data->my_ext_press[2] += packed_eng[10];
    data->my_en.e_vdW     += packed_eng[11];
    data->my_en.e_ele     += packed_eng[12];
  }
  else
  {
    //printf("vdw coulomb general\n");
    for (i = 0; i < ntypes; i++)//ntypes=4;
    {
      for (j = 0; j < ntypes; j++)
      {
        tbpb[i][j].powgi_vdW1     = pow(tbpr[i][j].gamma_w, -p_vdW1);
        tbpb[i][j].r_vdWinv       = 1.0 / tbpr[i][j].r_vdW;
        tbpb[i][j].alpha_div_rvdW = tbpr[i][j].alpha / tbpr[i][j].r_vdW;
        if(system->reax_param.gp.vdw_type==2 || system->reax_param.gp.vdw_type==3)
        { 
          tbpb[i][j].a_div_rcore = tbpr[i][j].acore / tbpr[i][j].rcore;
          tbpb[i][j].rcoreinv = 1.0 / tbpr[i][j].rcore;
        }
        double lgresq = tbpr[i][j].lgre * tbpr[i][j].lgre;
        tbpb[i][j].re6 = lgresq * lgresq * lgresq;

        pack_tbpr[i][j].alpha = tbpr[i][j].alpha;
        pack_tbpr[i][j].ecore = tbpr[i][j].ecore;
        pack_tbpr[i][j].acore = tbpr[i][j].acore;
        pack_tbpr[i][j].gamma = tbpr[i][j].gamma;
        pack_tbpr[i][j].lgcij = tbpr[i][j].lgcij;
        pack_tbpr[i][j].D     = tbpr[i][j].D;
      }
    }
    
    param->tbpb                 = tbpb;
    param->pack_tbpr            = pack_tbpr;
    param->packed_eng           = packed_eng;
    param->packed_atoms         = system->packed_atoms;
    param->f                    = workspace->fCdDelta;

    athread_spawn(vdW_Coulomb_Energy_Full_C_para_general, param);
  
    Validate_Lists_C(param);
    sort_bonds_C(param);

    athread_join();
  }

  //if(myrank == 0)
  //{
  //  lwpf_report_summary(stdout, &conf);
  //}

  param->system->vflag_global = tmp;

  //if(r == 10 && rank == 0)
  //{
  //  lwpf_finish(stdout);
  //}
  //r++;
}
#endif

#ifdef CPE
#include "STUBS/mpi.h"
#include <dma.h>
#include "dma_macros.h"
#define ISTEP 64
#define SNSTEP 32
#define NTYP 4
//#define LWPF_UNIT U(NONBOND)
//#define LWPF_KERNELS K(ALL) K(BEFORE)  K(CMP) K(SEC1) K(TRANS) K(VMAD1) K(MATH) K(POW) K(POW2) K(DIV) K(VMAD) K(REMAINDER) K(SUM)  K(REDUCE)  
//#include "lwpf2.h"
#include "poly_math.h"
//#define DMA_FAST
void vdW_Coulomb_Energy_Full_C_para_general(vdw_coulomb_pack_t * param)
{
  //lwpf_start(ALL);
  //lwpf_start(BEFORE);
  dma_init();
  vdw_coulomb_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(vdw_coulomb_pack_t));
  dma_syn();
  reax_system_c *system       = l_pm.system;
  control_params *control     = l_pm.control;
  simulation_data *data       = l_pm.data;
  storage *workspace          = l_pm.workspace;
  reax_list **lists           = l_pm.lists;
  
  reax_system_c   l_sys;
  control_params  l_ctrl;
  pe_get(system,   &l_sys,  sizeof(reax_system_c));
  pe_get(control,  &l_ctrl, sizeof(control_params));
  dma_syn();

  int i, j, pj, natoms, nlocal;
  int start_i, end_i, flag;
  rc_tagint orig_i, orig_j;
  double p_vdW1, p_vdW1i;
  double powr_vdW1, powgi_vdW1;
  double tmp, r_ij, fn13, exp1, exp2;
  double Tap, dTap, dfn13, CEvd, CEclmb, de_core;
  double dr3gamij_1, dr3gamij_3;
  double e_ele, e_vdW, e_core, SMALL = 0.0001;
  double e_lg, de_lg, r_ij5, r_ij6, re6;
  rvec temp, ext_press;
  two_body_parameters *twbp;
  far_neighbor_data_full *nbr_pj;
  reax_list *far_nbrs;
  reax_list l_fnbrs;

  // Tallying variables:
  double pe_vdw, f_tmp, delij[3];
  natoms = l_sys.N;
  nlocal = l_sys.n;
  far_nbrs = (*lists) + FAR_NBRS_FULL;
  p_vdW1 = l_sys.reax_param.gp.l[28];
  p_vdW1i = 1.0 / p_vdW1;
  e_core = 0;
  e_vdW = 0;
  e_lg = de_lg = 0.0;

  int ntypes = l_sys.reax_param.num_atom_types;
  tbp_boost_t tbpb[ntypes][ntypes];
  tbp_boost_t *tbpt;
  two_body_parameters **tbpr = l_sys.reax_param.tbp;
  two_body_parameters l_tbpr[ntypes][ntypes];
  int ist, isz, ied;
  int ii, ioff;
  int idx_st[ISTEP], idx_ed[ISTEP];
  double ei[ISTEP], vi[ISTEP][6];
  double fi[ISTEP][4];
  atom_pack_t packed_atoms[ISTEP];
  far_neighbor_data_full pj_nbr[SNSTEP];

  double wksp_tap[8];
  double packed_eng[16];
  pe_get(&(system->eng_vdwl), packed_eng, sizeof(double)*8);
  pe_get(&(data->my_ext_press), packed_eng+8, sizeof(double)*3);
  pe_get(&(data->my_en.e_vdW), packed_eng+11, sizeof(double)*2);
  pe_get(workspace->Tap, wksp_tap, sizeof(double)*8);
  pe_get(tbpr[0], l_tbpr, sizeof(two_body_parameters)*ntypes * ntypes);
  pe_get(far_nbrs, &l_fnbrs, sizeof(reax_list));
  dma_syn();


  doublev4 eng_virial[4];
  for(i = 0; i < 4; i++)
    eng_virial[i] = 0;
  double *eng_vdwl = (double*)(void*)eng_virial;
  double *eng_coul = eng_vdwl + 1; 
  double *virial   = eng_coul + 1;
  double *extp     = virial + 6;
  double *data_vdW = extp + 3;
  double *data_ele = data_vdW + 1;
  for (i = 0; i < ntypes; i ++)//ntypes=4;
  {
    for (j = 0; j < ntypes; j ++)
    {
      tbpb[i][j].powgi_vdW1     = p_powd(l_tbpr[i][j].gamma_w, -p_vdW1);
      tbpb[i][j].r_vdWinv       = 1.0 / l_tbpr[i][j].r_vdW;
      tbpb[i][j].alpha_div_rvdW = l_tbpr[i][j].alpha / l_tbpr[i][j].r_vdW;
      if(l_sys.reax_param.gp.vdw_type==2 || l_sys.reax_param.gp.vdw_type==3)
      { 
        tbpb[i][j].a_div_rcore = l_tbpr[i][j].acore / l_tbpr[i][j].rcore;
        tbpb[i][j].rcoreinv = 1.0 / l_tbpr[i][j].rcore;
      }
      double lgresq = l_tbpr[i][j].lgre * l_tbpr[i][j].lgre;
      tbpb[i][j].re6 = lgresq * lgresq * lgresq;
    }
  }
  volatile dma_desc get_desc = 0;
  volatile get_reply = 0;
  dma_set_mode(&get_desc, PE_MODE);
  dma_set_op(&get_desc, DMA_GET);
  dma_set_reply(&get_desc, &get_reply);
  //lwpf_stop(BEFORE);
  //for(ist = _MYID * ISTEP; ist < natoms; ist += ISTEP * 64)
  for(ist = _MYID * ISTEP; ist < nlocal; ist += ISTEP * 64)
  {
    ied = ist + ISTEP;
    if(ied > nlocal)
      ied = nlocal;
    isz = ied - ist;

    pe_get(l_fnbrs.index+ist, idx_st, sizeof(int)*isz);
    pe_get(l_fnbrs.end_index+ist, idx_ed, sizeof(int)*isz);
    pe_get(system->packed_atoms+ist, packed_atoms, sizeof(atom_pack_t)*isz);
    pe_get(workspace->fCdDelta[ist], fi[0], sizeof(double)*isz*4);
    if(l_sys.vflag_atom) 
      pe_get(system->vatom[ist], vi[0], sizeof(double)*isz*6);
    if(l_sys.eflag_atom) 
      pe_get(system->eatom+ist, ei, sizeof(double)*isz);
    dma_syn();

    for( ii = ist; ii < ied; ++ii) 
    {
      i = ii;
      int ioff = i - ist;
      if (packed_atoms[ioff].type < 0) continue;
      start_i = idx_st[ioff];
      end_i   = idx_ed[ioff];
      orig_i  = packed_atoms[ioff].orig_id;
      int jst, jed, jsz;
      for(jst = start_i; jst < end_i; jst+=SNSTEP)
      {
        jsz = SNSTEP;
        if(jst + SNSTEP > end_i)
          jsz = end_i - jst;

        get_reply = 0;
        dma_set_size(&get_desc, sizeof(far_neighbor_data_full)*jsz);
        dma(get_desc,l_fnbrs.select.far_nbr_list_full+jst, pj_nbr);
        while(get_reply != 1);
        
        //lwpf_start(PJ);
        for( pj = 0; pj < jsz; pj++ ) 
        {
          nbr_pj = &(pj_nbr[pj]);
          j = nbr_pj->nbr;
          if (nbr_pj->type < 0) continue;
          orig_j  = nbr_pj->orig_id;
          flag = 1;
          //if (i >= nlocal && j >= nlocal) continue;
          //if (i >= nlocal)
          //{
          //  if (orig_j > orig_i) continue;
          //  if (orig_j == orig_i && urb(nbr_pj->dvec)) continue;
          //}
          //if (j >= nlocal)
          //{
          //  if (orig_i > orig_j) continue;
          //  if (orig_i == orig_j && llf(nbr_pj->dvec)) continue;
          //}

          //lwpf_start(CEVD);
          r_ij = nbr_pj->d;
          double rijinv = 1 / r_ij;
          double r2ij = r_ij * r_ij;
          twbp = &(l_tbpr[packed_atoms[ioff].type][nbr_pj->type]);
          tbpt = tbpb[packed_atoms[ioff].type] + nbr_pj->type;
                  
          Tap = wksp_tap[7] * r_ij + wksp_tap[6];
          Tap = Tap * r_ij + wksp_tap[5];
          Tap = Tap * r_ij + wksp_tap[4];
          Tap = Tap * r_ij + wksp_tap[3];
          Tap = Tap * r_ij + wksp_tap[2];
          Tap = Tap * r_ij + wksp_tap[1];
          Tap = Tap * r_ij + wksp_tap[0];

          dTap = 7*wksp_tap[7] * r_ij + 6*wksp_tap[6];
          dTap = dTap * r_ij + 5*wksp_tap[5];
          dTap = dTap * r_ij + 4*wksp_tap[4];
          dTap = dTap * r_ij + 3*wksp_tap[3];
          dTap = dTap * r_ij + 2*wksp_tap[2];
          dTap += wksp_tap[1] * rijinv;
          //lwpf_start(CEIFS);
          if(l_sys.reax_param.gp.vdw_type==1 || l_sys.reax_param.gp.vdw_type==3)
          { 
            powr_vdW1 = p_powd(r_ij, p_vdW1);
            powgi_vdW1 = tbpt->powgi_vdW1;
            fn13 = p_powd( powr_vdW1 + powgi_vdW1, p_vdW1i );
            exp2 = p_expd( 0.5 * twbp->alpha * (1.0 - fn13 * tbpt->r_vdWinv) );
            exp1 = exp2 * exp2;
            e_vdW = twbp->D * (exp1 - 2.0 * exp2);
            *data_vdW += Tap * e_vdW;

            dfn13 = p_powd( powr_vdW1 + powgi_vdW1, p_vdW1i-1.0) 
                    * p_powd(r_ij, p_vdW1 - 2.0);
            CEvd=dTap*e_vdW -Tap*twbp->D * (tbpt->alpha_div_rvdW)*(exp1-exp2)*dfn13;
          }
          else
          { 
            exp2 = p_expd( 0.5 * twbp->alpha * (1.0 - r_ij * tbpt->r_vdWinv) );
            exp1 = exp2 * exp2;
            e_vdW = twbp->D * (exp1 - 2.0 * exp2);
            *data_vdW += Tap * e_vdW;
            CEvd = dTap*e_vdW -Tap*twbp->D*(tbpt->alpha_div_rvdW)*(exp1-exp2)*rijinv;
          }

          //lwpf_start(CEIFS2);
          if(l_sys.reax_param.gp.vdw_type==2 || l_sys.reax_param.gp.vdw_type==3)
          { 
            e_core = twbp->ecore*p_expd(twbp->acore*(1.0-(r_ij*tbpt->rcoreinv)));
            *data_vdW += Tap * e_core;
            de_core = -(tbpt->a_div_rcore) * e_core;
            CEvd += dTap * e_core + Tap * de_core  * rijinv;
            if (l_ctrl.lgflag) 
            {
              r_ij5 = r2ij * r2ij * r_ij;
              r_ij6 = r2ij * r2ij * r2ij;
              re6 = tbpt->re6;
              double rrinv = 1 / (r_ij6 + re6);
              e_lg = -(twbp->lgcij * rrinv);
              *data_vdW += Tap * e_lg;
              de_lg = -6.0 * e_lg *  r_ij5 * rrinv;
              CEvd += dTap * e_lg + Tap * de_lg  * rijinv;
            }
          }
          //lwpf_stop(CEIFS2);
          //lwpf_stop(CEIFS);
          dr3gamij_1 = ( r_ij * r_ij * r_ij + twbp->gamma );
          //lwpf_start(DIVF);
          dr3gamij_3 = p_powd( dr3gamij_1 , 0.33333333333333 );

          //lwpf_stop(DIVF);
          double dr3gamij_1inv = 1 / dr3gamij_1;
          double dr3gamij_3inv = 1 / dr3gamij_3;
          tmp = Tap * dr3gamij_3inv;

          *data_ele += e_ele = C_ele * packed_atoms[ioff].q * nbr_pj->q * tmp;
          CEclmb = C_ele * packed_atoms[ioff].q * nbr_pj->q *
                   (dTap -  Tap * r_ij * dr3gamij_1inv ) * dr3gamij_3inv;
          
          //lwpf_stop(CEVD);

          //lwpf_start(EV);
          if( l_sys.evflag || l_sys.vflag_atom || l_sys.vflag_global) 
          {
            pe_vdw = Tap * (e_vdW + e_core + e_lg);
            f_tmp = -(CEvd + CEclmb);
            double v[6];
            double evdwl = pe_vdw, ecoul = e_ele, fpair = f_tmp;
            double *del = nbr_pj->dvec;
            if (l_sys.eflag_global) 
            { 
              *eng_vdwl += 0.5*evdwl;
              *eng_coul += 0.5*ecoul;
            }
            if (l_sys.eflag_atom) 
              ei[ioff] += 0.5 * (evdwl + ecoul);
            if (l_sys.vflag_global || l_sys.vflag_atom) 
            {
              v[0] = 0.5*del[0]*del[0]*fpair;
              v[1] = 0.5*del[1]*del[1]*fpair;
              v[2] = 0.5*del[2]*del[2]*fpair;
              v[3] = 0.5*del[0]*del[1]*fpair;
              v[4] = 0.5*del[0]*del[2]*fpair;
              v[5] = 0.5*del[1]*del[2]*fpair;
              if (l_sys.vflag_global) 
              { 
                virial[0] += v[0];
                virial[1] += v[1];
                virial[2] += v[2];
                virial[3] += v[3];
                virial[4] += v[4];
                virial[5] += v[5];
              }
              if (l_sys.vflag_atom) 
              { 
                vi[ioff][0] += v[0];
                vi[ioff][1] += v[1];
                vi[ioff][2] += v[2];
                vi[ioff][3] += v[3];
                vi[ioff][4] += v[4];
                vi[ioff][5] += v[5];
              }
            }//if
          }//if

          double tmp_c = CEvd + CEclmb;
          if(l_ctrl.virial == 0 ) 
          {
            fi[ioff][0] += -tmp_c * nbr_pj->dvec[0];
            fi[ioff][1] += -tmp_c * nbr_pj->dvec[1];
            fi[ioff][2] += -tmp_c * nbr_pj->dvec[2];
          }
          else 
          { 
            temp[0] = (CEvd + CEclmb) *  nbr_pj->dvec[0];
            temp[1] = (CEvd + CEclmb) *  nbr_pj->dvec[1];
            temp[2] = (CEvd + CEclmb) *  nbr_pj->dvec[2];
            fi[ioff][0] += -temp[0];
            fi[ioff][1] += -temp[1];
            fi[ioff][2] += -temp[2];

            ext_press[0] = nbr_pj->rel_box[0] * temp[0];
            ext_press[1] = nbr_pj->rel_box[1] * temp[1];
            ext_press[2] = nbr_pj->rel_box[2] * temp[2];
            extp[0] += 0.5 * ext_press[0];
            extp[1] += 0.5 * ext_press[1];
            extp[2] += 0.5 * ext_press[2];
          }//else
          //lwpf_stop(EV);
        }//for-pj
        //lwpf_stop(PJ);
      }//for-jst 
      virial[0] += packed_atoms[ioff].x[0] * fi[ioff][0];
      virial[1] += packed_atoms[ioff].x[1] * fi[ioff][1];
      virial[2] += packed_atoms[ioff].x[2] * fi[ioff][2];
      virial[3] += packed_atoms[ioff].x[0] * fi[ioff][1];
      virial[4] += packed_atoms[ioff].x[0] * fi[ioff][2];
      virial[5] += packed_atoms[ioff].x[1] * fi[ioff][2];

    }//for-ii
    pe_put(workspace->fCdDelta[ist], fi[0], sizeof(double)*isz*4);
    if (l_sys.eflag_atom) 
      pe_put(system->eatom+ist, ei, sizeof(double)*isz);
    if (l_sys.vflag_atom) 
      pe_put(system->vatom[ist], vi[0], sizeof(double)*isz*6);
    dma_syn();
  }//for-ist
  reg_reduce_inplace_doublev4(eng_virial, 4);

  if(_MYID == 0)
  { 
    *data_vdW += packed_eng[11];
    *data_ele += packed_eng[12];
    pe_put(&(data->my_en.e_vdW), data_vdW, sizeof(double)*2);
    dma_syn();
    if(l_sys.eflag_global) 
    {
      *eng_vdwl += packed_eng[0];
      *eng_coul += packed_eng[1];
      pe_put(&(system->eng_vdwl), eng_vdwl, sizeof(double)*2);
      dma_syn();
    }
    if(l_sys.vflag_global) 
    {
      virial[0] += packed_eng[2];  
      virial[1] += packed_eng[3]; 
      virial[2] += packed_eng[4];  
      virial[3] += packed_eng[5];  
      virial[4] += packed_eng[6];  
      virial[5] += packed_eng[7];  
      pe_put(&(system->virial[0]), virial, sizeof(double)*6);
      dma_syn();
    }
    if(l_ctrl.virial != 0)
    { 
      extp[0] += packed_eng[8] ;
      extp[1] += packed_eng[9] ;
      extp[2] += packed_eng[10];
      pe_put(&(data->my_ext_press[0]), extp, sizeof(double)*3);
      dma_syn();
    }
  }
  //lwpf_stop(ALL);
}

#define vshuffd_rc(a, b, c, d) (d | (c << 2) | (b << 4) | (a << 6))
#define simd_vsumd(x){                                  \
  x += simd_vshff(x, x, vshuffd_rc(2, 3, 0, 1));        \
  x += simd_vshff(x, x, vshuffd_rc(1, 0, 3, 2));        \
}
#define transpose8x4( in0, in1, in2, in3, in4, in5, in6, in7,\
                      ot0, ot1, ot2, ot3, ot4, ot5, ot6, ot7)\
{\
  doublev4 o0 = simd_vshff(in1,in0,0x44);\
  doublev4 o1 = simd_vshff(in1,in0,0xEE);\
  doublev4 o2 = simd_vshff(in3,in2,0x44);\
  doublev4 o3 = simd_vshff(in3,in2,0xEE);\
  ot0 = simd_vshff(o2,o0,0x88); \
  ot1 = simd_vshff(o2,o0,0xDD); \
  ot2 = simd_vshff(o3,o1,0x88); \
  ot3 = simd_vshff(o3,o1,0xDD); \
  doublev4 o4 = simd_vshff(in5,in4,0x44); \
  doublev4 o5 = simd_vshff(in5,in4,0xEE); \
  doublev4 o6 = simd_vshff(in7,in6,0x44); \
  doublev4 o7 = simd_vshff(in7,in6,0xEE); \
  ot4 = simd_vshff(o6,o4,0x88); \
  ot5 = simd_vshff(o6,o4,0xDD); \
  ot6 = simd_vshff(o7,o5,0x88); \
  ot7 = simd_vshff(o7,o5,0xDD); \
}

#define transpose8x4( in0, in1, in2, in3, in4, in5, in6, in7,\
                      ot0, ot1, ot2, ot3, ot4, ot5, ot6, ot7)\
{\
  doublev4 o0 = simd_vshff(in1,in0,0x44);\
  doublev4 o1 = simd_vshff(in1,in0,0xEE);\
  doublev4 o2 = simd_vshff(in3,in2,0x44);\
  doublev4 o3 = simd_vshff(in3,in2,0xEE);\
  ot0 = simd_vshff(o2,o0,0x88); \
  ot3 = simd_vshff(o3,o1,0xDD); \
  doublev4 o4 = simd_vshff(in5,in4,0x44); \
  doublev4 o5 = simd_vshff(in5,in4,0xEE); \
  doublev4 o6 = simd_vshff(in7,in6,0x44); \
  doublev4 o7 = simd_vshff(in7,in6,0xEE); \
  ot4 = simd_vshff(o6,o4,0x88); \
  ot5 = simd_vshff(o6,o4,0xDD); \
  ot6 = simd_vshff(o7,o5,0x88); \
  ot7 = simd_vshff(o7,o5,0xDD); \
}

void vdW_Coulomb_Energy_Full_C_para_special(vdw_coulomb_pack_t * param)
{
  //lwpf_enter(NONBOND);
  //lwpf_start(ALL);
  //lwpf_start(BEFORE);
  dma_init();
  vdw_coulomb_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(vdw_coulomb_pack_t));
  dma_syn();
  reax_system_c *system       = l_pm.system;
  control_params *control     = l_pm.control;
  simulation_data *data       = l_pm.data;
  storage *workspace          = l_pm.workspace;
  reax_list **lists           = l_pm.lists;
  
  int i;
  //int j;
  int pj, natoms, nlocal;
  int start_i, end_i;//, flag;
  rc_tagint orig_i, orig_j;
  double p_vdW1;//, p_vdW1i;
  double powr_vdW1, powgi_vdW1;
  double tmp, r_ij, fn13, exp1, exp2;
  double Tap, dTap, dfn13, CEvd, CEclmb, de_core;
  double dr3gamij_1, dr3gamij_3;
  double e_ele, e_vdW, e_core, SMALL = 0.0001;
  double e_lg, de_lg, r_ij5, r_ij6, re6;
  rvec temp, ext_press;
  far_neighbor_data_full *nbr_pj;
  reax_list *far_nbrs;

  // Tallying variables:
  double pe_vdw, f_tmp, delij[3];
  natoms = l_pm.N;
  nlocal = l_pm.n;
  far_nbrs = (*lists) + FAR_NBRS_FULL;
  p_vdW1 = l_pm.l28;
  //p_vdW1i = 1.0 / p_vdW1;
  e_core = 0;
  e_vdW = 0;
  e_lg = de_lg = 0.0;

  int ntypes = l_pm.num_atom_types;
  tbp_boost_t tbpb[ntypes][ntypes];
  tbp_boost_t *tbpt;
  vdw_coulomb_tbp l_tbpr[ntypes][ntypes];
  vdw_coulomb_tbp *twbp;

  int ist, isz, ied;
  int ii, ioff;
  int idx_st[ISTEP], idx_ed[ISTEP];
  double ei[ISTEP], vi[ISTEP][6];
  double fi[ISTEP][4];
  atom_pack_t packed_atoms[ISTEP];
  //far_neighbor_data_full pj_nbr[SNSTEP];
  doublev4  nbr_v4[SNSTEP*2], *p_nbr_v4;
  far_neighbor_data_full *pj_nbr = nbr_v4;
  
  doublev4 t_powgi_vdW1_v4[NTYP];
  doublev4 t_r_vdWinv_v4[NTYP];
  doublev4 t_alpha_div_rvdW_v4[NTYP];
  doublev4 t_alpha_v4[NTYP];
  doublev4 t_gamma_v4[NTYP];
  doublev4 t_D_v4[NTYP];
  
  double *p_powgi_vdW1      = t_powgi_vdW1_v4;
  double *p_r_vdWinv        = t_r_vdWinv_v4;
  double *p_alpha_div_rvdW  = t_alpha_div_rvdW_v4;
  double *p_alpha           = t_alpha_v4;
  double *p_gamma           = t_gamma_v4;
  double *p_D               = t_D_v4;

  double *params = l_pm.params;

  //pe_get(l_pm.tbpb, tbpb, sizeof(tbp_boost_t) * ntypes * ntypes);
  //pe_get(l_pm.pack_tbpr, l_tbpr, sizeof(vdw_coulomb_tbp) * ntypes * ntypes);
  pe_get(l_pm.powgi_vdW1, p_powgi_vdW1, sizeof(double) * ntypes * ntypes);
  pe_get(l_pm.r_vdWinv      , p_r_vdWinv      , sizeof(double) * ntypes * ntypes);
  pe_get(l_pm.alpha_div_rvdW, p_alpha_div_rvdW, sizeof(double) * ntypes * ntypes);
  pe_get(l_pm.alpha         , p_alpha         , sizeof(double) * ntypes * ntypes);
  pe_get(l_pm.gamma         , p_gamma         , sizeof(double) * ntypes * ntypes);
  pe_get(l_pm.D             , p_D             , sizeof(double) * ntypes * ntypes);
  
  dma_syn();

  doublev4 eng_virial[4];
  for(i = 0; i < 4; i++)
    eng_virial[i] = 0;
  double *eng_vdwl = (double*)(void*)eng_virial;
  double *eng_coul = eng_vdwl + 1; 
  double *virial   = eng_coul + 1;
  double *extp     = virial + 6;
  double *data_vdW = extp + 3;
  double *data_ele = data_vdW + 1;
  
  volatile dma_desc get_desc = 0;
  volatile get_reply = 0;
  dma_set_mode(&get_desc, PE_MODE);
  dma_set_op(&get_desc, DMA_GET);
  dma_set_reply(&get_desc, &get_reply);

  //vectorize
  doublev4 r_ij_v4, rijinv_v4;
  doublev4 fn13_v4, dfn13_v4;
  doublev4 exp1_v4, exp2_v4;
  doublev4 evdW_v4;

  doublev4 dr3gamij_1_v4, dr3gamij_3_v4;
  doublev4 dr3gamij_1inv_v4, dr3gamij_3inv_v4;
  doublev4 CEvd_v4, CEclmb_v4;
  doublev4 tmp_v4, tmpc_v4, fpair_v4;
  doublev4 evdwl_v4, ecoul_v4;


  doublev4 alpha_v4;
  doublev4 D_v4, gamma_v4;
  doublev4 powgi_vdW1_v4;
  doublev4 r_vdWinv_v4;
  doublev4 alpha_div_rvdW_v4;
  doublev4 del0_v4, del1_v4, del2_v4;
  doublev4 qi_v4, qj_v4;
  doublev4 half_v4, one_v4, two_v4, thr_v4, four_v4;
  doublev4 five_v4, six_v4, seven_v4, thrinv_v4;

  //half_v4 = 0.5; one_v4 = 1.0; two_v4 = 2.0; thrinv_v4 = 0.3333333333;
  //thr_v4 = 3; four_v4 = 4; five_v4 = 5; six_v4 = 6; seven_v4 = 7;
  half_v4   =simd_set_doublev4(0.5, 0.5, 0.5, 0.5);
  one_v4    =simd_set_doublev4(1.0, 1.0, 1.0, 1.0);
  two_v4    =simd_set_doublev4(2.0, 2.0, 2.0, 2.0);
  thrinv_v4 =simd_set_doublev4(0.333333, 0.333333, 0.333333, 0.333333);
  thr_v4    =simd_set_doublev4(3.0, 3.0, 3.0, 3.0);
  four_v4   =simd_set_doublev4(4.0, 4.0, 4.0, 4.0);
  five_v4   =simd_set_doublev4(5.0, 5.0, 5.0, 5.0);
  six_v4    =simd_set_doublev4(6.0, 6.0, 6.0, 6.0);
  seven_v4  =simd_set_doublev4(7.0, 7.0, 7.0, 7.0);

  doublev4 remd_v4[4];//remainder;
  remd_v4[0] = simd_set_doublev4(1.0, 1.0, 1.0, 1.0);
  remd_v4[1] = simd_set_doublev4(1.0, 0.0, 0.0, 0.0);
  remd_v4[2] = simd_set_doublev4(1.0, 1.0, 0.0, 0.0);
  remd_v4[3] = simd_set_doublev4(1.0, 1.0, 1.0, 0.0);
  
  //doublev4 test_remd_v4[4];//remainder;
  //test_remd_v4[0] = simd_set_doublev4(0.0, 0.0, 0.0, 0.0);//(1.0, 1.0, 1.0, 1.0);
  //test_remd_v4[1] = simd_set_doublev4(0.0, 1.0, 1.0, 1.0);//(1.0, 0.0, 0.0, 0.0);
  //test_remd_v4[2] = simd_set_doublev4(0.0, 0.0, 1.0, 1.0);//(1.0, 1.0, 0.0, 0.0);
  //test_remd_v4[3] = simd_set_doublev4(0.0, 0.0, 0.0, 1.0);//(1.0, 1.0, 1.0, 0.0);

  doublev4 Cele_v4, p_vdW1_v4, p_vdW1i_v4, powr_vdW1_v4;
  Cele_v4 = simd_set_doublev4(332.06371, 332.06371, 332.06371, 332.06371);//C_ele;  
  p_vdW1_v4 = l_pm.l28;  
  //p_vdW1i_v4 = one_v4 / p_vdW1_v4; //1.0 / l_pm.l28;
  p_vdW1i_v4 = simd_vdivd(one_v4, p_vdW1_v4); //1.0 / l_pm.l28;
  doublev4 power_onev4 = p_vdW1i_v4-one_v4;
  doublev4 power_twov4 = p_vdW1_v4 -two_v4;

  doublev4 Tap_v4, dTap_v4;

  doublev4 Tap0_v4, Tap1_v4, Tap2_v4, Tap3_v4, Tap4_v4, Tap5_v4, Tap6_v4, Tap7_v4; 
  //Tap0_v4 = l_pm.wksp_tap0;
  //Tap1_v4 = l_pm.wksp_tap1;
  //Tap2_v4 = l_pm.wksp_tap2;
  //Tap3_v4 = l_pm.wksp_tap3;
  //Tap4_v4 = l_pm.wksp_tap4;
  //Tap5_v4 = l_pm.wksp_tap5;
  //Tap6_v4 = l_pm.wksp_tap6;
  //Tap7_v4 = l_pm.wksp_tap7;
  
  Tap0_v4 = simd_set_doublev4(l_pm.wksp_tap0,l_pm.wksp_tap0,l_pm.wksp_tap0,l_pm.wksp_tap0);
  Tap1_v4 = simd_set_doublev4(l_pm.wksp_tap1,l_pm.wksp_tap1,l_pm.wksp_tap1,l_pm.wksp_tap1);
  Tap2_v4 = simd_set_doublev4(l_pm.wksp_tap2,l_pm.wksp_tap2,l_pm.wksp_tap2,l_pm.wksp_tap2);
  Tap3_v4 = simd_set_doublev4(l_pm.wksp_tap3,l_pm.wksp_tap3,l_pm.wksp_tap3,l_pm.wksp_tap3);
  Tap4_v4 = simd_set_doublev4(l_pm.wksp_tap4,l_pm.wksp_tap4,l_pm.wksp_tap4,l_pm.wksp_tap4);
  Tap5_v4 = simd_set_doublev4(l_pm.wksp_tap5,l_pm.wksp_tap5,l_pm.wksp_tap5,l_pm.wksp_tap5);
  Tap6_v4 = simd_set_doublev4(l_pm.wksp_tap6,l_pm.wksp_tap6,l_pm.wksp_tap6,l_pm.wksp_tap6);
  Tap7_v4 = simd_set_doublev4(l_pm.wksp_tap7,l_pm.wksp_tap7,l_pm.wksp_tap7,l_pm.wksp_tap7);

  doublev4 eng_vdwl_v4, eng_coul_v4, data_vdW_v4, data_ele_v4;
  eng_vdwl_v4 = eng_coul_v4 = 0;
  data_vdW_v4 = data_ele_v4 = 0;
  doublev4 virial0_v4, virial1_v4, virial2_v4;
  doublev4 virial3_v4, virial4_v4, virial5_v4;
  virial0_v4 = virial1_v4 = virial2_v4 = 0;
  virial3_v4 = virial4_v4 = virial5_v4 = 0;
  //lwpf_stop(BEFORE);

  //lwpf_start(CMP);
  for(ist = _MYID * ISTEP; ist < nlocal; ist += ISTEP * 64)
  {
    ied = ist + ISTEP;
    if(ied > nlocal)
      ied = nlocal;
    isz = ied - ist;

    pe_get(l_pm.index+ist, idx_st, sizeof(int)*isz);
    pe_get(l_pm.end_index+ist, idx_ed, sizeof(int)*isz);
    pe_get(l_pm.packed_atoms+ist, packed_atoms, sizeof(atom_pack_t)*isz);
    pe_get(l_pm.f[ist], fi[0], sizeof(double)*isz*4);
    dma_syn();

    for( ii = ist; ii < ied; ++ii) 
    {
      i = ii;
      int ioff = i - ist;
      if (packed_atoms[ioff].type < 0) continue;
      start_i = idx_st[ioff];
      end_i   = idx_ed[ioff];
      orig_i  = packed_atoms[ioff].orig_id;
      int type_i = packed_atoms[ioff].type;
      int jst, jed, jsz, jsz_to;
      
      doublev4 fi0_v4, fi1_v4, fi2_v4;
      fi0_v4 = fi1_v4 = fi2_v4 = 0;
      
      //qi_v4 = packed_atoms[ioff].q; 
      qi_v4 = simd_set_doublev4(packed_atoms[ioff].q, packed_atoms[ioff].q, 
                                packed_atoms[ioff].q, packed_atoms[ioff].q);

      for(jst = start_i; jst < end_i; jst+=SNSTEP)
      {
        jsz = SNSTEP;
        if(jst + SNSTEP > end_i)
          jsz = end_i - jst;
        
        get_reply = 0;
        dma_set_size(&get_desc, sizeof(far_neighbor_data_full)*jsz);
        dma(get_desc,l_pm.far_nbr_list_full+jst, pj_nbr);
        while(get_reply != 1);
        jsz_to = jsz & ~3;
        int pjst, tt;
            
        //lwpf_start(SEC1);
        for(pjst = 0; pjst < jsz_to; pjst+=4)
        {
          //transpose
          //lwpf_start(TRANS);
          doublev4 trans_nbr[3];//trans_nbr[8];
          p_nbr_v4 = &nbr_v4[pjst * 2];

          transpose8x4( p_nbr_v4[0],p_nbr_v4[2],p_nbr_v4[4],p_nbr_v4[6],
                        p_nbr_v4[1],p_nbr_v4[3],p_nbr_v4[5],p_nbr_v4[7], 
                        trans_nbr[0],trans_nbr[1],trans_nbr[2],r_ij_v4,
                        qj_v4, del0_v4, del1_v4, del2_v4);
          int *pt = &trans_nbr[0];
          int trc = vshuffd_rc(pt[7], pt[5], pt[3], pt[1]);//(typej3, typej2, typej1, typej0);

          //vshuffle
          alpha_v4 = simd_vshff(t_alpha_v4[type_i],t_alpha_v4[type_i],trc); 
          D_v4     = simd_vshff(t_D_v4[type_i],t_D_v4[type_i],trc); 
          gamma_v4 = simd_vshff(t_gamma_v4[type_i],t_gamma_v4[type_i],trc); 
          powgi_vdW1_v4     = simd_vshff(t_powgi_vdW1_v4[type_i],t_powgi_vdW1_v4[type_i],trc); 
          r_vdWinv_v4       = simd_vshff(t_r_vdWinv_v4[type_i],t_r_vdWinv_v4[type_i],trc); 
          alpha_div_rvdW_v4 = simd_vshff(t_alpha_div_rvdW_v4[type_i],t_alpha_div_rvdW_v4[type_i],trc); 
      
          //lwpf_stop(TRANS);
          //lwpf_start(VMAD1);
          //vectorize
          //rijinv_v4 = one_v4 / r_ij_v4;
          rijinv_v4 = simd_vdivd(one_v4, r_ij_v4);
          
          //Tap_v4 = Tap7_v4 * r_ij_v4 + Tap6_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap5_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap4_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap3_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap2_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap1_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap0_v4;//l_pm.wksp_tap0;//Tap0_v4;
          
          Tap_v4 = simd_vmad(Tap7_v4, r_ij_v4, Tap6_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap5_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap4_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap3_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap2_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap1_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap0_v4);//l_pm.wksp_tap0;//Tap0_v4;

          //dTap_v4 = seven_v4 * Tap7_v4 * r_ij_v4 + six_v4  * Tap6_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + five_v4 * Tap5_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + four_v4 * Tap4_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + thr_v4  * Tap3_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + two_v4  * Tap2_v4;
          //dTap_v4 += Tap1_v4 * rijinv_v4;
          dTap_v4 = simd_vmad(seven_v4, Tap7_v4 * r_ij_v4, six_v4  * Tap6_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, five_v4 * Tap5_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, four_v4 * Tap4_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, thr_v4  * Tap3_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, two_v4  * Tap2_v4);
          dTap_v4 = simd_vmad(Tap1_v4, rijinv_v4, dTap_v4);

          //lwpf_stop(VMAD1);

          //lwpf_start(MATH);
         
          //powr_vdW1_v4 = simd_vpowd(r_ij_v4, p_vdW1_v4);
          doublev4 ln_rij_v4 = simd_vlnd(r_ij_v4);
          powr_vdW1_v4 = simd_vexpd(p_vdW1_v4 * ln_rij_v4);

          doublev4 base_vdW = powr_vdW1_v4 + powgi_vdW1_v4;
          
          //fn13_v4 = simd_vpowd(powr_vdW1_v4 + powgi_vdW1_v4, p_vdW1i_v4);
          //fn13_v4 = simd_vpowd(base_vdW, p_vdW1i_v4);
          doublev4 ln_vdW_v4 = simd_vlnd(base_vdW);
          fn13_v4 =  simd_vexpd(p_vdW1i_v4 * ln_vdW_v4);

          exp2_v4 = simd_vexpd(half_v4 * alpha_v4 * (one_v4 - fn13_v4 * r_vdWinv_v4));
          exp1_v4 = exp2_v4 * exp2_v4;
          evdW_v4 = D_v4 * (exp1_v4 - two_v4 * exp2_v4);
          
          //data_vdW_v4 += Tap_v4 * evdW_v4;
          data_vdW_v4 = simd_vmad(Tap_v4, evdW_v4, data_vdW_v4);
          
          //lwpf_start(POW);
          //dfn13_v4 = simd_vpowd(powr_vdW1_v4 + powgi_vdW1_v4, p_vdW1i_v4-one_v4) 
          //            * simd_vpowd(r_ij_v4, p_vdW1_v4 - two_v4);
          
          //dfn13_v4 = simd_vpowd(base_vdW, power_onev4)  * ans_twov4;
          //dfn13_v4 = simd_vpowd(base_vdW, power_onev4)* simd_vpowd(r_ij_v4, power_twov4);
          dfn13_v4 = simd_vexpd(power_onev4 *ln_vdW_v4)* simd_vexpd(power_twov4 * ln_rij_v4);

          //lwpf_stop(POW);
          //lwpf_start(POW2) ;
          CEvd_v4 = dTap_v4 * evdW_v4 - Tap_v4 * D_v4 * alpha_div_rvdW_v4 
                    * (exp1_v4-exp2_v4) *dfn13_v4;
          
          //dr3gamij_1_v4 = (r_ij_v4 * r_ij_v4 * r_ij_v4 + gamma_v4);
          dr3gamij_1_v4 = simd_vmad(r_ij_v4, r_ij_v4 * r_ij_v4, gamma_v4);
          dr3gamij_3_v4 = simd_vpowd(dr3gamij_1_v4, thrinv_v4);
          
          //lwpf_stop(POW2) ;
          //dr3gamij_1inv_v4 = one_v4 / dr3gamij_1_v4;
          //dr3gamij_3inv_v4 = one_v4 / dr3gamij_3_v4;
          //lwpf_start(DIV);
          dr3gamij_1inv_v4 = simd_vdivd(one_v4, dr3gamij_1_v4);
          dr3gamij_3inv_v4 = simd_vdivd(one_v4, dr3gamij_3_v4);
          //lwpf_stop(DIV);

          //lwpf_stop(MATH);
          //
          //lwpf_start(VMAD);
          tmp_v4 = Tap_v4 * dr3gamij_3inv_v4;
          doublev4 Cqij_v4 = Cele_v4 * qi_v4 * qj_v4;
          data_ele_v4 += ecoul_v4 = Cqij_v4 * tmp_v4;
          
          CEclmb_v4 = Cqij_v4 * (dTap_v4 -  Tap_v4 * r_ij_v4 * dr3gamij_1inv_v4) 
                     * dr3gamij_3inv_v4;
            
          evdwl_v4 = Tap_v4 * evdW_v4;
          fpair_v4 = -(CEvd_v4 + CEclmb_v4);
          //tmpc_v4 = CEvd_v4 + CEclmb_v4;
          
          //eng_vdwl_v4 += (half_v4 *evdwl_v4);
          //eng_coul_v4 += (half_v4 *ecoul_v4);
          eng_vdwl_v4 = simd_vmad(half_v4, evdwl_v4, eng_vdwl_v4);
          eng_coul_v4 = simd_vmad(half_v4, ecoul_v4, eng_coul_v4);

          //virial0_v4 += (half_v4 * del0_v4 * del0_v4 * fpair_v4);
          //virial1_v4 += (half_v4 * del1_v4 * del1_v4 * fpair_v4);
          //virial2_v4 += (half_v4 * del2_v4 * del2_v4 * fpair_v4);
          //virial3_v4 += (half_v4 * del0_v4 * del1_v4 * fpair_v4);
          //virial4_v4 += (half_v4 * del0_v4 * del2_v4 * fpair_v4);
          //virial5_v4 += (half_v4 * del1_v4 * del2_v4 * fpair_v4);
          doublev4 half_fpair = half_v4 * fpair_v4;
          virial0_v4 = simd_vmad(half_fpair, del0_v4 * del0_v4, virial0_v4);
          virial1_v4 = simd_vmad(half_fpair, del1_v4 * del1_v4, virial1_v4);
          virial2_v4 = simd_vmad(half_fpair, del2_v4 * del2_v4, virial2_v4);
          virial3_v4 = simd_vmad(half_fpair, del0_v4 * del1_v4, virial3_v4);
          virial4_v4 = simd_vmad(half_fpair, del0_v4 * del2_v4, virial4_v4);
          virial5_v4 = simd_vmad(half_fpair, del1_v4 * del2_v4, virial5_v4);

          //fi0_v4 += (-tmpc_v4 * del0_v4);
          //fi1_v4 += (-tmpc_v4 * del1_v4);
          //fi2_v4 += (-tmpc_v4 * del2_v4);
          fi0_v4 = simd_vmad(fpair_v4, del0_v4, fi0_v4);
          fi1_v4 = simd_vmad(fpair_v4, del1_v4, fi1_v4);
          fi2_v4 = simd_vmad(fpair_v4, del2_v4, fi2_v4);
          //lwpf_stop(VMAD);
        }//for-pjst
        //lwpf_stop(SEC1);

        //lwpf_start(REMAINDER);
        if(jsz_to < jsz)
        {
          /**********remainder*********/
          //transpose
          doublev4 trans_nbr[3];//trans_nbr[8];
          p_nbr_v4 = &nbr_v4[jsz_to * 2];

          transpose8x4( p_nbr_v4[0],p_nbr_v4[2],p_nbr_v4[4],p_nbr_v4[6],
                        p_nbr_v4[1],p_nbr_v4[3],p_nbr_v4[5],p_nbr_v4[7], 
                        trans_nbr[0],trans_nbr[1],trans_nbr[2],r_ij_v4,
                        qj_v4, del0_v4, del1_v4, del2_v4);
          int *pt = &trans_nbr[0];
          int trc = vshuffd_rc(pt[7], pt[5], pt[3], pt[1]);//(typej3, typej2, typej1, typej0);

          //vshuffle
          alpha_v4 = simd_vshff(t_alpha_v4[type_i],t_alpha_v4[type_i],trc); 
          D_v4     = simd_vshff(t_D_v4[type_i],t_D_v4[type_i],trc); 
          gamma_v4 = simd_vshff(t_gamma_v4[type_i],t_gamma_v4[type_i],trc); 
          powgi_vdW1_v4     = simd_vshff(t_powgi_vdW1_v4[type_i],t_powgi_vdW1_v4[type_i],trc); 
          r_vdWinv_v4       = simd_vshff(t_r_vdWinv_v4[type_i],t_r_vdWinv_v4[type_i],trc); 
          alpha_div_rvdW_v4 = simd_vshff(t_alpha_div_rvdW_v4[type_i],t_alpha_div_rvdW_v4[type_i],trc); 
                
          //vectorize
          //rijinv_v4 = one_v4 / r_ij_v4;
          rijinv_v4 = simd_vdivd(one_v4, r_ij_v4);
          //rijinv_v4 = simd_vsellt(test_remd_v4[jsz &3], one_v4 / r_ij_v4, test_remd_v4[0]);//* remd_v4[jsz & 3];

          //Tap_v4 = Tap7_v4 * r_ij_v4 + Tap6_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap5_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap4_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap3_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap2_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap1_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap0_v4;//l_pm.wksp_tap0;//Tap0_v4;
          //dTap_v4 = seven_v4 * Tap7_v4 * r_ij_v4 + six_v4  * Tap6_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + five_v4 * Tap5_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + four_v4 * Tap4_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + thr_v4  * Tap3_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + two_v4  * Tap2_v4;
          //dTap_v4 += Tap1_v4 * rijinv_v4;
          Tap_v4 = simd_vmad(Tap7_v4, r_ij_v4, Tap6_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap5_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap4_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap3_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap2_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap1_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap0_v4);//l_pm.wksp_tap0;//Tap0_v4;
          dTap_v4 = simd_vmad(seven_v4, Tap7_v4 * r_ij_v4, six_v4  * Tap6_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, five_v4 * Tap5_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, four_v4 * Tap4_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, thr_v4  * Tap3_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, two_v4  * Tap2_v4);
          dTap_v4 = simd_vmad(Tap1_v4, rijinv_v4, dTap_v4);
          
          

          //powr_vdW1_v4 = simd_vpowd(r_ij_v4, p_vdW1_v4);
          doublev4 ln_rij_v4 = simd_vlnd(r_ij_v4);
          powr_vdW1_v4 = simd_vexpd(p_vdW1_v4 * ln_rij_v4);

          doublev4 base_vdW = powr_vdW1_v4 + powgi_vdW1_v4;

          //fn13_v4 = simd_vpowd(powr_vdW1_v4 + powgi_vdW1_v4, p_vdW1i_v4);
          doublev4 ln_vdW_v4 = simd_vlnd(base_vdW);
          fn13_v4 =  simd_vexpd(p_vdW1i_v4 * ln_vdW_v4);


          exp2_v4 = simd_vexpd(half_v4 * alpha_v4 * (one_v4 - fn13_v4 * r_vdWinv_v4));
          exp1_v4 = exp2_v4 * exp2_v4;
          evdW_v4 = D_v4 * (exp1_v4 - two_v4 * exp2_v4);
          //data_vdW_v4 += (Tap_v4 * evdW_v4) * remd_v4[jsz & 3];
          data_vdW_v4 = simd_vmad(Tap_v4, evdW_v4 * remd_v4[jsz & 3], data_vdW_v4);

          //dfn13_v4 = simd_vpowd(powr_vdW1_v4 + powgi_vdW1_v4, p_vdW1i_v4-one_v4) 
          //            * simd_vpowd(r_ij_v4, p_vdW1_v4 - two_v4);
          dfn13_v4 = simd_vexpd(power_onev4 *ln_vdW_v4)* simd_vexpd(power_twov4 * ln_rij_v4);


          CEvd_v4 = dTap_v4 * evdW_v4 - Tap_v4 * D_v4 * alpha_div_rvdW_v4 
                    * (exp1_v4-exp2_v4) *dfn13_v4;
          //dr3gamij_1_v4 = (r_ij_v4 * r_ij_v4 * r_ij_v4 + gamma_v4);
          dr3gamij_1_v4 = simd_vmad(r_ij_v4, r_ij_v4 * r_ij_v4, gamma_v4);
          dr3gamij_3_v4 = simd_vpowd(dr3gamij_1_v4, thrinv_v4);
          //dr3gamij_1inv_v4 = one_v4 / dr3gamij_1_v4;
          //dr3gamij_3inv_v4 = one_v4 / dr3gamij_3_v4;
          dr3gamij_1inv_v4 = simd_vdivd(one_v4, dr3gamij_1_v4);
          dr3gamij_3inv_v4 = simd_vdivd(one_v4, dr3gamij_3_v4);

          
          tmp_v4 = Tap_v4 * dr3gamij_3inv_v4;
          doublev4 Cqij_v4 = Cele_v4 * qi_v4 * qj_v4;
          data_ele_v4 += ecoul_v4 = Cqij_v4 * tmp_v4 * remd_v4[jsz & 3];
          
          CEclmb_v4 = Cqij_v4 * (dTap_v4 -  Tap_v4 * r_ij_v4 * dr3gamij_1inv_v4) 
                     * dr3gamij_3inv_v4;
            
          evdwl_v4 = Tap_v4 * evdW_v4;
          
          fpair_v4 = -(CEvd_v4 + CEclmb_v4)* remd_v4[jsz & 3];
          //tmpc_v4 = (CEvd_v4 + CEclmb_v4)*remd_v4[jsz & 3];
          
          //eng_vdwl_v4 += (half_v4 *evdwl_v4)*remd_v4[jsz & 3];
          //eng_coul_v4 += (half_v4 *ecoul_v4)*remd_v4[jsz & 3];
          eng_vdwl_v4 = simd_vmad(half_v4, evdwl_v4*remd_v4[jsz & 3], eng_vdwl_v4);
          eng_coul_v4 = simd_vmad(half_v4, ecoul_v4*remd_v4[jsz & 3], eng_coul_v4);


          //virial0_v4 += (half_v4 * del0_v4 * del0_v4 * fpair_v4);// * remd_v4[jsz & 3];
          //virial1_v4 += (half_v4 * del1_v4 * del1_v4 * fpair_v4);// * remd_v4[jsz & 3];
          //virial2_v4 += (half_v4 * del2_v4 * del2_v4 * fpair_v4);// * remd_v4[jsz & 3];
          //virial3_v4 += (half_v4 * del0_v4 * del1_v4 * fpair_v4);// * remd_v4[jsz & 3];
          //virial4_v4 += (half_v4 * del0_v4 * del2_v4 * fpair_v4);// * remd_v4[jsz & 3];
          //virial5_v4 += (half_v4 * del1_v4 * del2_v4 * fpair_v4);// * remd_v4[jsz & 3];
          doublev4 half_fpair = half_v4 * fpair_v4;
          virial0_v4 = simd_vmad(half_fpair, del0_v4 * del0_v4, virial0_v4);
          virial1_v4 = simd_vmad(half_fpair, del1_v4 * del1_v4, virial1_v4);
          virial2_v4 = simd_vmad(half_fpair, del2_v4 * del2_v4, virial2_v4);
          virial3_v4 = simd_vmad(half_fpair, del0_v4 * del1_v4, virial3_v4);
          virial4_v4 = simd_vmad(half_fpair, del0_v4 * del2_v4, virial4_v4);
          virial5_v4 = simd_vmad(half_fpair, del1_v4 * del2_v4, virial5_v4);


          //fi0_v4 += (-tmpc_v4 * del0_v4);
          //fi1_v4 += (-tmpc_v4 * del1_v4);
          //fi2_v4 += (-tmpc_v4 * del2_v4);
          fi0_v4 = simd_vmad(fpair_v4, del0_v4, fi0_v4);
          fi1_v4 = simd_vmad(fpair_v4, del1_v4, fi1_v4);
          fi2_v4 = simd_vmad(fpair_v4, del2_v4, fi2_v4);

        }
        //lwpf_stop(REMAINDER);



        //for( pj = jsz_to; pj < (jsz+4)&~3; pj++ ) 
        //for( pj = jsz_to; pj < jsz; pj++ ) 
        //{
        //  nbr_pj = &(pj_nbr[pj]);
        //  //j = nbr_pj->nbr;
        //  //if (nbr_pj->type < 0) continue;
        //  //orig_j  = nbr_pj->orig_id;
        //  //flag = 1;
        //  
        //  //lwpf_start(CEVD);
        //  r_ij = nbr_pj->d;
        //  int type_j = nbr_pj->type;
        //  double rijinv = 1 / r_ij;
        //  double r2ij = r_ij * r_ij;
        //  //twbp = &(l_tbpr[packed_atoms[ioff].type][nbr_pj->type]);
        //  //tbpt = tbpb[packed_atoms[ioff].type] + nbr_pj->type;
        //          
        //  Tap = l_pm.wksp_tap7 * r_ij + l_pm.wksp_tap6;
        //  Tap = Tap * r_ij + l_pm.wksp_tap5;
        //  Tap = Tap * r_ij + l_pm.wksp_tap4;
        //  Tap = Tap * r_ij + l_pm.wksp_tap3;
        //  Tap = Tap * r_ij + l_pm.wksp_tap2;
        //  Tap = Tap * r_ij + l_pm.wksp_tap1;
        //  Tap = Tap * r_ij + l_pm.wksp_tap0;

        //  dTap = 7 * l_pm.wksp_tap7 * r_ij + 6 * l_pm.wksp_tap6;
        //  dTap = dTap * r_ij + 5 * l_pm.wksp_tap5;
        //  dTap = dTap * r_ij + 4 * l_pm.wksp_tap4;
        //  dTap = dTap * r_ij + 3 * l_pm.wksp_tap3;
        //  dTap = dTap * r_ij + 2 * l_pm.wksp_tap2;
        //  dTap += l_pm.wksp_tap1 * rijinv;

        //  //lwpf_start(CEIFS);
        //  powr_vdW1 = p_powd(r_ij, p_vdW1);
        //  //powgi_vdW1 = tbpt->powgi_vdW1;
        //  powgi_vdW1 = p_powgi_vdW1[type_i * ntypes+type_j];
        //  fn13 = p_powd( powr_vdW1 + powgi_vdW1, p_vdW1i );
        //  //exp2 = p_expd( 0.5 * twbp->alpha * (1.0 - fn13 * tbpt->r_vdWinv) );
        //  exp2 = p_expd( 0.5 * p_alpha[type_i*ntypes+type_j] * (1.0 - fn13 * p_r_vdWinv[type_i*ntypes+type_j]) );
        //  exp1 = exp2 * exp2;
        //  //e_vdW = twbp->D * (exp1 - 2.0 * exp2);
        //  e_vdW = p_D[type_i*ntypes+type_j] * (exp1 - 2.0 * exp2);
        //  *data_vdW += Tap * e_vdW;

        //  dfn13 = p_powd( powr_vdW1 + powgi_vdW1, p_vdW1i-1.0) 
        //          * p_powd(r_ij, p_vdW1 - 2.0);
        //  //CEvd=dTap*e_vdW -Tap*twbp->D * (tbpt->alpha_div_rvdW)*(exp1-exp2)*dfn13;
        //  CEvd=dTap*e_vdW -Tap*p_D[type_i*ntypes+type_j] * (p_alpha_div_rvdW[type_i * ntypes+type_j])*(exp1-exp2)*dfn13;
        //  
        //  //lwpf_stop(CEIFS);
        //  //dr3gamij_1 = ( r_ij * r_ij * r_ij + twbp->gamma );
        //  dr3gamij_1 = ( r_ij * r_ij * r_ij + p_gamma[type_i*ntypes+type_j]);
        //  //lwpf_start(DIVF);
        //  dr3gamij_3 = p_powd( dr3gamij_1 , 0.33333333333333 );

        //  //lwpf_stop(DIVF);
        //  double dr3gamij_1inv = 1 / dr3gamij_1;
        //  double dr3gamij_3inv = 1 / dr3gamij_3;
        //  tmp = Tap * dr3gamij_3inv;

        //  *data_ele += e_ele = C_ele * packed_atoms[ioff].q * nbr_pj->q * tmp;
        //  CEclmb = C_ele * packed_atoms[ioff].q * nbr_pj->q *
        //           (dTap -  Tap * r_ij * dr3gamij_1inv ) * dr3gamij_3inv;
        //  
        //  //lwpf_stop(CEVD);

        //  //lwpf_start(EV);
        //  pe_vdw = Tap * (e_vdW + e_core + e_lg);
        //  f_tmp = -(CEvd + CEclmb);
        //  double v[6];
        //  double evdwl = pe_vdw, ecoul = e_ele, fpair = f_tmp;
        //  double *del = nbr_pj->dvec;
        //  
        //  *eng_vdwl += 0.5*evdwl;
        //  *eng_coul += 0.5*ecoul;

        //  v[0] = 0.5*del[0]*del[0]*fpair;
        //  v[1] = 0.5*del[1]*del[1]*fpair;
        //  v[2] = 0.5*del[2]*del[2]*fpair;
        //  v[3] = 0.5*del[0]*del[1]*fpair;
        //  v[4] = 0.5*del[0]*del[2]*fpair;
        //  v[5] = 0.5*del[1]*del[2]*fpair;
        //  virial[0] += v[0];
        //  virial[1] += v[1];
        //  virial[2] += v[2];
        //  virial[3] += v[3];
        //  virial[4] += v[4];
        //  virial[5] += v[5];

        //  double tmp_c = CEvd + CEclmb;
        //  fi[ioff][0] += -tmp_c * nbr_pj->dvec[0];
        //  fi[ioff][1] += -tmp_c * nbr_pj->dvec[1];
        //  fi[ioff][2] += -tmp_c * nbr_pj->dvec[2];
        //  //lwpf_stop(EV);
        //}//for-pj



        //lwpf_start(PJ);
        //for( pj = 0; pj < jsz; pj++ ) 
        //{
        //  nbr_pj = &(pj_nbr[pj]);
        //  j = nbr_pj->nbr;
        //  if (nbr_pj->type < 0) continue;
        //  orig_j  = nbr_pj->orig_id;
        //  flag = 1;
        //  
        //  //lwpf_start(CEVD);
        //  r_ij = nbr_pj->d;
        //  double rijinv = 1 / r_ij;
        //  double r2ij = r_ij * r_ij;
        //  twbp = &(l_tbpr[packed_atoms[ioff].type][nbr_pj->type]);
        //  tbpt = tbpb[packed_atoms[ioff].type] + nbr_pj->type;
        //          
        //  Tap = l_pm.wksp_tap7 * r_ij + l_pm.wksp_tap6;
        //  Tap = Tap * r_ij + l_pm.wksp_tap5;
        //  Tap = Tap * r_ij + l_pm.wksp_tap4;
        //  Tap = Tap * r_ij + l_pm.wksp_tap3;
        //  Tap = Tap * r_ij + l_pm.wksp_tap2;
        //  Tap = Tap * r_ij + l_pm.wksp_tap1;
        //  Tap = Tap * r_ij + l_pm.wksp_tap0;

        //  dTap = 7 * l_pm.wksp_tap7 * r_ij + 6 * l_pm.wksp_tap6;
        //  dTap = dTap * r_ij + 5 * l_pm.wksp_tap5;
        //  dTap = dTap * r_ij + 4 * l_pm.wksp_tap4;
        //  dTap = dTap * r_ij + 3 * l_pm.wksp_tap3;
        //  dTap = dTap * r_ij + 2 * l_pm.wksp_tap2;
        //  dTap += l_pm.wksp_tap1 * rijinv;


        //  //lwpf_start(CEIFS);
        //  powr_vdW1 = p_powd(r_ij, p_vdW1);
        //  powgi_vdW1 = tbpt->powgi_vdW1;
        //  fn13 = p_powd( powr_vdW1 + powgi_vdW1, p_vdW1i );
        //  exp2 = p_expd( 0.5 * twbp->alpha * (1.0 - fn13 * tbpt->r_vdWinv) );
        //  exp1 = exp2 * exp2;
        //  e_vdW = twbp->D * (exp1 - 2.0 * exp2);
        //  *data_vdW += Tap * e_vdW;

        //  dfn13 = p_powd( powr_vdW1 + powgi_vdW1, p_vdW1i-1.0) 
        //          * p_powd(r_ij, p_vdW1 - 2.0);
        //  CEvd=dTap*e_vdW -Tap*twbp->D * (tbpt->alpha_div_rvdW)*(exp1-exp2)*dfn13;
        //  
        //  //lwpf_stop(CEIFS);
        //  dr3gamij_1 = ( r_ij * r_ij * r_ij + twbp->gamma );
        //  //lwpf_start(DIVF);
        //  dr3gamij_3 = p_powd( dr3gamij_1 , 0.33333333333333 );

        //  //lwpf_stop(DIVF);
        //  double dr3gamij_1inv = 1 / dr3gamij_1;
        //  double dr3gamij_3inv = 1 / dr3gamij_3;
        //  tmp = Tap * dr3gamij_3inv;

        //  *data_ele += e_ele = C_ele * packed_atoms[ioff].q * nbr_pj->q * tmp;
        //  CEclmb = C_ele * packed_atoms[ioff].q * nbr_pj->q *
        //           (dTap -  Tap * r_ij * dr3gamij_1inv ) * dr3gamij_3inv;
        //  
        //  //lwpf_stop(CEVD);

        //  //lwpf_start(EV);
        //  pe_vdw = Tap * (e_vdW + e_core + e_lg);
        //  f_tmp = -(CEvd + CEclmb);
        //  double v[6];
        //  double evdwl = pe_vdw, ecoul = e_ele, fpair = f_tmp;
        //  double *del = nbr_pj->dvec;
        //  
        //  *eng_vdwl += 0.5*evdwl;
        //  *eng_coul += 0.5*ecoul;

        //  v[0] = 0.5*del[0]*del[0]*fpair;
        //  v[1] = 0.5*del[1]*del[1]*fpair;
        //  v[2] = 0.5*del[2]*del[2]*fpair;
        //  v[3] = 0.5*del[0]*del[1]*fpair;
        //  v[4] = 0.5*del[0]*del[2]*fpair;
        //  v[5] = 0.5*del[1]*del[2]*fpair;
        //  virial[0] += v[0];
        //  virial[1] += v[1];
        //  virial[2] += v[2];
        //  virial[3] += v[3];
        //  virial[4] += v[4];
        //  virial[5] += v[5];

        //  double tmp_c = CEvd + CEclmb;
        //  fi[ioff][0] += -tmp_c * nbr_pj->dvec[0];
        //  fi[ioff][1] += -tmp_c * nbr_pj->dvec[1];
        //  fi[ioff][2] += -tmp_c * nbr_pj->dvec[2];
        //  //lwpf_stop(EV);
        //}//for-pj




        //lwpf_stop(PJ);
      }//for-jst 
      
      simd_vsumd(fi0_v4);
      simd_vsumd(fi1_v4);
      simd_vsumd(fi2_v4);

      fi[ioff][0] += fi0_v4;
      fi[ioff][1] += fi1_v4;
      fi[ioff][2] += fi2_v4;

      virial[0] += packed_atoms[ioff].x[0] * fi[ioff][0];
      virial[1] += packed_atoms[ioff].x[1] * fi[ioff][1];
      virial[2] += packed_atoms[ioff].x[2] * fi[ioff][2];
      virial[3] += packed_atoms[ioff].x[0] * fi[ioff][1];
      virial[4] += packed_atoms[ioff].x[0] * fi[ioff][2];
      virial[5] += packed_atoms[ioff].x[1] * fi[ioff][2];

    }//for-ii
    pe_put(l_pm.f[ist], fi[0], sizeof(double)*isz*4);
    dma_syn();
  }//for-ist
  //lwpf_stop(CMP);

  //lwpf_start(SUM);
  //sum
  simd_vsumd(data_vdW_v4);
  simd_vsumd(data_ele_v4);
  simd_vsumd(eng_vdwl_v4);
  simd_vsumd(eng_coul_v4);
  simd_vsumd(virial0_v4);
  simd_vsumd(virial1_v4);
  simd_vsumd(virial2_v4);
  simd_vsumd(virial3_v4);
  simd_vsumd(virial4_v4);
  simd_vsumd(virial5_v4);
  //reduce
  *data_vdW += data_vdW_v4;
  *data_ele += data_ele_v4;
  *eng_vdwl += eng_vdwl_v4;
  *eng_coul += eng_coul_v4;
  virial[0] += virial0_v4;
  virial[1] += virial1_v4;
  virial[2] += virial2_v4;
  virial[3] += virial3_v4;
  virial[4] += virial4_v4;
  virial[5] += virial5_v4;

  //lwpf_stop(SUM);
  //lwpf_start(REDUCE);
  reg_reduce_inplace_doublev4(eng_virial, 4);

  if(_MYID == 0)
  {
    pe_put(l_pm.packed_eng, eng_virial, sizeof(double)*16);
    dma_syn();
  }
  //lwpf_stop(REDUCE);
  //lwpf_stop(ALL);
  //lwpf_exit(NONBOND);
}

#define EVFLAG 1
#include"reaxc_nonbonded_cpe.h"
#undef EVFLAG

#define EVFLAG 0
#include"reaxc_nonbonded_cpe.h"
#undef EVFLAG

#endif
