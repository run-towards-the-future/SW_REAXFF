#include <math.h>
#include "reaxc_defs_sunway.h"
#include "reaxc_ctypes_sunway.h"
#include "reaxc_inlines_sw64.h"
#include "reaxc_multi_body_sw64.h"
#include "stdio.h"
#include "stdlib.h"
#include "gptl.h"
#include "sunway.h"
#include "simd.h"
#define ISTEP 32
#define MAX_LEN 20

#ifdef MPE
#include <mpi.h>
#include <athread.h>
//#define LWPF_UNITS U(BOFUNC)
//#include "lwpf2.h"


void Merge_Bonds_Atom_Energy_C_New(merge_bonds_eng_t *param) 
{
  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;

  int i, j, pj, type_i, type_j;
  double Delta_lpcorr, dfvl;
  double e_lp, expvd2, inv_expvd2, dElp, CElp, DlpVi;
  double e_lph, Di, vov3, deahu2dbo, deahu2dsbo;
  double e_ov, CEover1, CEover2, CEover3, CEover4;
  double exp_ovun1, exp_ovun2, sum_ovun1, sum_ovun2;
  double exp_ovun2n, exp_ovun6, exp_ovun8;
  double inv_exp_ovun1, inv_exp_ovun2, inv_exp_ovun2n, inv_exp_ovun8;
  double e_un, CEunder1, CEunder2, CEunder3, CEunder4;
  double p_lp2, p_lp3;
  double p_ovun2, p_ovun3, p_ovun4, p_ovun5, p_ovun6, p_ovun7, p_ovun8;
  double eng_tmp;
  int numbonds;

  single_body_parameters *sbp_i, *sbp_j;
  //single_body_parameters *sbp_i;
  two_body_parameters *twbp;
  bond_data *pbond;
  bond_order_data *bo_ij;
  reax_list *bonds = (*lists) + BONDS;
  
  double ebond, pow_BOs_be2, exp_be12, CEbo;
  double gp3, gp4, gp7, gp10, gp37;
  double exphu, exphua1, exphub1, exphuov, hulpov, estriph;
  double decobdbo, decobdboua, decobdboub;

  double *Cdbo_list     = bonds->Cdbo_list;
  double *Cdbopi_list   = bonds->Cdbopi_list;
  double *Cdbopi2_list  = bonds->Cdbopi2_list;
  double *BO_list       = bonds->BO_list;
  rvec2 *BOpi_list     = bonds->BOpi_list;

  /* Initialize parameters */
  p_lp3 = system->reax_param.gp.l[5];
  p_ovun3 = system->reax_param.gp.l[32];
  p_ovun4 = system->reax_param.gp.l[31];
  p_ovun6 = system->reax_param.gp.l[6];
  p_ovun7 = system->reax_param.gp.l[8];
  p_ovun8 = system->reax_param.gp.l[9];
  
  gp3 = system->reax_param.gp.l[3];
  gp4 = system->reax_param.gp.l[4];
  gp7 = system->reax_param.gp.l[7];
  gp10 = system->reax_param.gp.l[10];
  gp37 = (int) system->reax_param.gp.l[37];


  double ftmp, ftmp2;
  rvec4 fc_tmp;
  fc_tmp[0] = fc_tmp[1] = fc_tmp[2] = fc_tmp[3] = 0.0;
  double eng_vdwl = 0, eng_eov = 0, eng_eun = 0;
  double eng_ebond = 0, eng_elp = 0;
 
  e_lph   = 0.0;
  e_lp    = 0.0;
  e_un    = 0.0;
  e_ov    = 0.0;
  ebond   = 0.0;
  estriph = 0.0;

  double q, en_tmp = 0;
  data->my_en.e_pol = 0.0;

  for(i = 0; i < system->n; ++i) 
  {
    /* set the parameter pointer */
    type_i = system->packed_atoms[i].type;
    if (type_i < 0) continue;
    sbp_i = &(system->reax_param.sbp[ type_i ]);

    /*Compute_Polarization_Energy*/
    q = system->packed_atoms[i].q;
    en_tmp += KCALpMOL_to_EV * (sbp_i->chi * q + (sbp_i->eta / 2.) * SQR(q));
    data->my_en.e_pol += en_tmp;

    /* lone-pair Energy */
    p_lp2 = sbp_i->p_lp2;
    expvd2 = exp( -75 * workspace->Delta_lp[i] );
    inv_expvd2 = 1. / (1. + expvd2 );

    numbonds = 0;
    p_ovun2 = sbp_i->p_ovun2;
    sum_ovun1 = sum_ovun2 = 0;
    if( sbp_i->mass > 21.0 )
      dfvl = 0.0;
    else dfvl = 1.0; // only for 1st-row elements

    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
    {
      numbonds ++;
      j = bonds->select.bond_list[pj].nbr;
      type_j = system->packed_atoms[j].type;
	    if (type_j < 0) continue;
      bo_ij = &(bonds->bo_data_list[pj]);
      twbp = &(system->reax_param.tbp[ type_i ][ type_j ]);

      sum_ovun1 += twbp->p_ovun1 * twbp->De_s * BO_list[pj];
      sum_ovun2 += (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j])*
        (BOpi_list[pj][0] + BOpi_list[pj][1]);
    }
    
    
    dElp = p_lp2 * inv_expvd2 +
      75 * p_lp2 * workspace->Delta_lp[i] * expvd2 * SQR(inv_expvd2);
    CElp = dElp * workspace->dDelta_lp[i];

        
    if( p_lp3 > 0.001 && !strcmp(system->reax_param.sbp[type_i].name, "C") )
    {
      for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) 
      {
        j = bonds->select.bond_list[pj].nbr;
        type_j = system->packed_atoms[j].type;
	      if (type_j < 0) continue;

        if( !strcmp( system->reax_param.sbp[type_j].name, "C" ) ) 
        {
          twbp = &( system->reax_param.tbp[type_i][type_j]);
          bo_ij = &(bonds->bo_data_list[pj]);
          Di = workspace->Delta[i];
          vov3 = BO_list[pj] - Di - 0.040*pow(Di, 4.);

          if( vov3 > 3. ) 
          {
            e_lph += p_lp3 * SQR(vov3-3.0);

            deahu2dbo = 2.*p_lp3*(vov3 - 3.);
            deahu2dsbo = 2.*p_lp3*(vov3 - 3.)*(-1. - 0.16*pow(Di, 3.));

            Cdbo_list[pj] += deahu2dbo;
            workspace->fCdDelta[i][3] += deahu2dsbo;
          }
        }
      }//for
    }//if

    exp_ovun1 = p_ovun3 * exp( p_ovun4 * sum_ovun2 );
    inv_exp_ovun1 = 1.0 / (1 + exp_ovun1);
    Delta_lpcorr  = workspace->Delta[i] -
      (dfvl * workspace->Delta_lp_temp[i]) * inv_exp_ovun1;

    exp_ovun2 = exp( p_ovun2 * Delta_lpcorr );
    inv_exp_ovun2 = 1.0 / (1.0 + exp_ovun2);

    DlpVi = 1.0 / (Delta_lpcorr + sbp_i->valency + 1e-8);
    CEover1 = Delta_lpcorr * DlpVi * inv_exp_ovun2;

    e_ov += sum_ovun1 * CEover1;

    CEover2 = sum_ovun1 * DlpVi * inv_exp_ovun2 *
      (1.0 - Delta_lpcorr * ( DlpVi + p_ovun2 * exp_ovun2 * inv_exp_ovun2 ));

    CEover3 = CEover2 * (1.0 - dfvl * workspace->dDelta_lp[i] * inv_exp_ovun1 );

    CEover4 = CEover2 * (dfvl * workspace->Delta_lp_temp[i]) *
      p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1);


    /* under-coordination potential */
    p_ovun2 = sbp_i->p_ovun2;
    p_ovun5 = sbp_i->p_ovun5;

    exp_ovun2n = 1.0 / exp_ovun2;
    exp_ovun6 = exp( p_ovun6 * Delta_lpcorr );
    exp_ovun8 = p_ovun7 * exp(p_ovun8 * sum_ovun2);
    inv_exp_ovun2n = 1.0 / (1.0 + exp_ovun2n);
    inv_exp_ovun8 = 1.0 / (1.0 + exp_ovun8);

    e_un = 0.0;
    
    workspace->fCdDelta[i][3] += CEover3;   // OvCoor - 2nd term

    if (numbonds > 0 || control->enobondsflag)
    {
      data->my_en.e_un += e_un =
        -p_ovun5 * (1.0 - exp_ovun6) * inv_exp_ovun2n * inv_exp_ovun8;

      CEunder1 = inv_exp_ovun2n *
        ( p_ovun5 * p_ovun6 * exp_ovun6 * inv_exp_ovun8 +
          p_ovun2 * e_un * exp_ovun2n );
      CEunder2 = -e_un * p_ovun8 * exp_ovun8 * inv_exp_ovun8;
      CEunder3 = CEunder1 * (1.0 - dfvl*workspace->dDelta_lp[i]*inv_exp_ovun1);
      CEunder4 = CEunder1 * (dfvl*workspace->Delta_lp_temp[i]) *
        p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1) + CEunder2;
      
      workspace->fCdDelta[i][3] += CEunder3;  // UnCoor - 1st term
      workspace->fCdDelta[i][3] += CElp;  // lp - 1st term
      e_lp += p_lp2 * workspace->Delta_lp[i] * inv_expvd2;

      if( system->evflag)
      {
        system->eng_vdwl += e_un;
      }
    }
    else
    {
      CEunder1 = inv_exp_ovun2n *
        ( p_ovun5 * p_ovun6 * exp_ovun6 * inv_exp_ovun8 +
          p_ovun2 * e_un * exp_ovun2n );
      CEunder2 = -e_un * p_ovun8 * exp_ovun8 * inv_exp_ovun8;
      CEunder3 = CEunder1 * (1.0 - dfvl*workspace->dDelta_lp[i]*inv_exp_ovun1);
      CEunder4 = CEunder1 * (dfvl*workspace->Delta_lp_temp[i]) *
        p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1) + CEunder2;

    }
    
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) 
    {
      pbond = &(bonds->select.bond_list[pj]);
      j = pbond->nbr;
      type_j = system->packed_atoms[j].type;
      twbp = &( system->reax_param.tbp[type_i][type_j] );
      bo_ij= &(bonds->bo_data_list[pj]);

      Cdbo_list[pj] += CEover1 * twbp->p_ovun1 * twbp->De_s;// OvCoor-1st
      ftmp2 = (1.0 - dfvl*workspace->dDelta_lp[j]) *(BOpi_list[pj][0] + BOpi_list[pj][1]);
      workspace->fCdDelta[j][3] += CEover4 * ftmp2; 
      ftmp = CEover4 *(workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]); // OvCoor-3b
      Cdbopi_list[pj] += ftmp;
      Cdbopi2_list[pj] += ftmp;
      workspace->fCdDelta[j][3] += CEunder4 * ftmp2; 
      fc_tmp[3] = CEunder4 *(workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]);  
      Cdbopi_list[pj] += fc_tmp[3];
      Cdbopi2_list[pj] += fc_tmp[3];

      /////////////for Bonds/////////////
      int flag = 1;
      if( system->packed_atoms[i].orig_id > system->packed_atoms[j].orig_id )
	      flag = 0;
      if( system->packed_atoms[i].orig_id == system->packed_atoms[j].orig_id ) 
      {
        if (system->packed_atoms[j].x[2] <  system->packed_atoms[i].x[2]) flag = 0;;
      	if (system->packed_atoms[j].x[2] == system->packed_atoms[i].x[2] &&
      	    system->packed_atoms[j].x[1] <  system->packed_atoms[i].x[1]) flag = 0;;
        if (system->packed_atoms[j].x[2] == system->packed_atoms[i].x[2] &&
      	    system->packed_atoms[j].x[1] == system->packed_atoms[i].x[1] &&
      	    system->packed_atoms[j].x[0] <  system->packed_atoms[i].x[0]) flag = 0;
      }
      if(flag)
      {
        /* set the pointers */
        sbp_i = &( system->reax_param.sbp[type_i] );
        sbp_j = &( system->reax_param.sbp[type_j] );

        /* calculate the constants */
        if (bo_ij->BO_s == 0.0) pow_BOs_be2 = 0.0;
        else pow_BOs_be2 = pow( bo_ij->BO_s, twbp->p_be2 );
        exp_be12 = exp( twbp->p_be1 * ( 1.0 - pow_BOs_be2 ) );
        CEbo = -twbp->De_s * exp_be12 *
	          ( 1.0 - twbp->p_be1 * twbp->p_be2 * pow_BOs_be2 );

        /* calculate the Bond Energy */
        ebond +=
	        -twbp->De_s * bo_ij->BO_s * exp_be12
	        -twbp->De_p * BOpi_list[pj][0]
	        -twbp->De_pp * BOpi_list[pj][1];

        /* calculate derivatives of Bond Orders */
        Cdbo_list[pj]     += CEbo;
        Cdbopi_list[pj]   -= (CEbo + twbp->De_p);
        Cdbopi2_list[pj]  -= (CEbo + twbp->De_pp);

        /* Stabilisation terminal triple bond */
        if(BO_list[pj] >= 1.00 ) 
        {
	        if( gp37 == 2 ||
	        (sbp_i->mass == 12.0000 && sbp_j->mass == 15.9990) ||
	        (sbp_j->mass == 12.0000 && sbp_i->mass == 15.9990) ) 
          {
	          exphu = exp( -gp7 * SQR(BO_list[pj] - 2.50) );
	          exphua1 = exp(-gp3 * (workspace->bo_dboc[i][0]-BO_list[pj]));
	          exphub1 = exp(-gp3 * (workspace->bo_dboc[j][0]-BO_list[pj]));
	          exphuov = exp(gp4 * (workspace->Delta[i] + workspace->Delta[j]));
	          hulpov = 1.0 / (1.0 + 25.0 * exphuov);

	          estriph += gp10 * exphu * hulpov * (exphua1 + exphub1);

	          decobdbo = gp10 * exphu * hulpov * (exphua1 + exphub1) *
	            ( gp3 - 2.0 * gp7 * (BO_list[pj]-2.50) );

	          decobdboua = -gp10 * exphu * hulpov *
	            (gp3*exphua1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));
	          decobdboub = -gp10 * exphu * hulpov *
	            (gp3*exphub1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));

	          Cdbo_list[pj] += decobdbo;
	          workspace->fCdDelta[i][3] += decobdboua;
	          workspace->fCdDelta[j][3] += decobdboub;
	        }//if
        }//if
      }//if-flag
    }//for
  }//for-i

  if(system->evflag)
  {

    system->eng_vdwl += (e_lph+e_lp+e_ov+ebond+estriph);
    system->eng_coul += en_tmp;
  }
  data->my_en.e_lp += (e_lph+e_lp);
  data->my_en.e_ov += e_ov ;
  data->my_en.e_bond += (ebond+estriph);
}


void Atom_Energy_C(merge_bonds_eng_t *param) 
{
  Merge_Bonds_Atom_Energy_C_New(param);
  return;

  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;

  int i, j, pj, type_i, type_j;
  double Delta_lpcorr, dfvl;
  double e_lp, expvd2, inv_expvd2, dElp, CElp, DlpVi;
  double e_lph, Di, vov3, deahu2dbo, deahu2dsbo;
  double e_ov, CEover1, CEover2, CEover3, CEover4;
  double exp_ovun1, exp_ovun2, sum_ovun1, sum_ovun2;
  double exp_ovun2n, exp_ovun6, exp_ovun8;
  double inv_exp_ovun1, inv_exp_ovun2, inv_exp_ovun2n, inv_exp_ovun8;
  double e_un, CEunder1, CEunder2, CEunder3, CEunder4;
  double p_lp2, p_lp3;
  double p_ovun2, p_ovun3, p_ovun4, p_ovun5, p_ovun6, p_ovun7, p_ovun8;
  double eng_tmp;
  int numbonds;

  single_body_parameters *sbp_i;
  two_body_parameters *twbp;
  bond_data *pbond;
  bond_order_data *bo_ij;
  reax_list *bonds = (*lists) + BONDS;
  
  double *Cdbo_list     = bonds->Cdbo_list;
  double *Cdbopi_list   = bonds->Cdbopi_list;
  double *Cdbopi2_list  = bonds->Cdbopi2_list;
  double *BO_list       = bonds->BO_list;
  rvec2 *BOpi_list     = bonds->BOpi_list;

  /* Initialize parameters */
  p_lp3 = system->reax_param.gp.l[5];
  p_ovun3 = system->reax_param.gp.l[32];
  p_ovun4 = system->reax_param.gp.l[31];
  p_ovun6 = system->reax_param.gp.l[6];
  p_ovun7 = system->reax_param.gp.l[8];
  p_ovun8 = system->reax_param.gp.l[9];


  for( i = 0; i < system->n; ++i ) 
  {
    /* set the parameter pointer */
    type_i = system->packed_atoms[i].type;
    if (type_i < 0) continue;
    sbp_i = &(system->reax_param.sbp[ type_i ]);

    /* lone-pair Energy */
    p_lp2 = sbp_i->p_lp2;
    expvd2 = exp( -75 * workspace->Delta_lp[i] );
    inv_expvd2 = 1. / (1. + expvd2 );

    numbonds = 0;
    e_lp = 0.0;
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
      numbonds ++;

    /* calculate the energy */
    if (numbonds > 0 || control->enobondsflag)
      data->my_en.e_lp += e_lp =
        p_lp2 * workspace->Delta_lp[i] * inv_expvd2;

    dElp = p_lp2 * inv_expvd2 +
      75 * p_lp2 * workspace->Delta_lp[i] * expvd2 * SQR(inv_expvd2);
    CElp = dElp * workspace->dDelta_lp[i];

    if (numbonds > 0 || control->enobondsflag)
      workspace->fCdDelta[i][3] += CElp;  // lp - 1st term

    /* tally into per-atom energy */
    if( system->evflag)
    {
      //system->pair_ptr->ev_tally(i,i,system->n,1,e_lp,0.0,0.0,0.0,0.0,0.0);
      system->eng_vdwl += e_lp;
    }

    /* correction for C2 */
    if( p_lp3 > 0.001 && !strcmp(system->reax_param.sbp[type_i].name, "C") )
    {
      for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) 
      {
        j = bonds->select.bond_list[pj].nbr;
        type_j = system->packed_atoms[j].type;
	      if (type_j < 0) continue;

        if( !strcmp( system->reax_param.sbp[type_j].name, "C" ) ) 
        {
          twbp = &( system->reax_param.tbp[type_i][type_j]);
          //bo_ij = &( bonds->select.bond_list[pj].bo_data );
          bo_ij = &(bonds->bo_data_list[pj]);
          Di = workspace->Delta[i];
          //vov3 = bo_ij->BO - Di - 0.040*pow(Di, 4.);
          vov3 = BO_list[pj] - Di - 0.040*pow(Di, 4.);

          if( vov3 > 3. ) 
          {
            data->my_en.e_lp += e_lph = p_lp3 * SQR(vov3-3.0);

            deahu2dbo = 2.*p_lp3*(vov3 - 3.);
            deahu2dsbo = 2.*p_lp3*(vov3 - 3.)*(-1. - 0.16*pow(Di, 3.));

            //bo_ij->Cdbo += deahu2dbo;
            Cdbo_list[pj] += deahu2dbo;
            workspace->fCdDelta[i][3] += deahu2dsbo;

            /* tally into per-atom energy */
            if( system->evflag)
            {
              //system->pair_ptr->ev_tally(i,j,system->n,1,e_lph,0.0,0.0,0.0,0.0,0.0);
              system->eng_vdwl += e_lph;
            }

          }
        }
      }
    }//if
  }

    for( i = 0; i < system->n; ++i ) 
    {
      type_i = system->packed_atoms[i].type;
      if (type_i < 0) continue;
      sbp_i = &(system->reax_param.sbp[ type_i ]);

      /* over-coordination energy */
      if( sbp_i->mass > 21.0 )
        dfvl = 0.0;
      else dfvl = 1.0; // only for 1st-row elements

      p_ovun2 = sbp_i->p_ovun2;
      sum_ovun1 = sum_ovun2 = 0;
      for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) 
      {
        j = bonds->select.bond_list[pj].nbr;
        type_j = system->packed_atoms[j].type;
	      if (type_j < 0) continue;
        //bo_ij = &(bonds->select.bond_list[pj].bo_data);
        bo_ij = &(bonds->bo_data_list[pj]);
        twbp = &(system->reax_param.tbp[ type_i ][ type_j ]);

        //sum_ovun1 += twbp->p_ovun1 * twbp->De_s * bo_ij->BO;
        sum_ovun1 += twbp->p_ovun1 * twbp->De_s * BO_list[pj];
        sum_ovun2 += (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j])*
          ( BOpi_list[pj][0] + BOpi_list[pj][1]);
          //( bo_ij->BO_pi + bo_ij->BO_pi2 );

      }

    exp_ovun1 = p_ovun3 * exp( p_ovun4 * sum_ovun2 );
    inv_exp_ovun1 = 1.0 / (1 + exp_ovun1);
    Delta_lpcorr  = workspace->Delta[i] -
      (dfvl * workspace->Delta_lp_temp[i]) * inv_exp_ovun1;

    exp_ovun2 = exp( p_ovun2 * Delta_lpcorr );
    inv_exp_ovun2 = 1.0 / (1.0 + exp_ovun2);

    DlpVi = 1.0 / (Delta_lpcorr + sbp_i->valency + 1e-8);
    CEover1 = Delta_lpcorr * DlpVi * inv_exp_ovun2;

    data->my_en.e_ov += e_ov = sum_ovun1 * CEover1;

    CEover2 = sum_ovun1 * DlpVi * inv_exp_ovun2 *
      (1.0 - Delta_lpcorr * ( DlpVi + p_ovun2 * exp_ovun2 * inv_exp_ovun2 ));

    CEover3 = CEover2 * (1.0 - dfvl * workspace->dDelta_lp[i] * inv_exp_ovun1 );

    CEover4 = CEover2 * (dfvl * workspace->Delta_lp_temp[i]) *
      p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1);


    /* under-coordination potential */
    p_ovun2 = sbp_i->p_ovun2;
    p_ovun5 = sbp_i->p_ovun5;

    exp_ovun2n = 1.0 / exp_ovun2;
    exp_ovun6 = exp( p_ovun6 * Delta_lpcorr );
    exp_ovun8 = p_ovun7 * exp(p_ovun8 * sum_ovun2);
    inv_exp_ovun2n = 1.0 / (1.0 + exp_ovun2n);
    inv_exp_ovun8 = 1.0 / (1.0 + exp_ovun8);

    numbonds = 0;
    e_un = 0.0;
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
      numbonds ++;

    if (numbonds > 0 || control->enobondsflag)
      data->my_en.e_un += e_un =
        -p_ovun5 * (1.0 - exp_ovun6) * inv_exp_ovun2n * inv_exp_ovun8;

    CEunder1 = inv_exp_ovun2n *
      ( p_ovun5 * p_ovun6 * exp_ovun6 * inv_exp_ovun8 +
        p_ovun2 * e_un * exp_ovun2n );
    CEunder2 = -e_un * p_ovun8 * exp_ovun8 * inv_exp_ovun8;
    CEunder3 = CEunder1 * (1.0 - dfvl*workspace->dDelta_lp[i]*inv_exp_ovun1);
    CEunder4 = CEunder1 * (dfvl*workspace->Delta_lp_temp[i]) *
      p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1) + CEunder2;

    /* tally into per-atom energy */
    if( system->evflag) 
    {
      eng_tmp = e_ov;
      if (numbonds > 0 || control->enobondsflag)
        eng_tmp += e_un;
      //system->pair_ptr->ev_tally(i,i,system->n,1,eng_tmp,0.0,0.0,0.0,0.0,0.0);
      system->eng_vdwl += eng_tmp;
    }

    /* forces */
    workspace->fCdDelta[i][3] += CEover3;   // OvCoor - 2nd term
    if (numbonds > 0 || control->enobondsflag)
      workspace->fCdDelta[i][3] += CEunder3;  // UnCoor - 1st term

    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) 
    {
      pbond = &(bonds->select.bond_list[pj]);
      j = pbond->nbr;
      //bo_ij = &(pbond->bo_data);
      bo_ij= &(bonds->bo_data_list[pj]);
      twbp  = &(system->reax_param.tbp[ system->packed_atoms[i].type ]
                [system->packed_atoms[pbond->nbr].type]);


      //bo_ij->Cdbo += CEover1 * twbp->p_ovun1 * twbp->De_s;// OvCoor-1st
      Cdbo_list[pj] += CEover1 * twbp->p_ovun1 * twbp->De_s;// OvCoor-1st

      workspace->fCdDelta[j][3] += CEover4 * (1.0 - dfvl*workspace->dDelta_lp[j]) *
        (BOpi_list[pj][0] + BOpi_list[pj][1]); // OvCoor-3a
        //(bo_ij->BO_pi + bo_ij->BO_pi2); // OvCoor-3a

      //bo_ij->Cdbopi += CEover4 *
      Cdbopi_list[pj] += CEover4 *
        (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]); // OvCoor-3b
      //bo_ij->Cdbopi2 += CEover4 *
      Cdbopi2_list[pj] += CEover4 *
        (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]);  // OvCoor-3b


      workspace->fCdDelta[j][3] += CEunder4 * (1.0 - dfvl*workspace->dDelta_lp[j]) *
        (BOpi_list[pj][0]+BOpi_list[pj][1]);   // UnCoor - 2a
        //(bo_ij->BO_pi + bo_ij->BO_pi2);   // UnCoor - 2a

      //bo_ij->Cdbopi += CEunder4 *
      Cdbopi_list[pj] += CEunder4 *
        (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]);  // UnCoor-2b
      //bo_ij->Cdbopi2 += CEunder4 *
      Cdbopi2_list[pj] += CEunder4 *
        (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]);  // UnCoor-2b

    }

  }
}
#endif

#ifdef CPE
#include "STUBS/mpi.h"
#include <dma.h>
#include "slave.h"
#define DMA_FAST
#include "poly_math.h"
#include "dma_macros.h"
//#define LWPF_UNIT U(BOFUNC)
//#define LWPF_KERNELS K(ALL) K(BEFORE) K(CMP) K(LOG)   
//#include "lwpf2.h"
#endif
