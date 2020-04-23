#include <math.h>
#include "reaxc_defs_sunway.h"
#include "reaxc_ctypes_sunway.h"
#include "reaxc_inlines_sw64.h"
#include "reaxc_torsion_angles_sw64.h"
#include "stdio.h"
#include "stdlib.h"
#include "gptl.h"
#include "sunway.h"
#include "simd.h"
#include "sleef_math.h"
#define MIN_SINE 1e-10
#define JSTEP 8
#define BLEN  20

#ifdef MPE
#include <mpi.h>
#include <athread.h>
//#define LWPF_UNITS U(TORSION)
//#include "lwpf2.h"
double Calculate_Omega( rvec dvec_ij, double r_ij,
                      rvec dvec_jk, double r_jk,
                      rvec dvec_kl, double r_kl,
                      rvec dvec_li, double r_li,
                      three_body_interaction_data *p_ijk,
                      three_body_interaction_data *p_jkl,
                      rvec dcos_omega_di, rvec dcos_omega_dj,
                      rvec dcos_omega_dk, rvec dcos_omega_dl,
                      output_controls *out_control )
{
  double unnorm_cos_omega, unnorm_sin_omega, omega;
  double sin_ijk, cos_ijk, sin_jkl, cos_jkl;
  double htra, htrb, htrc, hthd, hthe, hnra, hnrc, hnhd, hnhe;
  double arg, poem, tel;
  rvec cross_jk_kl;

  sin_ijk = sin(p_ijk->theta);
  cos_ijk = cos(p_ijk->theta);
  sin_jkl = sin(p_jkl->theta);
  cos_jkl = cos(p_jkl->theta);

  /* omega */
  unnorm_cos_omega = -rvec_Dot(dvec_ij, dvec_jk) * rvec_Dot(dvec_jk, dvec_kl) +
                    SQR(r_jk) *  rvec_Dot(dvec_ij, dvec_kl);

  rvec_Cross(cross_jk_kl, dvec_jk, dvec_kl);
  unnorm_sin_omega = -r_jk * rvec_Dot( dvec_ij, cross_jk_kl );

  omega = atan2( unnorm_sin_omega, unnorm_cos_omega );

  htra = r_ij + cos_ijk * (r_kl * cos_jkl - r_jk);
  htrb = r_jk - r_ij * cos_ijk - r_kl * cos_jkl;
  htrc = r_kl + cos_jkl * (r_ij * cos_ijk - r_jk);
  hthd = r_ij * sin_ijk * (r_jk - r_kl * cos_jkl);
  hthe = r_kl * sin_jkl * (r_jk - r_ij * cos_ijk);
  hnra = r_kl * sin_ijk * sin_jkl;
  hnrc = r_ij * sin_ijk * sin_jkl;
  hnhd = r_ij * r_kl * cos_ijk * sin_jkl;
  hnhe = r_ij * r_kl * sin_ijk * cos_jkl;

  poem = 2.0 * r_ij * r_kl * sin_ijk * sin_jkl;
  if( poem < 1e-20 ) poem = 1e-20;

  tel  = SQR( r_ij ) + SQR( r_jk ) + SQR( r_kl ) - SQR( r_li ) -
    2.0 * ( r_ij * r_jk * cos_ijk - r_ij * r_kl * cos_ijk * cos_jkl +
            r_jk * r_kl * cos_jkl );

  arg  = tel / poem;
  if( arg >  1.0 ) arg =  1.0;
  if( arg < -1.0 ) arg = -1.0;

  if( sin_ijk >= 0 && sin_ijk <= MIN_SINE ) sin_ijk = MIN_SINE;
  else if( sin_ijk <= 0 && sin_ijk >= -MIN_SINE ) sin_ijk = -MIN_SINE;
  if( sin_jkl >= 0 && sin_jkl <= MIN_SINE ) sin_jkl = MIN_SINE;
  else if( sin_jkl <= 0 && sin_jkl >= -MIN_SINE ) sin_jkl = -MIN_SINE;

  // dcos_omega_di
  rvec_ScaledSum( dcos_omega_di, (htra-arg*hnra)/r_ij, dvec_ij, -1., dvec_li );
  rvec_ScaledAdd( dcos_omega_di,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_dk );
  rvec_Scale( dcos_omega_di, 2.0 / poem, dcos_omega_di );

  // dcos_omega_dj
  rvec_ScaledSum( dcos_omega_dj,-(htra-arg*hnra)/r_ij, dvec_ij,
                  -htrb / r_jk, dvec_jk );
  rvec_ScaledAdd( dcos_omega_dj,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_dj );
  rvec_ScaledAdd( dcos_omega_dj,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_di );
  rvec_Scale( dcos_omega_dj, 2.0 / poem, dcos_omega_dj );

  // dcos_omega_dk
  rvec_ScaledSum( dcos_omega_dk,-(htrc-arg*hnrc)/r_kl, dvec_kl,
                  htrb / r_jk, dvec_jk );
  rvec_ScaledAdd( dcos_omega_dk,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_di );
  rvec_ScaledAdd( dcos_omega_dk,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_dj );
  rvec_Scale( dcos_omega_dk, 2.0 / poem, dcos_omega_dk );

  // dcos_omega_dl
  rvec_ScaledSum( dcos_omega_dl, (htrc-arg*hnrc)/r_kl, dvec_kl, 1., dvec_li );
  rvec_ScaledAdd( dcos_omega_dl,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_dk );
  rvec_Scale( dcos_omega_dl, 2.0 / poem, dcos_omega_dl );

  return omega;
}


extern SLAVE_FUN(Merge_Torsion_Valence_Angles_ev0_cpe)(merge_param_pack_t *);
extern SLAVE_FUN(Merge_Torsion_Valence_Angles_ev1_cpe)(merge_param_pack_t *);

/////////////version1/////////////
void Torsion_Angles_C(merge_param_pack_t *param) 
{
  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;
  //pack_nbody_atoms *atoms_in    = param->atoms_in;

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  

  int i, j, k, l, pi, pj, pk, pl, pij, plk, natoms;
  int type_i, type_j, type_k, type_l;
  int start_j, end_j;
  int start_pj, end_pj, start_pk, end_pk;

  double Delta_j, Delta_k;
  double r_ij, r_jk, r_kl, r_li;
  double BOA_ij, BOA_jk, BOA_kl;

  double exp_tor2_ij, exp_tor2_jk, exp_tor2_kl;
  double exp_tor1, exp_tor3_DjDk, exp_tor4_DjDk, exp_tor34_inv;
  double exp_cot2_jk, exp_cot2_ij, exp_cot2_kl;
  double fn10, f11_DjDk, dfn11, fn12;
  double theta_ijk, theta_jkl;
  double sin_ijk, sin_jkl;
  double cos_ijk, cos_jkl;
  double tan_ijk_i, tan_jkl_i;
  double omega, cos_omega, cos2omega, cos3omega;
  rvec dcos_omega_di, dcos_omega_dj, dcos_omega_dk, dcos_omega_dl;
  double CV, cmn, CEtors1, CEtors2, CEtors3, CEtors4;
  double CEtors5, CEtors6, CEtors7, CEtors8, CEtors9;
  double Cconj, CEconj1, CEconj2, CEconj3;
  double CEconj4, CEconj5, CEconj6;
  double e_tor, e_con;
  rvec dvec_li;
  rvec force, ext_press;
  ivec rel_box_jl;
  four_body_header *fbh;
  four_body_parameters *fbp;
  bond_data *pbond_ij, *pbond_jk, *pbond_kl;
  bond_order_data *bo_ij, *bo_jk, *bo_kl;
  three_body_interaction_data *p_ijk, *p_jkl;
  double p_tor2 = system->reax_param.gp.l[23];
  double p_tor3 = system->reax_param.gp.l[24];
  double p_tor4 = system->reax_param.gp.l[25];
  double p_cot2 = system->reax_param.gp.l[27];
  reax_list *bonds = (*lists) + BONDS;
  double *Cdbo_list     = bonds->Cdbo_list;
  double *Cdbopi_list   = bonds->Cdbopi_list;
  double *Cdbopi2_list  = bonds->Cdbopi2_list;
  double *BO_list       = bonds->BO_list;
  rvec2  *BOpi_list     = bonds->BOpi_list;

  reax_list *thb_intrs = (*lists) + THREE_BODIES;

  // Virial tallying variables
  double delil[3], deljl[3], delkl[3];
  double eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];

  natoms = system->n;
  

  int tot_torsion = 0;
  for( j = 0; j < natoms; ++j ) 
  {
    type_j = system->packed_atoms[j].type;
    //type_j = atoms_in[j].type;

    //Delta_j = workspace->Delta_boc[j];
    Delta_j = workspace->bo_dboc[j][1];
    //Delta_j = atoms_in[j].Delta_boc;

    start_j = bonds->index[j];
    end_j = bonds->end_index[j];

    for(pk = start_j; pk < end_j; ++pk) 
    {
      pbond_jk = &( bonds->select.bond_list[pk] );
      k = pbond_jk->nbr;
      //bo_jk = &( pbond_jk->bo_data );
      bo_jk = &(bonds->bo_data_list[pk]);
      //BOA_jk = bo_jk->BO - control->thb_cut;
      BOA_jk = BO_list[pk] - control->thb_cut;

      if( system->packed_atoms[j].orig_id > system->packed_atoms[k].orig_id )
	      continue;
      if(system->packed_atoms[j].orig_id == system->packed_atoms[k].orig_id) 
      {
        if (system->packed_atoms[k].x[2] <  system->packed_atoms[j].x[2]) continue;
      	if (system->packed_atoms[k].x[2] == system->packed_atoms[j].x[2] &&
      	    system->packed_atoms[k].x[1] <  system->packed_atoms[j].x[1]) continue;
        if (system->packed_atoms[k].x[2] == system->packed_atoms[j].x[2] &&
      	    system->packed_atoms[k].x[1] == system->packed_atoms[j].x[1] &&
      	    system->packed_atoms[k].x[0] <  system->packed_atoms[j].x[0]) continue;
      }


      //if(atoms_in[j].orig_id > atoms_in[k].orig_id )
	    //  continue;
      //if(atoms_in[j].orig_id == atoms_in[k].orig_id) 
      //{
      //  if (atoms_in[k].x[2] <  atoms_in[j].x[2]) continue;
      //	if (atoms_in[k].x[2] == atoms_in[j].x[2] &&
      //	    atoms_in[k].x[1] <  atoms_in[j].x[1]) continue;
      //  if (atoms_in[k].x[2] == atoms_in[j].x[2] &&
      //	    atoms_in[k].x[1] == atoms_in[j].x[1] &&
      //	    atoms_in[k].x[0] <  atoms_in[j].x[0]) continue;
      //}

      pj = pbond_jk->sym_index; // pj points to j on k's list
     
      if (BOA_jk <= 0)  continue;
      
      type_k  = system->packed_atoms[k].type;
      //type_k  = atoms_in[k].type;

      //Delta_k = workspace->Delta_boc[k];
      Delta_k = workspace->bo_dboc[k][1];
      //Delta_k = atoms_in[k].Delta_boc;

      r_jk    = pbond_jk->d;

      start_pk  = thb_intrs->index[pk]; 
      end_pk    = thb_intrs->end_index[pk]; 
      start_pj  = thb_intrs->index[pj]; 
      end_pj    = thb_intrs->end_index[pj];  

      exp_tor2_jk = exp( -p_tor2 * BOA_jk );
      exp_cot2_jk = exp( -p_cot2 * SQR(BOA_jk - 1.5) );
      exp_tor3_DjDk = exp( -p_tor3 * (Delta_j + Delta_k) );
      exp_tor4_DjDk = exp( p_tor4  * (Delta_j + Delta_k) );
      exp_tor34_inv = 1.0 / (1.0 + exp_tor3_DjDk + exp_tor4_DjDk);
      f11_DjDk = (2.0 + exp_tor3_DjDk) * exp_tor34_inv;

      for( pi = start_pk; pi < end_pk; ++pi ) 
      {
        p_ijk = &( thb_intrs->select.three_body_list[pi] );
        pij = p_ijk->pthb; // pij is pointer to i on j's bond_list
        pbond_ij = &( bonds->select.bond_list[pij] );
        //bo_ij = &( pbond_ij->bo_data );
        bo_ij = &(bonds->bo_data_list[pij]);
        //BOA_ij = bo_ij->BO - control->thb_cut;
        BOA_ij = BO_list[pij] - control->thb_cut;

        //if(BOA_ij <= 0 || (bo_ij->BO*bo_jk->BO <= control->thb_cutsq)) continue;
        if(BOA_ij <= 0 || (BO_list[pij]*BO_list[pk] <= control->thb_cutsq)) continue;
        
        i = p_ijk->thb;
        type_i = system->packed_atoms[i].type;
        //type_i = atoms_in[i].type;
        r_ij = pbond_ij->d;

        theta_ijk = p_ijk->theta;
        sin_ijk = sin( theta_ijk );
        cos_ijk = cos( theta_ijk );
        if( sin_ijk >= 0 && sin_ijk <= MIN_SINE )
          tan_ijk_i = cos_ijk / MIN_SINE;
        else if( sin_ijk <= 0 && sin_ijk >= -MIN_SINE )
          tan_ijk_i = cos_ijk / -MIN_SINE;
        else tan_ijk_i = cos_ijk / sin_ijk;

        exp_tor2_ij = exp( -p_tor2 * BOA_ij );
        exp_cot2_ij = exp( -p_cot2 * SQR(BOA_ij -1.5) );
        
        for( pl = start_pj; pl < end_pj; ++pl ) 
        {

          
          p_jkl = &( thb_intrs->select.three_body_list[pl] );
          l = p_jkl->thb;
          plk = p_jkl->pthb; //pointer to l on k's bond_list!
          pbond_kl = &( bonds->select.bond_list[plk] );
          //bo_kl = &( pbond_kl->bo_data );
          bo_kl = &(bonds->bo_data_list[plk]);

          //BOA_kl = bo_kl->BO - control->thb_cut;
          BOA_kl = BO_list[plk] - control->thb_cut;
          type_l = system->packed_atoms[l].type;
          //type_l = atoms_in[l].type;
          fbh = &(system->reax_param.fbp[type_i][type_j]
                  [type_k][type_l]);
          fbp = &(system->reax_param.fbp[type_i][type_j]
                  [type_k][type_l].prm[0]);
          

          if( i == l || !fbh->cnt || BOA_kl <= 0 ||
              BO_list[pij] * BO_list[pk] * BO_list[plk] <= control->thb_cut) continue;
              //bo_ij->BO * bo_jk->BO * bo_kl->BO <= control->thb_cut) continue;

          tot_torsion ++;

          r_kl = pbond_kl->d;
          //BOA_kl = bo_kl->BO - control->thb_cut;
          BOA_kl = BO_list[plk] - control->thb_cut;

          theta_jkl = p_jkl->theta;
          sin_jkl = sin( theta_jkl );
          cos_jkl = cos( theta_jkl );
          if( sin_jkl >= 0 && sin_jkl <= MIN_SINE )
            tan_jkl_i = cos_jkl / MIN_SINE;
          else if( sin_jkl <= 0 && sin_jkl >= -MIN_SINE )
            tan_jkl_i = cos_jkl / -MIN_SINE;
          else tan_jkl_i = cos_jkl /sin_jkl;
          dvec_li[0] = system->packed_atoms[i].x[0] - system->packed_atoms[l].x[0];
          dvec_li[1] = system->packed_atoms[i].x[1] - system->packed_atoms[l].x[1];
          dvec_li[2] = system->packed_atoms[i].x[2] - system->packed_atoms[l].x[2];

          //dvec_li[0] = atoms_in[i].x[0] - system->my_atoms[l].x[0];
          //dvec_li[1] = atoms_in[i].x[1] - system->my_atoms[l].x[1];
          //dvec_li[2] = atoms_in[i].x[2] - system->my_atoms[l].x[2];
          r_li = sqrt( SQR(dvec_li[0]) + SQR(dvec_li[1]) + SQR(dvec_li[2]));

          /* omega and its derivative */
          omega = Calculate_Omega( pbond_ij->dvec, r_ij,
                                   pbond_jk->dvec, r_jk,
                                   pbond_kl->dvec, r_kl,
                                   dvec_li, r_li,
                                   p_ijk, p_jkl,
                                   dcos_omega_di, dcos_omega_dj,
                                   dcos_omega_dk, dcos_omega_dl,
                                   out_control );

          cos_omega = cos( omega );
          cos2omega = cos( 2. * omega );
          cos3omega = cos( 3. * omega );
          /* end omega calculations */

          /* torsion energy */
          //exp_tor1 = exp(fbp->p_tor1 * SQR(2.0 - bo_jk->BO_pi - f11_DjDk));
          exp_tor1 = exp(fbp->p_tor1 * SQR(2.0 - BOpi_list[pk][0] - f11_DjDk));
          exp_tor2_kl = exp( -p_tor2 * BOA_kl );
          exp_cot2_kl = exp( -p_cot2 * SQR(BOA_kl - 1.5) );
          fn10 = (1.0-exp_tor2_ij) * (1.0-exp_tor2_jk) * (1.0-exp_tor2_kl);

          CV = 0.5 * ( fbp->V1 * (1.0 + cos_omega) +
                       fbp->V2 * exp_tor1 * (1.0 - cos2omega) +
                       fbp->V3 * (1.0 + cos3omega) );

          data->my_en.e_tor += e_tor = fn10 * sin_ijk * sin_jkl * CV;

          dfn11 = (-p_tor3 * exp_tor3_DjDk +
                   (p_tor3 * exp_tor3_DjDk - p_tor4 * exp_tor4_DjDk) *
                   (2.0+exp_tor3_DjDk) * exp_tor34_inv) * exp_tor34_inv;

          CEtors1 = sin_ijk * sin_jkl * CV;

          //CEtors2 = -fn10 * 2.0 * fbp->p_tor1 * fbp->V2 * exp_tor1 *
          //          (2.0 - bo_jk->BO_pi - f11_DjDk) * (1.0-SQR(cos_omega)) *
          //          sin_ijk * sin_jkl;
          CEtors2 = -fn10 * 2.0 * fbp->p_tor1 * fbp->V2 * exp_tor1 *
                          (2.0 - BOpi_list[pk][0] -f11_DjDk) * (1.0 -SQR(cos_omega)) *
                          sin_ijk * sin_jkl;

          CEtors3 = CEtors2 * dfn11;

          CEtors4 = CEtors1 * p_tor2 * exp_tor2_ij *
                    (1.0-exp_tor2_jk) * (1.0-exp_tor2_kl);
          CEtors5 = CEtors1 * p_tor2 *
                    (1.0-exp_tor2_ij) * exp_tor2_jk * (1.0-exp_tor2_kl);
          CEtors6 = CEtors1 * p_tor2 *
                    (1.0-exp_tor2_ij) * (1.0-exp_tor2_jk) * exp_tor2_kl;

          cmn = -fn10 * CV;
          CEtors7 = cmn * sin_jkl * tan_ijk_i;
          CEtors8 = cmn * sin_ijk * tan_jkl_i;

          CEtors9 = fn10 * sin_ijk * sin_jkl *
            (0.5 * fbp->V1 - 2.0 * fbp->V2 * exp_tor1 * cos_omega +
             1.5 * fbp->V3 * (cos2omega + 2.0 * SQR(cos_omega)));
          /* end  of torsion energy */

          /* 4-body conjugation energy */
          fn12 = exp_cot2_ij * exp_cot2_jk * exp_cot2_kl;
          data->my_en.e_con += e_con = 
                            fbp->p_cot1 * fn12 *
                            (1.0+(SQR(cos_omega)-1.0) * sin_ijk * sin_jkl);

          Cconj = -2.0 * fn12 * fbp->p_cot1 * p_cot2 *
                  (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);

          CEconj1 = Cconj * (BOA_ij - 1.5e0);
          CEconj2 = Cconj * (BOA_jk - 1.5e0);
          CEconj3 = Cconj * (BOA_kl - 1.5e0);

          CEconj4 = -fbp->p_cot1 * fn12 *
                    (SQR(cos_omega) - 1.0) * sin_jkl * tan_ijk_i;
          CEconj5 = -fbp->p_cot1 * fn12 *
                    (SQR(cos_omega) - 1.0) * sin_ijk * tan_jkl_i;
          CEconj6 = 2.0 * fbp->p_cot1 * fn12 * cos_omega * sin_ijk * sin_jkl;
          /* end 4-body conjugation energy */

          /* forces */
          //bo_jk->Cdbopi += CEtors2;
          Cdbopi_list[pk] += CEtors2;

          workspace->fCdDelta[j][3] += CEtors3;
          workspace->fCdDelta[k][3] += CEtors3;
          //bo_ij->Cdbo += (CEtors4 + CEconj1);
          //bo_jk->Cdbo += (CEtors5 + CEconj2);
          //bo_kl->Cdbo += (CEtors6 + CEconj3);
          Cdbo_list[pij] += (CEtors4 + CEconj1);
          Cdbo_list[pk] += (CEtors5 + CEconj2);
          Cdbo_list[plk] += (CEtors6 + CEconj3);


          double c74 = CEtors7+CEconj4;
          double c85 = CEtors8+CEconj5;
          double c96 = CEtors9+CEconj6;
          
          ////////////// control->virial == 0 ////////////
          /* dcos_theta_ijk */
          workspace->fCdDelta[i][0] += c74 * p_ijk->dcos_dk[0];
          workspace->fCdDelta[i][1] += c74 * p_ijk->dcos_dk[1];
          workspace->fCdDelta[i][2] += c74 * p_ijk->dcos_dk[2];

          workspace->fCdDelta[j][0] += c74 * p_ijk->dcos_dj[0];
          workspace->fCdDelta[j][1] += c74 * p_ijk->dcos_dj[1];
          workspace->fCdDelta[j][2] += c74 * p_ijk->dcos_dj[2];
          
          workspace->fCdDelta[k][0] += c74 * p_ijk->dcos_di[0];
          workspace->fCdDelta[k][1] += c74 * p_ijk->dcos_di[1];
          workspace->fCdDelta[k][2] += c74 * p_ijk->dcos_di[2];

          /* dcos_theta_jkl */
          workspace->fCdDelta[j][0] += c85 * p_jkl->dcos_di[0];
          workspace->fCdDelta[j][1] += c85 * p_jkl->dcos_di[1];
          workspace->fCdDelta[j][2] += c85 * p_jkl->dcos_di[2];

          workspace->fCdDelta[k][0] += c85 * p_jkl->dcos_dj[0];
          workspace->fCdDelta[k][1] += c85 * p_jkl->dcos_dj[1];
          workspace->fCdDelta[k][2] += c85 * p_jkl->dcos_dj[2];
          
          workspace->fCdDelta[l][0] += c85 * p_jkl->dcos_dk[0];
          workspace->fCdDelta[l][1] += c85 * p_jkl->dcos_dk[1];
          workspace->fCdDelta[l][2] += c85 * p_jkl->dcos_dk[2];

          /* dcos_omega */
          workspace->fCdDelta[i][0] += c96 * dcos_omega_di[0];
          workspace->fCdDelta[i][1] += c96 * dcos_omega_di[1];
          workspace->fCdDelta[i][2] += c96 * dcos_omega_di[2];

          workspace->fCdDelta[j][0] += c96 * dcos_omega_dj[0];
          workspace->fCdDelta[j][1] += c96 * dcos_omega_dj[1];
          workspace->fCdDelta[j][2] += c96 * dcos_omega_dj[2];
          
          workspace->fCdDelta[k][0] += c96 * dcos_omega_dk[0];
          workspace->fCdDelta[k][1] += c96 * dcos_omega_dk[1];
          workspace->fCdDelta[k][2] += c96 * dcos_omega_dk[2];
          
          workspace->fCdDelta[l][0] += c96 * dcos_omega_dl[0];
          workspace->fCdDelta[l][1] += c96 * dcos_omega_dl[1];
          workspace->fCdDelta[l][2] += c96 * dcos_omega_dl[2];
          
          

          /* tally into per-atom virials */
          if( system->vflag_atom || system->evflag) 
          {
            // acquire vectors
            delil[0] = system->packed_atoms[l].x[0] - system->packed_atoms[i].x[0];
            delil[1] = system->packed_atoms[l].x[1] - system->packed_atoms[i].x[1];
            delil[2] = system->packed_atoms[l].x[2] - system->packed_atoms[i].x[2];
            deljl[0] = system->packed_atoms[l].x[0] - system->packed_atoms[j].x[0];
            deljl[1] = system->packed_atoms[l].x[1] - system->packed_atoms[j].x[1];
            deljl[2] = system->packed_atoms[l].x[2] - system->packed_atoms[j].x[2];
            delkl[0] = system->packed_atoms[l].x[0] - system->packed_atoms[k].x[0];
            delkl[1] = system->packed_atoms[l].x[1] - system->packed_atoms[k].x[1];
            delkl[2] = system->packed_atoms[l].x[2] - system->packed_atoms[k].x[2];

            // dcos_theta_ijk
            fi_tmp[0] = c74 * p_ijk->dcos_dk[0];
            fi_tmp[1] = c74 * p_ijk->dcos_dk[1];
            fi_tmp[2] = c74 * p_ijk->dcos_dk[2];
            fj_tmp[0] = c74 * p_ijk->dcos_dj[0];
            fj_tmp[1] = c74 * p_ijk->dcos_dj[1];
            fj_tmp[2] = c74 * p_ijk->dcos_dj[2];
            fk_tmp[0] = c74 * p_ijk->dcos_di[0];
            fk_tmp[1] = c74 * p_ijk->dcos_di[1];
            fk_tmp[2] = c74 * p_ijk->dcos_di[2];

            // dcos_theta_jkl
            fj_tmp[0] += c85 * p_jkl->dcos_di[0];
            fj_tmp[1] += c85 * p_jkl->dcos_di[1];
            fj_tmp[2] += c85 * p_jkl->dcos_di[2];
            fk_tmp[0] += c85 * p_jkl->dcos_dj[0];
            fk_tmp[1] += c85 * p_jkl->dcos_dj[1];
            fk_tmp[2] += c85 * p_jkl->dcos_dj[2];


            // dcos_omega
            fi_tmp[0] += c96 * dcos_omega_di[0];
            fi_tmp[1] += c96 * dcos_omega_di[1];
            fi_tmp[2] += c96 * dcos_omega_di[2];
            fj_tmp[0] += c96 * dcos_omega_dj[0];
            fj_tmp[1] += c96 * dcos_omega_dj[1];
            fj_tmp[2] += c96 * dcos_omega_dj[2];
            fk_tmp[0] += c96 * dcos_omega_dk[0];
            fk_tmp[1] += c96 * dcos_omega_dk[1];
            fk_tmp[2] += c96 * dcos_omega_dk[2];

            // tally
            eng_tmp = e_tor + e_con;
            if(system->evflag)
            {
              if (system->eflag_global) 
              {
                  system->eng_vdwl += eng_tmp;
              }
              if (system->eflag_atom) 
              {
                system->eatom[j] += 0.5 * eng_tmp;//evdwl;
                system->eatom[k] += 0.5 * eng_tmp;//evdwl;
              }
            }//if
            if(system->vflag_atom)
            {
              double v[6];
              double *fi = fi_tmp, *fj = fj_tmp, *fk = fk_tmp;
              double *dril = delil, *drjl = deljl, *drkl = delkl; 
              v[0] = 0.25*(dril[0]*fi[0] + drjl[0]*fj[0] + drkl[0]*fk[0]);
              v[1] = 0.25*(dril[1]*fi[1] + drjl[1]*fj[1] + drkl[1]*fk[1]);
              v[2] = 0.25*(dril[2]*fi[2] + drjl[2]*fj[2] + drkl[2]*fk[2]);
              v[3] = 0.25*(dril[0]*fi[1] + drjl[0]*fj[1] + drkl[0]*fk[1]);
              v[4] = 0.25*(dril[0]*fi[2] + drjl[0]*fj[2] + drkl[0]*fk[2]);
              v[5] = 0.25*(dril[1]*fi[2] + drjl[1]*fj[2] + drkl[1]*fk[2]);

              system->vatom[i][0] += v[0]; 
              system->vatom[i][1] += v[1]; 
              system->vatom[i][2] += v[2];
              system->vatom[i][3] += v[3]; 
              system->vatom[i][4] += v[4]; 
              system->vatom[i][5] += v[5];

              system->vatom[j][0] += v[0]; 
              system->vatom[j][1] += v[1]; 
              system->vatom[j][2] += v[2];
              system->vatom[j][3] += v[3]; 
              system->vatom[j][4] += v[4]; 
              system->vatom[j][5] += v[5];

              system->vatom[k][0] += v[0]; 
              system->vatom[k][1] += v[1]; 
              system->vatom[k][2] += v[2];
              system->vatom[k][3] += v[3]; 
              system->vatom[k][4] += v[4]; 
              system->vatom[k][5] += v[5];

              system->vatom[l][0] += v[0]; 
              system->vatom[l][1] += v[1]; 
              system->vatom[l][2] += v[2];
              system->vatom[l][3] += v[3]; 
              system->vatom[l][4] += v[4]; 
              system->vatom[l][5] += v[5];
            }
          }
        } // pl loop ends
      } // pi loop ends
    } // pk loop ends
  } // j loop
}


//swcache_lock_t *locks = NULL;
//swcache_lock_t *locks_frc = NULL;
void Merge_Torsion_Valence_Angles(merge_param_pack_t *param)
{
  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;
 
  int i, j, k, l, pi, pj, pk, pl, pij, plk, natoms, pw, w, h, ph, t_tmp, t;
  int type_i, type_j, type_k, type_l, type_w, type_h;
  int start_j, end_j, start_k, end_k;
  int start_pj, end_pj, start_pk, end_pk;

  double Delta_j, Delta_k;
  double r_ij, r_jk, r_kl, r_li;
  double BOA_ij, BOA_jk, BOA_kl, BOA_kw, BOA_hj;

  double exp_tor2_ij, exp_tor2_jk, exp_tor2_kl;
  double exp_tor1, exp_tor3_DjDk, exp_tor4_DjDk, exp_tor34_inv;
  double exp_cot2_jk, exp_cot2_ij, exp_cot2_kl;
  double fn10, f11_DjDk, dfn11, fn12;
  double theta_ijk, theta_jkl, theta_jkw, theta_hjk;
  double cos_theta_jkw, cos_theta_hjk;
  double sin_ijk, sin_jkl;
  double cos_ijk, cos_jkl;
  double tan_ijk_i, tan_jkl_i;
  double omega, cos_omega, cos2omega, cos3omega;
  rvec dcos_omega_di, dcos_omega_dj, dcos_omega_dk, dcos_omega_dl;
  double CV, cmn, CEtors1, CEtors2, CEtors3, CEtors4;
  double CEtors5, CEtors6, CEtors7, CEtors8, CEtors9;
  double Cconj, CEconj1, CEconj2, CEconj3;
  double CEconj4, CEconj5, CEconj6;
  double e_tor, e_con;
  rvec dvec_li;
  rvec force, ext_press;
  ivec rel_box_jl;
  four_body_header *fbh;
  four_body_parameters *fbp;
  bond_data *pbond_ij, *pbond_jk, *pbond_kl, *pbond_kw, *pbond_hj, *pbond_pj, *pbond_jt;
  bond_order_data *bo_ij, *bo_jk, *bo_kl, *bo_kw, *bo_hj, *bo_jt;
  three_body_interaction_data *p_ijk, *p_jkl, *p_jkw, *p_tmp, *p_hjk;
  three_body_interaction_data thb_intrs_data;
  three_body_interaction_data thb_intrs_data_hjk;
  
    
  reax_list *bonds = (*lists) + BONDS;
  double *Cdbo_list     = bonds->Cdbo_list;
  double *Cdbopi_list   = bonds->Cdbopi_list;
  double *Cdbopi2_list  = bonds->Cdbopi2_list;
  double *BO_list       = bonds->BO_list;
  rvec2  *BOpi_list     = bonds->BOpi_list;

  // Virial tallying variables
  double delil[3], deljl[3], delkl[3];
  double eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];

  //valence
  int cnt;
  double temp, temp_bo_jt, pBOjt7;
  double p_val1, p_val2, p_val3, p_val4, p_val5;
  double p_val6, p_val7, p_val8, p_val9, p_val10;
  double p_pen1, p_pen2, p_pen3, p_pen4;
  double p_coa1, p_coa2, p_coa3, p_coa4;
  double trm8, expval6, expval7, expval2theta, expval12theta, exp3hj, exp3jk;
  double exp_pen2hj, exp_pen2jk, exp_pen3, exp_pen4, trm_pen34, exp_coa2;
  double dSBO1, dSBO2, SBO, SBO2, CSBO2, SBOp, prod_SBO, vlpadj;
  double CEval1, CEval2, CEval3, CEval4, CEval5, CEval6, CEval7, CEval8;
  double CEpen1, CEpen2, CEpen3;
  double e_ang, e_coa, e_pen;
  double CEcoa1, CEcoa2, CEcoa3, CEcoa4, CEcoa5;
  double Cf7hj, Cf7jk, Cf8j, Cf9j;
  double f7_hj, f7_jk, f8_Dj, f9_Dj;
  double Ctheta_0, theta_0, theta_00, theta, cos_theta, sin_theta_hjk;
  double fh_tmp[3];
  double deljk[3], delhj[3];
  three_body_header *thbh;
  three_body_parameters *thbp;
  double p_tor2, p_tor3, p_tor4, p_cot2;
  
  #define SUNWAY
  #ifdef SUNWAY
  int nall = bonds->num_intrs;
  //if(locks == NULL)
  //  locks = swcache_u_prepare_locks(nall, 3);
  //if(locks_frc == NULL)
  //  locks_frc = swcache_u_prepare_locks(nall, 3);
  
  merge_param_pack_t pm;
  //pm.locks        = locks;
  pm.locks        = bonds->locks_cdbo;
  //pm.locks_frc    = locks_frc;
  pm.locks_frc    = bonds->locks_frc;

  pm.index        = bonds->index;
  pm.end_index    = bonds->end_index;

  pm.system       = system;
  pm.control      = control;
  pm.data         = data;
  pm.workspace    = workspace;
  pm.lists        = lists;
  //pm.out_control  = out_control;
  pm.packed_atoms   = system->packed_atoms;
  pm.bond_list      = bonds->select.bond_list;
  pm.Cdbo_list      = bonds->Cdbo_list;
  pm.Cdbopi_list    = bonds->Cdbopi_list;
  pm.Cdbopi2_list   = bonds->Cdbopi2_list;
  pm.BO_list        = bonds->BO_list;
  pm.BOpi_list      = bonds->BOpi_list;
  pm.fCdDelta       = workspace->fCdDelta;
  pm.bo_dboc        = workspace->bo_dboc;
  pm.vlpex          = workspace->vlpex;
  pm.nlp            = workspace->nlp;
  pm.dDelta_lp      = workspace->dDelta_lp;
  pm.Delta          = workspace->Delta;
  pm.Delta_val      = workspace->Delta_val;
  pm.bo_dboc        = workspace->bo_dboc;


  pm.p_tor2   = system->reax_param.gp.l[23];
  pm.p_tor3   = system->reax_param.gp.l[24];
  pm.p_tor4   = system->reax_param.gp.l[25];
  pm.p_cot2   = system->reax_param.gp.l[27];
  pm.p_val6   = system->reax_param.gp.l[14];
  pm.p_val8   = system->reax_param.gp.l[33];
  pm.p_val9   = system->reax_param.gp.l[16];
  pm.p_val10  = system->reax_param.gp.l[17];
  pm.p_pen2   = system->reax_param.gp.l[19];
  pm.p_pen3   = system->reax_param.gp.l[20];
  pm.p_pen4   = system->reax_param.gp.l[21];
  pm.p_coa2   = system->reax_param.gp.l[2];
  pm.p_coa3   = system->reax_param.gp.l[38];
  pm.p_coa4   = system->reax_param.gp.l[30];
  pm.n          = system->n;
  pm.thb_cut    = control->thb_cut;
  pm.thb_cutsq  = control->thb_cutsq;
  pm.ntypes     = system->reax_param.num_atom_types;

  double packed_eng[8];
  for(i = 0; i < 8; i++)
    packed_eng[i] = 0;
  pm.packed_eng = packed_eng;

  //perf_config_t conf;
  //conf.pcrc = PCRC_ALL;
  //conf.pcr0 = PC0_CYCLE;
  //conf.pcr1 = PC1_CYCLE;
  //conf.pcr2 = PC2_N_GLD;
  ////conf.pcr2 = PC2_N_GF_AND_A;
  ////conf.pcr2 = PC2_N_DMA_REQ;
  //lwpf_init(&conf);
  //int myrank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  //NT = 4(ntypes = 4);
  int fbp_cnt[256]; 
  four_body_parameters fbp_prm[256];
  int thbp_cnt[128];
  three_body_parameters thbp_prm[100];
  
  int psum = 0, thsum = 0;
  int ntypes = system->reax_param.num_atom_types;
  for(i = 0; i < ntypes; i++)
  {
    pm.sbp[i].p_val3 = system->reax_param.sbp[i].p_val3;
    pm.sbp[i].p_val5 = system->reax_param.sbp[i].p_val5;
    for(j = 0; j < ntypes; j++)
    {
      for(k = 0; k < ntypes; k++)
      {
        thbh = &(system->reax_param.thbp[i][j][k]); 
        thbp_cnt[NT3+i*ntypes*ntypes+j*ntypes+k] = thbh->cnt;//count
        thbp_cnt[i*ntypes*ntypes+j*ntypes+k]    = thsum;//index

        for(t = 0; t < thbh->cnt; t++)
        {
          thbp_prm[thsum] = thbh->prm[t];
          thsum++;
        }
        for(l = 0; l < ntypes; l++)
        {
          
          fbh = &(system->reax_param.fbp[i][j][k][l]);
          fbp_prm[psum].p_tor1 = fbh->prm[0].p_tor1;
          fbp_prm[psum].p_cot1 = fbh->prm[0].p_cot1;
          fbp_prm[psum].V1     = fbh->prm[0].V1;
          fbp_prm[psum].V2     = fbh->prm[0].V2;
          fbp_prm[psum].V3     = fbh->prm[0].V3;
          fbp_cnt[psum]        = fbh->cnt;
          psum++;
        }
      }
    }
  }
  pm.fbp_prm  = fbp_prm;
  pm.fbp_cnt  = fbp_cnt;
  pm.fbp_tot  = psum;

  pm.thbp_prm = thbp_prm;
  pm.thbp_cnt = thbp_cnt;
  pm.thbp_tot = thsum;


  if(athread_idle() == 0)
    athread_init();
  if (system->evflag)
    athread_spawn(Merge_Torsion_Valence_Angles_ev1_cpe, &pm);
  else
    athread_spawn(Merge_Torsion_Valence_Angles_ev0_cpe, &pm);
  athread_join();
  
  //if(myrank == 0)
  //{
  //  lwpf_report_summary(stdout, &conf);
  //}

  system->eng_vdwl  += pm.packed_eng[0];
  data->my_en.e_ang += pm.packed_eng[1];
  data->my_en.e_pen += pm.packed_eng[2];
  data->my_en.e_coa += pm.packed_eng[3];
  data->my_en.e_tor += pm.packed_eng[4];
  data->my_en.e_con += pm.packed_eng[5];
 
  //if(locks != NULL)
  //{
  //  free(locks); 
  //  locks = NULL;
  //}
  //if(locks_frc != NULL)
  //{
  //  free(locks_frc); 
  //  locks_frc = NULL;
  //}
  return;
  
#else
    
  p_tor2 = system->reax_param.gp.l[23];
  p_tor3 = system->reax_param.gp.l[24];
  p_tor4 = system->reax_param.gp.l[25];
  p_cot2 = system->reax_param.gp.l[27];
  p_val6 = system->reax_param.gp.l[14];
  p_val8 = system->reax_param.gp.l[33];
  p_val9 = system->reax_param.gp.l[16];
  p_val10 = system->reax_param.gp.l[17];
  p_pen2 = system->reax_param.gp.l[19];
  p_pen3 = system->reax_param.gp.l[20];
  p_pen4 = system->reax_param.gp.l[21];
  p_coa2 = system->reax_param.gp.l[2];
  p_coa3 = system->reax_param.gp.l[38];
  p_coa4 = system->reax_param.gp.l[30];

  natoms = system->n;

  for( j = 0; j < system->N; ++j ) 
  {
    type_j = system->packed_atoms[j].type;
    //Delta_j = workspace->Delta_boc[j];
    Delta_j = workspace->bo_dboc[j][1];
    start_j = Start_Index(j, bonds);
    end_j = End_Index(j, bonds);

    if (type_j < 0 || end_j <= start_j) continue;
    p_val3 = system->reax_param.sbp[ type_j ].p_val3;
    p_val5 = system->reax_param.sbp[ type_j ].p_val5;

    SBOp = 0, prod_SBO = 1;
    for( t = start_j; t < end_j; ++t ) 
    {
      //bo_jt = &(bonds->select.bond_list[t].bo_data);
      bo_jt = &(bonds->bo_data_list[t]);
      //SBOp += (bo_jt->BO_pi + bo_jt->BO_pi2);
      SBOp += (BOpi_list[t][0] + BOpi_list[t][1]);
      //temp = SQR( bo_jt->BO );
      temp = SQR(BO_list[t]);
      temp *= temp;
      temp *= temp;
      prod_SBO *= exp( -temp );
    }

    if (workspace->vlpex[j] >= 0) 
    {
      vlpadj = 0;
      dSBO2 = prod_SBO - 1;
    } 
    else 
    {
      vlpadj = workspace->nlp[j];
      dSBO2 = (prod_SBO - 1) * (1 - p_val8 * workspace->dDelta_lp[j]);
    }

    //SBO = SBOp + (1 - prod_SBO) * (-workspace->Delta_boc[j] - p_val8 * vlpadj);
    //dSBO1 = -8 * prod_SBO * ( workspace->Delta_boc[j] + p_val8 * vlpadj );
    SBO = SBOp + (1 - prod_SBO) * (-workspace->bo_dboc[j][1] - p_val8 * vlpadj);
    dSBO1 = -8 * prod_SBO * ( workspace->bo_dboc[j][1] + p_val8 * vlpadj );

    if (SBO <= 0)
      SBO2 = 0, CSBO2 = 0;
    else if (SBO > 0 && SBO <= 1) 
    {
        SBO2 = pow( SBO, p_val9 );
        CSBO2 = p_val9 * pow( SBO, p_val9 - 1 );
    }
    else if (SBO > 1 && SBO < 2) 
    {
      SBO2 = 2 - pow( 2-SBO, p_val9 );
      CSBO2 = p_val9 * pow( 2 - SBO, p_val9 - 1 );
    }
    else
      SBO2 = 2, CSBO2 = 0;

    //expval6 = exp( p_val6 * workspace->Delta_boc[j] );
    expval6 = exp( p_val6 * workspace->bo_dboc[j][1] );

    for( pk = start_j; pk < end_j; ++pk ) 
    {
      pbond_jk = &( bonds->select.bond_list[pk] );
      k = pbond_jk->nbr;
      //bo_jk = &( pbond_jk->bo_data );
      bo_jk = &(bonds->bo_data_list[pk]);
      //BOA_jk = bo_jk->BO - control->thb_cut;
      BOA_jk = BO_list[pk] - control->thb_cut;

      if(BOA_jk > 0.0 && (j < system->n || pbond_jk->nbr < system->n)) 
      {
        pj = pbond_jk->sym_index; // pj points to j on k's list
        pbond_pj = &(bonds->select.bond_list[pj]);
        start_k = Start_Index(k, bonds);
        end_k = End_Index(k, bonds);
        if(end_k > start_k)
        //if (Num_Entries(pj, thb_intrs)) 
        {
          type_k = system->packed_atoms[k].type;
          //Delta_k = workspace->Delta_boc[k];
          Delta_k = workspace->bo_dboc[k][1];
          r_jk = pbond_jk->d;


          exp_tor2_jk = exp( -p_tor2 * BOA_jk );
          exp_cot2_jk = exp( -p_cot2 * SQR(BOA_jk - 1.5) );
          exp_tor3_DjDk = exp( -p_tor3 * (Delta_j + Delta_k) );
          exp_tor4_DjDk = exp( p_tor4  * (Delta_j + Delta_k) );
          exp_tor34_inv = 1.0 / (1.0 + exp_tor3_DjDk + exp_tor4_DjDk);
          f11_DjDk = (2.0 + exp_tor3_DjDk) * exp_tor34_inv;

          for(ph = start_j; ph < end_j; ++ph) 
          {
            pbond_hj = &(bonds->select.bond_list[ph]);
            //bo_hj = &(pbond_hj->bo_data);
            bo_hj = &(bonds->bo_data_list[ph]);
            //BOA_hj = bo_hj->BO - control->thb_cut;
            BOA_hj = BO_list[ph] - control->thb_cut;
            h = pbond_hj->nbr;

            if(ph != pk)
            {
              type_h = system->packed_atoms[h].type;
              p_hjk   = &thb_intrs_data_hjk;
              
              Calculate_Theta( pbond_jk->dvec, pbond_jk->d,
                               pbond_hj->dvec, pbond_hj->d,
                               &theta_hjk, &cos_theta_hjk);

              Calculate_dCos_Theta( pbond_jk->dvec, pbond_jk->d,
                                    pbond_hj->dvec, pbond_hj->d,
                                    &(p_hjk->dcos_di), &(p_hjk->dcos_dj),
                                    &(p_hjk->dcos_dk));
                  
              p_hjk->thb = h;
              p_hjk->pthb = ph;
              p_hjk->theta = theta_hjk;
              
              sin_theta_hjk = sin(theta_hjk);
              if (sin_theta_hjk < 1.0e-5)
                sin_theta_hjk = 1.0e-5;


              i = h;
              type_i= type_h;
              theta_hjk = p_hjk->theta;
              pbond_ij = &(bonds->select.bond_list[ph]);
              //bo_ij = &(pbond_ij->bo_data);
              bo_ij = &(bonds->bo_data_list[ph]);
              p_ijk = p_hjk;
             
              //printf("valence:(j, k, h, i) = (%d, %d, %d, %d)\n", j, k, h, i);

              //Valence
              //if((ph > pk) && (j < system->n) &&
              //  (bo_jk->BO > control->thb_cut) &&
              //  (bo_hj->BO > control->thb_cut) &&
              //  (bo_jk->BO * bo_hj->BO > control->thb_cutsq) ) 
              if((ph > pk) && (j < system->n) &&
                  (BO_list[pk] > control->thb_cut) &&
                  (BO_list[ph] > control->thb_cut) &&
                  (BO_list[pk] * BO_list[ph] > control->thb_cutsq) ) 

              {
                thbh = &( system->reax_param.thbp[ type_k][ type_j][ type_h]);

                for( cnt = 0; cnt < thbh->cnt; ++cnt ) 
                {
                  if (fabs(thbh->prm[cnt].p_val1) > 0.001) 
                  {
                    thbp = &( thbh->prm[cnt] );

                    /* ANGLE ENERGY */
                    p_val1 = thbp->p_val1;
                    p_val2 = thbp->p_val2;
                    p_val4 = thbp->p_val4;
                    p_val7 = thbp->p_val7;
                    theta_00 = thbp->theta_00;

                    exp3jk = exp( -p_val3 * pow( BOA_jk, p_val4 ) );
                    f7_jk = 1.0 - exp3jk;
                    Cf7jk = p_val3 * p_val4 * pow( BOA_jk, p_val4 - 1.0 ) * exp3jk;

                    exp3hj = exp( -p_val3 * pow( BOA_hj, p_val4 ) );
                    f7_hj = 1.0 - exp3hj;
                    Cf7hj = p_val3 * p_val4 * pow( BOA_hj, p_val4 - 1.0 ) * exp3hj;

                    //expval7 = exp( -p_val7 * workspace->Delta_boc[j] );
                    expval7 = exp( -p_val7 * workspace->bo_dboc[j][1]);
                    trm8 = 1.0 + expval6 + expval7;
                    f8_Dj = p_val5 - ( (p_val5 - 1.0) * (2.0 + expval6) / trm8 );
                    Cf8j = ( (1.0 - p_val5) / SQR(trm8) ) *
                      ( p_val6 * expval6 * trm8 -
                        (2.0 + expval6) * ( p_val6*expval6 - p_val7*expval7 ) );

                    theta_0 = 180.0 - theta_00 * (1.0 -
                                                  exp(-p_val10 * (2.0 - SBO2)));
                    theta_0 = DEG2RAD( theta_0 );

                    expval2theta  = exp( -p_val2 * SQR(theta_0 - theta_hjk) );
                    if (p_val1 >= 0)
                      expval12theta = p_val1 * (1.0 - expval2theta);
                    else // To avoid linear Me-H-Me angles (6/6/06)
                      expval12theta = p_val1 * -expval2theta;

                    CEval1 = Cf7jk * f7_hj * f8_Dj * expval12theta;
                    CEval2 = Cf7hj * f7_jk * f8_Dj * expval12theta;
                    CEval3 = Cf8j  * f7_hj * f7_jk * expval12theta;
                    CEval4 = -2.0 * p_val1 * p_val2 * f7_jk * f7_hj * f8_Dj *
                      expval2theta * (theta_0 - theta_hjk);

                    Ctheta_0 = p_val10 * DEG2RAD(theta_00) *
                      exp( -p_val10 * (2.0 - SBO2) );

                    CEval5 = -CEval4 * Ctheta_0 * CSBO2;
                    CEval6 = CEval5 * dSBO1;
                    CEval7 = CEval5 * dSBO2;
                    CEval8 = -CEval4 / sin_theta_hjk;

                    data->my_en.e_ang += e_ang =
                      f7_jk * f7_hj * f8_Dj * expval12theta;
                    /* END ANGLE ENERGY*/

                    /* PENALTY ENERGY */
                    p_pen1 = thbp->p_pen1;
                    
                    exp_pen2jk = exp( -p_pen2 * SQR( BOA_jk - 2.0 ) );
                    exp_pen2hj = exp( -p_pen2 * SQR( BOA_hj - 2.0 ) );
                    exp_pen3 = exp( -p_pen3 * workspace->Delta[j] );
                    exp_pen4 = exp(  p_pen4 * workspace->Delta[j] );
                    trm_pen34 = 1.0 + exp_pen3 + exp_pen4;
                    f9_Dj = ( 2.0 + exp_pen3 ) / trm_pen34;
                    Cf9j = ( -p_pen3 * exp_pen3 * trm_pen34 -
                             (2.0 + exp_pen3) * ( -p_pen3 * exp_pen3 +
                                                  p_pen4 * exp_pen4 ) ) /
                      SQR( trm_pen34 );

                    data->my_en.e_pen += e_pen =
                      p_pen1 * f9_Dj * exp_pen2jk * exp_pen2hj;

                    CEpen1 = e_pen * Cf9j / f9_Dj;
                    temp   = -2.0 * p_pen2 * e_pen;
                    CEpen2 = temp * (BOA_jk - 2.0);
                    CEpen3 = temp * (BOA_hj - 2.0);
                    /* END PENALTY ENERGY */

                    /* COALITION ENERGY */
                    p_coa1 = thbp->p_coa1;

                    exp_coa2 = exp( p_coa2 * workspace->Delta_val[j] );
                    data->my_en.e_coa += e_coa =
                      p_coa1 / (1. + exp_coa2) *
                      exp( -p_coa3 * SQR(workspace->bo_dboc[k][0]-BOA_jk) ) *
                      exp( -p_coa3 * SQR(workspace->bo_dboc[h][0]-BOA_hj) ) *
                      exp( -p_coa4 * SQR(BOA_jk - 1.5) ) *
                      exp( -p_coa4 * SQR(BOA_hj - 1.5) );

                    CEcoa1 = -2 * p_coa4 * (BOA_jk - 1.5) * e_coa;
                    CEcoa2 = -2 * p_coa4 * (BOA_hj - 1.5) * e_coa;
                    CEcoa3 = -p_coa2 * exp_coa2 * e_coa / (1 + exp_coa2);
                    CEcoa4 = -2 * p_coa3 *
                      (workspace->bo_dboc[k][0]-BOA_jk) * e_coa;
                    CEcoa5 = -2 * p_coa3 *
                      (workspace->bo_dboc[h][0]-BOA_hj) * e_coa;
                    /* END COALITION ENERGY */

                    /* FORCES */
                    //bo_jk->Cdbo += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));
                    //bo_hj->Cdbo += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));
                    Cdbo_list[pk] += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));
                    Cdbo_list[ph] += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));

                    workspace->fCdDelta[j][3] += ((CEval3 + CEval7) + CEpen1 + CEcoa3);
                    workspace->fCdDelta[k][3] += CEcoa4;
                    workspace->fCdDelta[h][3] += CEcoa5;

                    for( t = start_j; t < end_j; ++t ) 
                    {
                      pbond_jt = &( bonds->select.bond_list[t] );
                      //bo_jt = &(pbond_jt->bo_data);
                      bo_jt = &(bonds->bo_data_list[t]);
                      //temp_bo_jt = bo_jt->BO;
                      temp_bo_jt = BO_list[t];
                      temp = CUBE( temp_bo_jt );
                      pBOjt7 = temp * temp * temp_bo_jt;

                      //bo_jt->Cdbo += (CEval6 * pBOjt7);
                      Cdbo_list[t] += (CEval6 * pBOjt7);
                      //bo_jt->Cdbopi += CEval5;
                      Cdbopi_list[t] += CEval5;
                      //bo_jt->Cdbopi2 += CEval5;
                      Cdbopi2_list[t] += CEval5;
                    }

                   rvec4_ScaledAdd( workspace->fCdDelta[k], CEval8, p_hjk->dcos_di );
                   rvec4_ScaledAdd( workspace->fCdDelta[j], CEval8, p_hjk->dcos_dj );
                   rvec4_ScaledAdd( workspace->fCdDelta[h], CEval8, p_hjk->dcos_dk );
                    
                    /* tally into per-atom virials */
                    if(system->evflag)
                    {
                      eng_tmp = e_ang + e_pen + e_coa;
                      system->eng_vdwl += eng_tmp;
                    }
                  }
                }
              }

              //Torsion
              if (system->packed_atoms[j].orig_id >  system->packed_atoms[k].orig_id)
                continue;
              if (system->packed_atoms[j].orig_id == system->packed_atoms[k].orig_id) 
              {
                if (system->packed_atoms[k].x[2] <  system->packed_atoms[j].x[2]) continue;
                if (system->packed_atoms[k].x[2] == system->packed_atoms[j].x[2] &&
                    system->packed_atoms[k].x[1] <  system->packed_atoms[j].x[1]) continue;
                if (system->packed_atoms[k].x[2] == system->packed_atoms[j].x[2] &&
                    system->packed_atoms[k].x[1] == system->packed_atoms[j].x[1] &&
                    system->packed_atoms[k].x[0] <  system->packed_atoms[j].x[0]) continue;
              }


              //if (bo_ij->BO > control->thb_cut /*0*/ && j < system->n && end_k > start_k) 
              if (BO_list[ph] > control->thb_cut /*0*/ && j < system->n && end_k > start_k) 
              {
                i = p_ijk->thb;
                type_i = system->packed_atoms[i].type;
                r_ij = pbond_ij->d;
                //BOA_ij = bo_ij->BO - control->thb_cut;
                BOA_ij = BO_list[ph] - control->thb_cut;

                theta_ijk = p_ijk->theta;
                sin_ijk = sin( theta_ijk );
                cos_ijk = cos( theta_ijk );
                //tan_ijk_i = 1. / tan( theta_ijk );
                if (sin_ijk >= 0 && sin_ijk <= MIN_SINE)
                  tan_ijk_i = cos_ijk / MIN_SINE;
                else if( sin_ijk <= 0 && sin_ijk >= -MIN_SINE )
                  tan_ijk_i = cos_ijk / -MIN_SINE;
                else tan_ijk_i = cos_ijk / sin_ijk;

                exp_tor2_ij = exp( -p_tor2 * BOA_ij );
                exp_cot2_ij = exp( -p_cot2 * SQR(BOA_ij -1.5) );
                
                for(pw = start_k; pw < end_k; pw++)
                {
                  pbond_kw = &(bonds->select.bond_list[pw]);
                  //bo_kw = &(pbond_kw->bo_data);
                  bo_kw = &(bonds->bo_data_list[pw]);
                  //BOA_kw = bo_kw->BO - control->thb_cut;
                  BOA_kw = BO_list[pw] - control->thb_cut;
                  w = pbond_kw->nbr;
                  if(pw != pj)
                  {
                    type_w = system->packed_atoms[w].type;
                    p_jkw   = &thb_intrs_data;
                    
                    Calculate_Theta( pbond_pj->dvec, pbond_pj->d,
                                     pbond_kw->dvec, pbond_kw->d,
                                     &theta_jkw, &cos_theta_jkw);

                    Calculate_dCos_Theta( pbond_pj->dvec, pbond_pj->d,
                                          pbond_kw->dvec, pbond_kw->d,
                                          &(p_jkw->dcos_di), &(p_jkw->dcos_dj),
                                          &(p_jkw->dcos_dk));
                        
                    p_jkw->thb = w;
                    p_jkw->pthb = pw;
                    p_jkw->theta = theta_jkw;
                    

                    l = w;
                    type_l = type_w;
                    theta_jkl = p_jkw->theta;
                    pbond_kl = &(bonds->select.bond_list[pw]);
                    //bo_kl = &(pbond_kl->bo_data);
                    bo_kl = &(bonds->bo_data_list[pw]);
                    p_jkl = p_jkw;
                    BOA_kl = BOA_kw;
                    
                    fbh = &(system->reax_param.fbp[type_i][type_j]
                            [type_k][type_l]);
                    fbp = &(system->reax_param.fbp[type_i][type_j]
                            [type_k][type_l].prm[0]);

              
                    //printf("torsion:(j, k, i, l) = (%d, %d, %d, %d)\n", j, k, i, l);

                    //if( i != l && fbh->cnt &&
                    //    bo_kl->BO > control->thb_cut/*0*/ &&
                    //    bo_ij->BO * bo_jk->BO * bo_kl->BO > control->thb_cut/*0*/ )
                    if( i != l && fbh->cnt &&
                          BO_list[pw] > control->thb_cut/*0*/ &&
                          BO_list[ph] * BO_list[pk] * BO_list[pw] > control->thb_cut/*0*/ )

                    {
                      //++num_frb_intrs;
                      r_kl = pbond_kl->d;
                      //BOA_kl = bo_kl->BO - control->thb_cut;
                      BOA_kl = BO_list[pw] - control->thb_cut;

                      //theta_jkl = p_jkl->theta;
                      sin_jkl = sin( theta_jkl );
                      cos_jkl = cos( theta_jkl );
                      //tan_jkl_i = 1. / tan( theta_jkl );
                      if (sin_jkl >= 0 && sin_jkl <= MIN_SINE)
                        tan_jkl_i = cos_jkl / MIN_SINE;
                      else if( sin_jkl <= 0 && sin_jkl >= -MIN_SINE )
                        tan_jkl_i = cos_jkl / -MIN_SINE;
                      else tan_jkl_i = cos_jkl /sin_jkl;

                      rvec_ScaledSum( dvec_li, 1., system->packed_atoms[i].x,
                                      -1., system->packed_atoms[l].x );
                      r_li = sqrt(SQR(dvec_li[0]) + SQR(dvec_li[1]) + SQR(dvec_li[2]));


                      /* omega and its derivative */
                      omega = Calculate_Omega( pbond_ij->dvec, r_ij,
                                               pbond_jk->dvec, r_jk,
                                               pbond_kl->dvec, r_kl,
                                               dvec_li, r_li,
                                               p_ijk, p_jkl,
                                               dcos_omega_di, dcos_omega_dj,
                                               dcos_omega_dk, dcos_omega_dl,
                                               out_control );

                      cos_omega = cos( omega );
                      cos2omega = cos( 2. * omega );
                      cos3omega = cos( 3. * omega );
                      /* end omega calculations */

                      /* torsion energy */
                      exp_tor1 = exp( fbp->p_tor1 *
                                      SQR(2.0 - BOpi_list[pk][0] - f11_DjDk));
                                      //SQR(2.0 - bo_jk->BO_pi - f11_DjDk) );
                      exp_tor2_kl = exp( -p_tor2 * BOA_kl );
                      exp_cot2_kl = exp( -p_cot2 * SQR(BOA_kl - 1.5) );
                      fn10 = (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk) *
                        (1.0 - exp_tor2_kl);

                      CV = 0.5 * ( fbp->V1 * (1.0 + cos_omega) +
                                   fbp->V2 * exp_tor1 * (1.0 - cos2omega) +
                                   fbp->V3 * (1.0 + cos3omega) );

                      data->my_en.e_tor += e_tor = fn10 * sin_ijk * sin_jkl * CV;

                      dfn11 = (-p_tor3 * exp_tor3_DjDk +
                               (p_tor3 * exp_tor3_DjDk - p_tor4 * exp_tor4_DjDk) *
                               (2.0 + exp_tor3_DjDk) * exp_tor34_inv) *
                        exp_tor34_inv;

                      CEtors1 = sin_ijk * sin_jkl * CV;

                      //CEtors2 = -fn10 * 2.0 * fbp->p_tor1 * fbp->V2 * exp_tor1 *
                      //  (2.0 - bo_jk->BO_pi - f11_DjDk) * (1.0 - SQR(cos_omega)) *
                      //  sin_ijk * sin_jkl;
                      CEtors2 = -fn10 * 2.0 * fbp->p_tor1 * fbp->V2 * exp_tor1 *
                          (2.0 - BOpi_list[pk][0] -f11_DjDk) * (1.0 -SQR(cos_omega)) *
                          sin_ijk * sin_jkl;

                      CEtors3 = CEtors2 * dfn11;

                      CEtors4 = CEtors1 * p_tor2 * exp_tor2_ij *
                        (1.0 - exp_tor2_jk) * (1.0 - exp_tor2_kl);
                      CEtors5 = CEtors1 * p_tor2 *
                        (1.0 - exp_tor2_ij) * exp_tor2_jk * (1.0 - exp_tor2_kl);
                      CEtors6 = CEtors1 * p_tor2 *
                        (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk) * exp_tor2_kl;

                      cmn = -fn10 * CV;
                      CEtors7 = cmn * sin_jkl * tan_ijk_i;
                      CEtors8 = cmn * sin_ijk * tan_jkl_i;

                      CEtors9 = fn10 * sin_ijk * sin_jkl *
                        (0.5 * fbp->V1 - 2.0 * fbp->V2 * exp_tor1 * cos_omega +
                         1.5 * fbp->V3 * (cos2omega + 2.0 * SQR(cos_omega)));
                      /* end  of torsion energy */

                      /* 4-body conjugation energy */
                      fn12 = exp_cot2_ij * exp_cot2_jk * exp_cot2_kl;
                      data->my_en.e_con += e_con =
                        fbp->p_cot1 * fn12 *
                        (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);

                      Cconj = -2.0 * fn12 * fbp->p_cot1 * p_cot2 *
                        (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);

                      CEconj1 = Cconj * (BOA_ij - 1.5e0);
                      CEconj2 = Cconj * (BOA_jk - 1.5e0);
                      CEconj3 = Cconj * (BOA_kl - 1.5e0);

                      CEconj4 = -fbp->p_cot1 * fn12 *
                        (SQR(cos_omega) - 1.0) * sin_jkl * tan_ijk_i;
                      CEconj5 = -fbp->p_cot1 * fn12 *
                        (SQR(cos_omega) - 1.0) * sin_ijk * tan_jkl_i;
                      CEconj6 = 2.0 * fbp->p_cot1 * fn12 *
                        cos_omega * sin_ijk * sin_jkl;
                      /* end 4-body conjugation energy */

                      /* forces */
                      //bo_jk->Cdbopi += CEtors2;
                      Cdbopi_list[pk] += CEtors2;
                      workspace->fCdDelta[j][3] += CEtors3;
                      workspace->fCdDelta[k][3] += CEtors3;
                      //bo_ij->Cdbo += (CEtors4 + CEconj1);
                      //bo_jk->Cdbo += (CEtors5 + CEconj2);
                      //bo_kl->Cdbo += (CEtors6 + CEconj3);
                      Cdbo_list[ph] += (CEtors4 + CEconj1);
                      Cdbo_list[pk] += (CEtors5 + CEconj2);
                      Cdbo_list[pw] += (CEtors6 + CEconj3);


                      /* dcos_theta_ijk */
                      rvec4_ScaledAdd( workspace->fCdDelta[i],
                                      CEtors7 + CEconj4, p_ijk->dcos_dk );
                      rvec4_ScaledAdd( workspace->fCdDelta[j],
                                      CEtors7 + CEconj4, p_ijk->dcos_dj );
                      rvec4_ScaledAdd( workspace->fCdDelta[k],
                                      CEtors7 + CEconj4, p_ijk->dcos_di );

                      /* dcos_theta_jkl */
                      rvec4_ScaledAdd( workspace->fCdDelta[j],
                                      CEtors8 + CEconj5, p_jkl->dcos_di );
                      rvec4_ScaledAdd( workspace->fCdDelta[k],
                                      CEtors8 + CEconj5, p_jkl->dcos_dj );
                      rvec4_ScaledAdd( workspace->fCdDelta[l],
                                      CEtors8 + CEconj5, p_jkl->dcos_dk );

                      /* dcos_omega */
                      rvec4_ScaledAdd( workspace->fCdDelta[i],
                                      CEtors9 + CEconj6, dcos_omega_di );
                      rvec4_ScaledAdd( workspace->fCdDelta[j],
                                      CEtors9 + CEconj6, dcos_omega_dj );
                      rvec4_ScaledAdd( workspace->fCdDelta[k],
                                      CEtors9 + CEconj6, dcos_omega_dk );
                      rvec4_ScaledAdd( workspace->fCdDelta[l],
                                      CEtors9 + CEconj6, dcos_omega_dl );
                      
                      /* tally into per-atom virials */
                      if(system->evflag)
                      {
                        eng_tmp = e_tor + e_con;
                        system->eng_vdwl += eng_tmp;
                      }

                    } // pl check ends
                  }//if pw!=pj
                }//for-pw
              } // pi check ends
            } // if ph!=pi
          }//for-ph
        } // k-j neighbor check ends
      } // j-k neighbor check ends
    } // pk loop ends
  } // j loop
#endif
}
#endif



#ifdef CPE
#include "STUBS/mpi.h"
#include <dma.h>
#include "slave.h"

#define DMA_FAST
#include "poly_math.h"
#include "dma_macros.h"

//#define LWPF_UNIT U(TORSION)
//#define LWPF_KERNELS K(ALL)//  K(CMP) K(CMP_TOR)  
//#define DMA_FAST
//#include "lwpf2.h"

//read cache;
#define READ_C_H    7
#define READ_C_S    4
#define READ_C_LSZ  (1 << READ_C_S)
#define READ_C_LCNT (1 << (READ_C_H - READ_C_S))
#define READ_C_MM   (READ_C_LSZ - 1)
#define READ_C_LM   (READ_C_LCNT - 1)

//bo_dboc;
#define DT_C_H    5
#define DT_C_S    3
#define DT_C_LSZ  (1 << DT_C_S)
#define DT_C_LCNT (1 << (DT_C_H - DT_C_S))
#define DT_C_MM   (DT_C_LSZ - 1)
#define DT_C_LM   (DT_C_LCNT - 1)

//index, end_index for atom k;
#define K_C_H    5
#define K_C_S    3
#define K_C_LSZ  (1 << K_C_S)
#define K_C_LCNT (1 << (K_C_H - K_C_S))
#define K_C_MM   (K_C_LSZ - 1)
#define K_C_LM   (K_C_LCNT - 1)

//read bond_data for pj;
#define BOND_C_H    5
#define BOND_C_S    3
#define BOND_C_LSZ  (1 << BOND_C_S)
#define BOND_C_LCNT (1 << (BOND_C_H - BOND_C_S))
#define BOND_C_MM   (BOND_C_LSZ - 1)
#define BOND_C_LM   (BOND_C_LCNT - 1)

void read_cache(int i,
                atom_pack_t atoms_cache[][READ_C_LSZ], 
                int *read_ctag,
                atom_pack_t *packed_atoms,
                atom_pack_t *atom_i)
{
  dma_init();
  if (read_ctag[(i >> READ_C_S) & READ_C_LM] != i >> READ_C_S)
  {
    pe_get(packed_atoms + (i & ~READ_C_MM), 
          atoms_cache[(i >> READ_C_S) & READ_C_LM], 
          sizeof(atom_pack_t) * READ_C_LSZ);
    dma_syn();
    read_ctag[(i >> READ_C_S) & READ_C_LM] = i >> READ_C_S;
  }
  *atom_i = atoms_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
}

void read_cache_bo_dboc(int k, 
                        rvec2 bo_dboc_cache [][DT_C_LSZ], 
                        int *bo_dboc_tag,
                        rvec2 *bo_dboc,
                        rvec2 bo_dboc_k)
{
  dma_init();
  if (bo_dboc_tag[(k >> DT_C_S) & DT_C_LM] != k >> DT_C_S)
  {
    pe_get(bo_dboc + (k & ~DT_C_MM), 
          bo_dboc_cache[(k >> DT_C_S) & DT_C_LM], 
          sizeof(rvec2) * DT_C_LSZ);
    dma_syn();
    bo_dboc_tag[(k >> DT_C_S) & DT_C_LM] = k >> DT_C_S;
  }
  (bo_dboc_k)[0] = bo_dboc_cache[(k >> DT_C_S) & DT_C_LM][k & DT_C_MM][0];
  (bo_dboc_k)[1] = bo_dboc_cache[(k >> DT_C_S) & DT_C_LM][k & DT_C_MM][1];
}

void read_cache_k_index(int k, 
                        int k_idx_cache [][K_C_LSZ], 
                        int k_ed_idx_cache [][K_C_LSZ], 
                        int *k_tag,
                        int *k_index,
                        int *k_end_index,
                        int *k_st, int *k_ed)
{
  dma_init();
  if (k_tag[(k >> K_C_S) & K_C_LM] != k >> K_C_S)
  {
    pe_get(k_index + (k & ~K_C_MM), 
          k_idx_cache[(k >> K_C_S) & K_C_LM], 
          sizeof(int) * K_C_LSZ);
    pe_get(k_end_index + (k & ~K_C_MM), 
          k_ed_idx_cache[(k >> K_C_S) & K_C_LM], 
          sizeof(int) * K_C_LSZ);

    dma_syn();
    k_tag[(k >> K_C_S) & K_C_LM] = k >> K_C_S;
  }
  *k_st = k_idx_cache[   (k >> K_C_S) & K_C_LM][k & K_C_MM];
  *k_ed = k_ed_idx_cache[(k >> K_C_S) & K_C_LM][k & K_C_MM];
}

void read_cache_bond(int pj,
                bond_data bond_cache[][BOND_C_LSZ], 
                int *bond_tag,
                bond_data *bond_list,
                bond_data *bond_pj)
{
  dma_init();
  if (bond_tag[(pj >> BOND_C_S) & BOND_C_LM] != pj >> BOND_C_S)
  {
    pe_get(bond_list + (pj & ~BOND_C_MM), 
          bond_cache[(pj >> BOND_C_S) & BOND_C_LM], 
          sizeof(bond_data) * BOND_C_LSZ);
    dma_syn();
    bond_tag[(pj >> BOND_C_S) & BOND_C_LM] = pj >> BOND_C_S;
  }
  *bond_pj = bond_cache[(pj >> BOND_C_S) & BOND_C_LM][pj & BOND_C_MM];
}


double Dot( double* v1, double* v2, int k )
{
  double ret = 0.0;
  int i;
  for( i=0; i < k; ++i )
    ret +=  v1[i] * v2[i];

  return ret;
}

void Calculate_Theta( rvec dvec_ji, double d_ji, rvec dvec_jk, double d_jk,
                      double *theta, double *cos_theta )
{
  (*cos_theta) = Dot( dvec_ji, dvec_jk, 3 ) / ( d_ji * d_jk );
  if( *cos_theta > 1. ) *cos_theta  = 1.0;
  if( *cos_theta < -1. ) *cos_theta  = -1.0;

  //(*theta) = acos( *cos_theta );
  (*theta) = xacos( *cos_theta );
  //(*theta) = acos_sw( *cos_theta );
  //if(myrank == 0)
  //  printf("acos(%.10f) = %.10f, %.10f\n", *cos_theta, *theta, acos(*cos_theta));
  
}

void Calculate_dCos_Theta( rvec dvec_ji, double d_ji, rvec dvec_jk, double d_jk,
                           rvec* dcos_theta_di,
                           rvec* dcos_theta_dj,
                           rvec* dcos_theta_dk )
{
  int t;
  double sqr_d_ji = SQR(d_ji);
  double sqr_d_jk = SQR(d_jk);
  double inv_dists = 1.0 / (d_ji * d_jk);
  double inv_dists3 = p_powd( inv_dists, 3.0 );
  double dot_dvecs = Dot( dvec_ji, dvec_jk, 3 );
  double Cdot_inv3 = dot_dvecs * inv_dists3;

  for( t = 0; t < 3; ++t ) 
  {
    (*dcos_theta_di)[t] = dvec_jk[t] * inv_dists -
      Cdot_inv3 * sqr_d_jk * dvec_ji[t];
    (*dcos_theta_dj)[t] = -(dvec_jk[t] + dvec_ji[t]) * inv_dists +
      Cdot_inv3 * ( sqr_d_jk * dvec_ji[t] + sqr_d_ji * dvec_jk[t] );
    (*dcos_theta_dk)[t] = dvec_ji[t] * inv_dists -
      Cdot_inv3 * sqr_d_ji * dvec_jk[t];
  }
}

double Calculate_Omega_cpe( rvec dvec_ij, double r_ij,
                      rvec dvec_jk, double r_jk,
                      rvec dvec_kl, double r_kl,
                      rvec dvec_li, double r_li,
                      three_body_interaction_data *p_ijk,
                      three_body_interaction_data *p_jkl,
                      rvec dcos_omega_di, rvec dcos_omega_dj,
                      rvec dcos_omega_dk, rvec dcos_omega_dl)
{
  double unnorm_cos_omega, unnorm_sin_omega, omega;
  double sin_ijk, cos_ijk, sin_jkl, cos_jkl;
  double htra, htrb, htrc, hthd, hthe, hnra, hnrc, hnhd, hnhe;
  double arg, poem, tel;
  rvec cross_jk_kl;

  sin_ijk = p_sind(p_ijk->theta);
  cos_ijk = p_cosd(p_ijk->theta);
  sin_jkl = p_sind(p_jkl->theta);
  cos_jkl = p_cosd(p_jkl->theta);

  /* omega */
  unnorm_cos_omega = -rvec_Dot(dvec_ij, dvec_jk) * rvec_Dot(dvec_jk, dvec_kl) +
                    SQR(r_jk) *  rvec_Dot(dvec_ij, dvec_kl);

  rvec_Cross(cross_jk_kl, dvec_jk, dvec_kl);
  unnorm_sin_omega = -r_jk * rvec_Dot( dvec_ij, cross_jk_kl );

  //omega = atan2( unnorm_sin_omega, unnorm_cos_omega );
  omega = xatan2( unnorm_sin_omega, unnorm_cos_omega );

  htra = r_ij + cos_ijk * (r_kl * cos_jkl - r_jk);
  htrb = r_jk - r_ij * cos_ijk - r_kl * cos_jkl;
  htrc = r_kl + cos_jkl * (r_ij * cos_ijk - r_jk);
  hthd = r_ij * sin_ijk * (r_jk - r_kl * cos_jkl);
  hthe = r_kl * sin_jkl * (r_jk - r_ij * cos_ijk);
  hnra = r_kl * sin_ijk * sin_jkl;
  hnrc = r_ij * sin_ijk * sin_jkl;
  hnhd = r_ij * r_kl * cos_ijk * sin_jkl;
  hnhe = r_ij * r_kl * sin_ijk * cos_jkl;

  poem = 2.0 * r_ij * r_kl * sin_ijk * sin_jkl;
  if( poem < 1e-20 ) poem = 1e-20;

  tel  = SQR( r_ij ) + SQR( r_jk ) + SQR( r_kl ) - SQR( r_li ) -
    2.0 * ( r_ij * r_jk * cos_ijk - r_ij * r_kl * cos_ijk * cos_jkl +
            r_jk * r_kl * cos_jkl );

  arg  = tel / poem;
  if( arg >  1.0 ) arg =  1.0;
  if( arg < -1.0 ) arg = -1.0;

  if( sin_ijk >= 0 && sin_ijk <= MIN_SINE ) sin_ijk = MIN_SINE;
  else if( sin_ijk <= 0 && sin_ijk >= -MIN_SINE ) sin_ijk = -MIN_SINE;
  if( sin_jkl >= 0 && sin_jkl <= MIN_SINE ) sin_jkl = MIN_SINE;
  else if( sin_jkl <= 0 && sin_jkl >= -MIN_SINE ) sin_jkl = -MIN_SINE;

  // dcos_omega_di
  rvec_ScaledSum( dcos_omega_di, (htra-arg*hnra)/r_ij, dvec_ij, -1., dvec_li );
  rvec_ScaledAdd( dcos_omega_di,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_dk );
  rvec_Scale( dcos_omega_di, 2.0 / poem, dcos_omega_di );

  // dcos_omega_dj
  rvec_ScaledSum( dcos_omega_dj,-(htra-arg*hnra)/r_ij, dvec_ij,
                  -htrb / r_jk, dvec_jk );
  rvec_ScaledAdd( dcos_omega_dj,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_dj );
  rvec_ScaledAdd( dcos_omega_dj,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_di );
  rvec_Scale( dcos_omega_dj, 2.0 / poem, dcos_omega_dj );

  // dcos_omega_dk
  rvec_ScaledSum( dcos_omega_dk,-(htrc-arg*hnrc)/r_kl, dvec_kl,
                  htrb / r_jk, dvec_jk );
  rvec_ScaledAdd( dcos_omega_dk,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_di );
  rvec_ScaledAdd( dcos_omega_dk,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_dj );
  rvec_Scale( dcos_omega_dk, 2.0 / poem, dcos_omega_dk );

  // dcos_omega_dl
  rvec_ScaledSum( dcos_omega_dl, (htrc-arg*hnrc)/r_kl, dvec_kl, 1., dvec_li );
  rvec_ScaledAdd( dcos_omega_dl,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_dk );
  rvec_Scale( dcos_omega_dl, 2.0 / poem, dcos_omega_dl );

  return omega;
}



#define DT_CACHE Cdbo_list,double,3,3,6,FADDD
#define DT_CACHE_F fCdDelta,rvec4,5,3,6,FADDD
#define EVFLAG 1
#include "reaxc_torsion_angles_cpe.h"
#undef EVFLAG
#define EVFLAG 0
#include "reaxc_torsion_angles_cpe.h"
#undef EVFLAG
#endif
