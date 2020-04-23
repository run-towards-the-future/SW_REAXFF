/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_reaxc_sunway.h"
#include "reaxc_hydrogen_bonds_sunway.h"
#include "reaxc_bond_orders_sunway.h"
#include "reaxc_list_sunway.h"
#include "reaxc_valence_angles_sunway.h"
#include "reaxc_vector_sunway.h"
#include "gptl.h"
#include "reaxc_ctypes_sunway.h"
#include "reaxc_hydrogen_bonds_sw64.h"


namespace REAXC_SUNWAY_NS{
extern"C"
{
  void Hydrogen_Bonds_C(void *param);
}


//void Hydrogen_Bonds_org( reax_system *system, control_params *control,
//                     simulation_data *data, storage *workspace,
//                     reax_list **lists, output_controls *out_control )
//{
////  reax_system_c csys;
////  system->to_c_sys(&csys);
////
////  //Hydrogen_Bonds_C(&csys, control, data, workspace, lists, out_control);
////  
////  param_pack_t param;
////  param.system      = &csys;
////  param.control     = control;
////  param.data        = data;
////  param.workspace   = workspace;
////  param.lists       = lists;
////  param.out_control = out_control;
////  GPTLstart("reaxc hydrogen bonds c");
////  Hydrogen_Bonds_C(&param);
////  GPTLstop("reaxc hydrogen bonds c");
////
////  
////  system->from_c_sys(&csys);
////  return;
//
//  int  i, j, k, pi, pk;
//  int  type_i, type_j, type_k;
//  int  start_j, end_j, hb_start_j, hb_end_j;
//  int  hblist[MAX_BONDS];
//  int  itr, top;
//  int  num_hb_intrs = 0;
//  ivec rel_jk;
//  double r_jk, theta, cos_theta, sin_xhz4, cos_xhz1, sin_theta2;
//  double e_hb, exp_hb2, exp_hb3, CEhb1, CEhb2, CEhb3;
//  rvec dcos_theta_di, dcos_theta_dj, dcos_theta_dk;
//  rvec dvec_jk, force, ext_press;
//  hbond_parameters *hbp;
//  bond_order_data *bo_ij;
//  bond_data *pbond_ij;
//  far_neighbor_data_full *nbr_jk;
//  reax_list *bonds, *hbonds;
//  bond_data *bond_list;
//  hbond_data *hbond_list;
//
//  // tally variables
//  double fi_tmp[3], fk_tmp[3], delij[3], delkj[3];
//
//  bonds = (*lists) + BONDS;
//  bond_list = bonds->select.bond_list;
//  double *Cdbo_list = bonds->Cdbo_list;
//  double *BO_list   = bonds->BO_list;
//
//  hbonds = (*lists) + HBONDS;
//  hbond_list = hbonds->select.hbond_list;
//
//  int myrank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//
//  for( j = 0; j < system->n; ++j )
//  {
//    if( system->reax_param.sbp[system->packed_atoms[j].type].p_hbond == 1 ) 
//    {
//      //if(myrank == 0)
//      //  printf("j = %d, n= %d, N= %d\n", j, system->n, system->N);
//      type_j     = system->packed_atoms[j].type;
//      if (type_j < 0) continue;
//      start_j    = Start_Index(j, bonds);
//      end_j      = End_Index(j, bonds);
//      //hb_start_j = Start_Index( system->my_atoms[j].Hindex, hbonds );
//      //hb_end_j   = End_Index( system->my_atoms[j].Hindex, hbonds );
//      
//      hb_start_j = Start_Index( system->Hindex[j], hbonds );
//      hb_end_j   = End_Index  ( system->Hindex[j], hbonds );
//
//
//      top = 0;
//      for( pi = start_j; pi < end_j; ++pi )  
//      {
//        pbond_ij = &( bond_list[pi] );
//        i = pbond_ij->nbr;
//        type_i = system->packed_atoms[i].type;
//	    if (type_i < 0) continue;
//        //bo_ij = &(pbond_ij->bo_data);
//        bo_ij = &(bonds->bo_data_list[pi]);
//
//        //if( system->reax_param.sbp[type_i].p_hbond == 2 && bo_ij->BO >= HB_THRESHOLD )
//        if(system->reax_param.sbp[type_i].p_hbond == 2 && BO_list[pi] >= HB_THRESHOLD)
//          hblist[top++] = pi;
//      }
//
//      for( pk = hb_start_j; pk < hb_end_j; ++pk ) 
//      {
//        /* set k's varibles */
//        k = hbond_list[pk].nbr;
//        type_k = system->packed_atoms[k].type;
//	      if (type_k < 0) continue;
//        nbr_jk = hbond_list[pk].ptr;
//        r_jk = nbr_jk->d;
//        rvec_Scale( dvec_jk, hbond_list[pk].scl, nbr_jk->dvec );
//
//        for( itr = 0; itr < top; ++itr ) 
//        {
//          pi = hblist[itr];
//          pbond_ij = &( bonds->select.bond_list[pi] );
//          i = pbond_ij->nbr;
//
//          if( system->packed_atoms[i].orig_id != system->packed_atoms[k].orig_id ) 
//          {
//            //bo_ij = &(pbond_ij->bo_data);
//            bo_ij = &(bonds->bo_data_list[pi]);
//            type_i = system->packed_atoms[i].type;
//	        if (type_i < 0) continue;
//                hbp = &(system->reax_param.hbp[ type_i ][ type_j ][ type_k ]);
//	        if (hbp->r0_hb <= 0.0) continue;
//                ++num_hb_intrs;
//
//            Calculate_Theta( pbond_ij->dvec, pbond_ij->d, dvec_jk, r_jk,
//                             &theta, &cos_theta );
//            /* the derivative of cos(theta) */
//            Calculate_dCos_Theta( pbond_ij->dvec, pbond_ij->d, dvec_jk, r_jk,
//                                  &dcos_theta_di, &dcos_theta_dj,
//                                  &dcos_theta_dk );
//
//            /* hyrogen bond energy*/
//            sin_theta2 = sin( theta/2.0 );
//            sin_xhz4 = SQR(sin_theta2);
//            sin_xhz4 *= sin_xhz4;
//            cos_xhz1 = ( 1.0 - cos_theta );
//            //exp_hb2 = exp( -hbp->p_hb2 * bo_ij->BO );
//            exp_hb2 = exp( -hbp->p_hb2 * BO_list[pi]);
//            exp_hb3 = exp( -hbp->p_hb3 * ( hbp->r0_hb / r_jk +
//                                           r_jk / hbp->r0_hb - 2.0 ) );
//
//            data->my_en.e_hb += e_hb =
//              hbp->p_hb1 * (1.0 - exp_hb2) * exp_hb3 * sin_xhz4;
//
//            CEhb1 = hbp->p_hb1 * hbp->p_hb2 * exp_hb2 * exp_hb3 * sin_xhz4;
//            CEhb2 = -hbp->p_hb1/2.0 * (1.0 - exp_hb2) * exp_hb3 * cos_xhz1;
//            CEhb3 = -hbp->p_hb3 *
//              (-hbp->r0_hb / SQR(r_jk) + 1.0 / hbp->r0_hb) * e_hb;
//
//            /* hydrogen bond forces */
//            //bo_ij->Cdbo += CEhb1; // dbo term
//            Cdbo_list[pi] += CEhb1; // dbo term
//
//            if( control->virial == 0 ) 
//            {
//              // dcos terms
//              rvec4_ScaledAdd( workspace->fCdDelta[i], +CEhb2, dcos_theta_di );
//              rvec4_ScaledAdd( workspace->fCdDelta[j], +CEhb2, dcos_theta_dj );
//              rvec4_ScaledAdd( workspace->fCdDelta[k], +CEhb2, dcos_theta_dk );
//              // dr terms
//              rvec4_ScaledAdd( workspace->fCdDelta[j], -CEhb3/r_jk, dvec_jk );
//              rvec4_ScaledAdd( workspace->fCdDelta[k], +CEhb3/r_jk, dvec_jk );
//            }
//            else 
//            {
//              rvec_Scale( force, +CEhb2, dcos_theta_di ); // dcos terms
//              rvec4_Add( workspace->fCdDelta[i], force );
//              rvec_iMultiply( ext_press, pbond_ij->rel_box, force );
//              rvec_ScaledAdd( data->my_ext_press, 1.0, ext_press );
//
//              rvec4_ScaledAdd( workspace->fCdDelta[j], +CEhb2, dcos_theta_dj );
//
//              ivec_Scale( rel_jk, hbond_list[pk].scl, nbr_jk->rel_box );
//              rvec_Scale( force, +CEhb2, dcos_theta_dk );
//              rvec4_Add( workspace->fCdDelta[k], force );
//              rvec_iMultiply( ext_press, rel_jk, force );
//              rvec_ScaledAdd( data->my_ext_press, 1.0, ext_press );
//              // dr terms
//              rvec4_ScaledAdd( workspace->fCdDelta[j], -CEhb3/r_jk, dvec_jk );
//
//              rvec_Scale( force, CEhb3/r_jk, dvec_jk );
//              rvec4_Add( workspace->fCdDelta[k], force );
//              rvec_iMultiply( ext_press, rel_jk, force );
//              rvec_ScaledAdd( data->my_ext_press, 1.0, ext_press );
//            }
//
//            /* tally into per-atom virials */
//            if (system->pair_ptr->vflag_atom || system->pair_ptr->evflag) 
//            {
//              rvec_ScaledSum( delij, 1., system->packed_atoms[j].x,
//                                    -1., system->packed_atoms[i].x );
//              rvec_ScaledSum( delkj, 1., system->packed_atoms[j].x,
//                                     -1.,system->packed_atoms[k].x );
//
//              rvec_Scale(fi_tmp, CEhb2, dcos_theta_di);
//              rvec_Scale(fk_tmp, CEhb2, dcos_theta_dk);
//              rvec_ScaledAdd(fk_tmp, CEhb3/r_jk, dvec_jk);
//
//              system->pair_ptr->ev_tally3(i,j,k,e_hb,0.0,fi_tmp,fk_tmp,delij,delkj);
//            }
//          }
//        }
//      }
//    }//if
//  }//for
//}

void Hydrogen_Bonds( reax_system *system, control_params *control,
                     simulation_data *data, storage *workspace,
                     reax_list **lists, output_controls *out_control )
{
  int eflag_either  = system->pair_ptr->eflag_either;
  int eflag_atom    = system->pair_ptr->eflag_atom;
  int eflag_global  = system->pair_ptr->eflag_global;
  int vflag_either  = system->pair_ptr->vflag_either;
  int vflag_atom    = system->pair_ptr->vflag_atom;
  int vflag_global  = system->pair_ptr->vflag_global;

  if(control->virial == 0)
  {
    reax_system_c csys;
    system->to_c_sys(&csys);
    
    hb_param_pack_t param;
    param.system      = &csys;
    param.control     = control;
    param.data        = data;
    param.workspace   = workspace;
    param.lists       = lists;
    param.out_control = out_control;

    Hydrogen_Bonds_C(&param);
          
    system->from_c_sys(&csys);

    return;
  }
  else
  {/*****NPT******/
    printf("Hydrogen Bonds on MPE\n");
    printf("evflag =%d, eflag_either=%d, eflag_atom=%d, eflag_global=%d, vflag_either=%d,vflag_atom=%d, vflag_global=%d\n", system->pair_ptr->evflag, system->pair_ptr->eflag_either, system->pair_ptr->eflag_atom, system->pair_ptr->eflag_global, system->pair_ptr->vflag_either, system->pair_ptr->vflag_atom, system->pair_ptr->vflag_global);

    int  i, j, k, pi, pk;
    int  type_i, type_j, type_k;
    int  start_j, end_j, hb_start_j, hb_end_j;
    int  hblist[MAX_BONDS];
    int  itr, top;
    ivec rel_jk;
    double r_jk, theta, cos_theta, sin_xhz4, cos_xhz1, sin_theta2;
    double e_hb, exp_hb2, exp_hb3, CEhb1, CEhb2, CEhb3;
    rvec dcos_theta_di, dcos_theta_dj, dcos_theta_dk;
    rvec dvec_jk, force, ext_press;
    hbond_parameters *hbp;
    bond_order_data *bo_ij;
    bond_data *pbond_ij;
    far_neighbor_data_full *nbr_jk;
    reax_list *bonds, *hbonds;
    bond_data *bond_list;
    hbond_data *hbond_list;

    // tally variables
    double fi_tmp[3], fk_tmp[3], delij[3], delkj[3];

    bonds = (*lists) + BONDS;
    bond_list = bonds->select.bond_list;
    double *Cdbo_list = bonds->Cdbo_list;
    double *BO_list   = bonds->BO_list;

    hbonds = (*lists) + HBONDS;
    hbond_list = hbonds->select.hbond_list;

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    for( j = 0; j < system->n; ++j )
    {
      if( system->reax_param.sbp[system->packed_atoms[j].type].p_hbond == 1 ) 
      {
        type_j     = system->packed_atoms[j].type;
        if (type_j < 0) continue;
        start_j    = Start_Index(j, bonds);
        end_j      = End_Index(j, bonds);
        hb_start_j = Start_Index( system->Hindex[j], hbonds );
        hb_end_j   = End_Index  ( system->Hindex[j], hbonds );

        top = 0;
        for( pi = start_j; pi < end_j; ++pi )  
        {
          pbond_ij = &( bond_list[pi] );
          i = pbond_ij->nbr;
          type_i = system->packed_atoms[i].type;
	        if (type_i < 0) continue;
          bo_ij = &(bonds->bo_data_list[pi]);

          if(system->reax_param.sbp[type_i].p_hbond == 2 && BO_list[pi] >= HB_THRESHOLD)
            hblist[top++] = pi;
        }

        for( pk = hb_start_j; pk < hb_end_j; ++pk ) 
        {
          /* set k's varibles */
          k = hbond_list[pk].nbr;
          type_k = system->packed_atoms[k].type;
	        if (type_k < 0) continue;
          nbr_jk = hbond_list[pk].ptr;
          r_jk = nbr_jk->d;
          rvec_Scale( dvec_jk, hbond_list[pk].scl, nbr_jk->dvec );

          for( itr = 0; itr < top; ++itr ) 
          {
            pi = hblist[itr];
            pbond_ij = &( bonds->select.bond_list[pi] );
            i = pbond_ij->nbr;

            if( system->packed_atoms[i].orig_id != system->packed_atoms[k].orig_id ) 
            {
              bo_ij = &(bonds->bo_data_list[pi]);
              type_i = system->packed_atoms[i].type;
	            if (type_i < 0) continue;
                hbp = &(system->reax_param.hbp[ type_i ][ type_j ][ type_k ]);
	            if (hbp->r0_hb <= 0.0) continue;

              Calculate_Theta( pbond_ij->dvec, pbond_ij->d, dvec_jk, r_jk,
                               &theta, &cos_theta );
              /* the derivative of cos(theta) */
              Calculate_dCos_Theta( pbond_ij->dvec, pbond_ij->d, dvec_jk, r_jk,
                                    &dcos_theta_di, &dcos_theta_dj,
                                    &dcos_theta_dk );

              /* hyrogen bond energy*/
              sin_theta2 = sin( theta/2.0 );
              sin_xhz4 = SQR(sin_theta2);
              sin_xhz4 *= sin_xhz4;
              cos_xhz1 = ( 1.0 - cos_theta );
              exp_hb2 = exp( -hbp->p_hb2 * BO_list[pi]);
              exp_hb3 = exp( -hbp->p_hb3 * ( hbp->r0_hb / r_jk +
                                             r_jk / hbp->r0_hb - 2.0 ) );

              data->my_en.e_hb += e_hb =
                hbp->p_hb1 * (1.0 - exp_hb2) * exp_hb3 * sin_xhz4;

              CEhb1 = hbp->p_hb1 * hbp->p_hb2 * exp_hb2 * exp_hb3 * sin_xhz4;
              CEhb2 = -hbp->p_hb1/2.0 * (1.0 - exp_hb2) * exp_hb3 * cos_xhz1;
              CEhb3 = -hbp->p_hb3 *
                (-hbp->r0_hb / SQR(r_jk) + 1.0 / hbp->r0_hb) * e_hb;

              /* hydrogen bond forces */
              Cdbo_list[pi] += CEhb1; // dbo term

              if( control->virial == 0 ) 
              {
                // dcos terms
                rvec4_ScaledAdd( workspace->fCdDelta[i], +CEhb2, dcos_theta_di );
                rvec4_ScaledAdd( workspace->fCdDelta[j], +CEhb2, dcos_theta_dj );
                rvec4_ScaledAdd( workspace->fCdDelta[k], +CEhb2, dcos_theta_dk );
                // dr terms
                rvec4_ScaledAdd( workspace->fCdDelta[j], -CEhb3/r_jk, dvec_jk );
                rvec4_ScaledAdd( workspace->fCdDelta[k], +CEhb3/r_jk, dvec_jk );
              }
              else 
              {
                rvec_Scale( force, +CEhb2, dcos_theta_di ); // dcos terms
                rvec4_Add( workspace->fCdDelta[i], force );
                rvec_iMultiply( ext_press, pbond_ij->rel_box, force );
                rvec_ScaledAdd( data->my_ext_press, 1.0, ext_press );

                rvec4_ScaledAdd( workspace->fCdDelta[j], +CEhb2, dcos_theta_dj );

                ivec_Scale( rel_jk, hbond_list[pk].scl, nbr_jk->rel_box );
                rvec_Scale( force, +CEhb2, dcos_theta_dk );
                rvec4_Add( workspace->fCdDelta[k], force );
                rvec_iMultiply( ext_press, rel_jk, force );
                rvec_ScaledAdd( data->my_ext_press, 1.0, ext_press );
                // dr terms
                rvec4_ScaledAdd( workspace->fCdDelta[j], -CEhb3/r_jk, dvec_jk );

                rvec_Scale( force, CEhb3/r_jk, dvec_jk );
                rvec4_Add( workspace->fCdDelta[k], force );
                rvec_iMultiply( ext_press, rel_jk, force );
                rvec_ScaledAdd( data->my_ext_press, 1.0, ext_press );
              }

              /* tally into per-atom virials */
              if (system->pair_ptr->vflag_atom || system->pair_ptr->evflag) 
              {
                rvec_ScaledSum( delij, 1., system->packed_atoms[j].x,
                                      -1., system->packed_atoms[i].x );
                rvec_ScaledSum( delkj, 1., system->packed_atoms[j].x,
                                       -1.,system->packed_atoms[k].x );

                rvec_Scale(fi_tmp, CEhb2, dcos_theta_di);
                rvec_Scale(fk_tmp, CEhb2, dcos_theta_dk);
                rvec_ScaledAdd(fk_tmp, CEhb3/r_jk, dvec_jk);

                system->pair_ptr->ev_tally3(i,j,k,e_hb,0.0,fi_tmp,fk_tmp,delij,delkj);
              }
            }
          }
        }
      }//if
    }//for
  }
}
}
