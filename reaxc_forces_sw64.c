#include "sunway.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "reaxc_defs_sunway.h"
#include "reaxc_ctypes_sunway.h"
#include "reaxc_inlines_sw64.h"
#include "reaxc_forces_sw64.h"
#include "simd.h"
#include "sleef_math.h"
#define MAX_LEN 20

#define MAX(A,B) ((A) > (B) ? (A) : (B))
//void my_bubble_sort(bond_data *arr, int n)
//{
//  int i, j;
//  for(i = 0; i < n-1; i++)
//  {
//    for(j = 0; j < n-1-i; j++)
//    {
//      if(arr[j].nbr > arr[j+1].nbr)
//      {
//        bond_data tmp = arr[j];
//        arr[j] = arr[j+1];
//        arr[j+1] = tmp;
//      }
//    }
//  }
//}

void my_bubble_sort(bond_data *arr, double *Cdbo, double *Cdbopi, double *Cdbopi2, double *BO, rvec2 *BOpi, bond_order_data *bo_data, int n)
{
  int i, j;
  for(i = 0; i < n-1; i++)
  {
    for(j = 0; j < n-1-i; j++)
    {
      if(arr[j].nbr > arr[j+1].nbr)
      {
        bond_data tmp = arr[j];
        arr[j] = arr[j+1];
        arr[j+1] = tmp;
        
        double tcdbo = Cdbo[j];
        Cdbo[j] = Cdbo[j+1];
        Cdbo[j+1] = tcdbo;
        
        double tcdbopi  = Cdbopi[j];
        Cdbopi[j]     = Cdbopi[j+1];
        Cdbopi[j+1]   = tcdbopi;
        
        double tcdbopi2 = Cdbopi2[j];
        Cdbopi2[j]    = Cdbopi2[j+1];
        Cdbopi2[j+1]  = tcdbopi2;
        
        double tbo = BO[j];
        BO[j]    = BO[j+1];
        BO[j+1]  = tbo;
        
        double tbopi = BOpi[j][0];
        BOpi[j][0]    = BOpi[j+1][0];
        BOpi[j+1][0]  = tbopi;

        double tbopi2 = BOpi[j][1];
        BOpi[j][1]    = BOpi[j+1][1];
        BOpi[j+1][1]  = tbopi2;

        bond_order_data tbo_data = bo_data[j];
        bo_data[j] = bo_data[j+1];
        bo_data[j+1] = tbo_data;

      }
    }
  }
}

void Add_All_dBond_to_Forces_C_org( reax_system_c *system, control_params *control,
                                simulation_data *data, storage *workspace,
                                reax_list **lists, mpi_datatypes *mpi_data)
{
  int i, pj;
  reax_list *bonds = (*lists) + BONDS;
  double *Cdbo_list = bonds->Cdbo_list;
  double *Cdbopi_list = bonds->Cdbopi_list;
  double *Cdbopi2_list = bonds->Cdbopi2_list;
  
  bond_order_data *bo_t;

  printf("in Add dBond Forces C org\n");

  for( i = 0; i < system->N; ++i )
  {
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
    {
      if( i < bonds->select.bond_list[pj].nbr ) 
      {

        bond_data *nbr_j, *nbr_k;
        bond_order_data *bo_ij, *bo_ji;
        dbond_coefficients coef;
        int pk, k, j;

        /* Virial Tallying variables */
        rvec fi_tmp, fj_tmp, fk_tmp, delij, delji, delki, delkj, temp;

        /* Initializations */
        nbr_j = &(bonds->select.bond_list[pj]);
        j = nbr_j->nbr;
        //bo_ij = &(nbr_j->bo_data);
        bo_ij = &(bonds->bo_data_list[pj]);
        //bo_ji = &(bonds->select.bond_list[ nbr_j->sym_index ].bo_data);
        bo_ji = &(bonds->bo_data_list[ nbr_j->sym_index]);

        //coef.C1dbo = bo_ij->C1dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
        //coef.C2dbo = bo_ij->C2dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
        //coef.C3dbo = bo_ij->C3dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
        coef.C1dbo = bo_ij->C1dbo * (Cdbo_list[pj] + Cdbo_list[nbr_j->sym_index]);
        coef.C2dbo = bo_ij->C2dbo * (Cdbo_list[pj] + Cdbo_list[nbr_j->sym_index]);
        coef.C3dbo = bo_ij->C3dbo * (Cdbo_list[pj] + Cdbo_list[nbr_j->sym_index]);


        //coef.C1dbopi = bo_ij->C1dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
        //coef.C2dbopi = bo_ij->C2dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
        //coef.C3dbopi = bo_ij->C3dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
        //coef.C4dbopi = bo_ij->C4dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
        coef.C1dbopi = bo_ij->C1dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);
        coef.C2dbopi = bo_ij->C2dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);
        coef.C3dbopi = bo_ij->C3dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);
        coef.C4dbopi = bo_ij->C4dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);


        //coef.C1dbopi2 = bo_ij->C1dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
        //coef.C2dbopi2 = bo_ij->C2dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
        //coef.C3dbopi2 = bo_ij->C3dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
        //coef.C4dbopi2 = bo_ij->C4dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
        coef.C1dbopi2 = bo_ij->C1dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
        coef.C2dbopi2 = bo_ij->C2dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
        coef.C3dbopi2 = bo_ij->C3dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
        coef.C4dbopi2 = bo_ij->C4dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);


        coef.C1dDelta = bo_ij->C1dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);
        coef.C2dDelta = bo_ij->C2dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);
        coef.C3dDelta = bo_ij->C3dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);

        // forces on i
        rvec_Scale(           temp, coef.C1dbo,    bo_ij->dBOp );
        rvec_ScaledAdd( temp, coef.C2dbo,    workspace->dDeltap_self[i] );
        rvec_ScaledAdd( temp, coef.C1dDelta, bo_ij->dBOp );
        rvec_ScaledAdd( temp, coef.C2dDelta, workspace->dDeltap_self[i] );
        rvec_ScaledAdd( temp, coef.C1dbopi,  bo_ij->dln_BOp_pi );
        rvec_ScaledAdd( temp, coef.C2dbopi,  bo_ij->dBOp );
        rvec_ScaledAdd( temp, coef.C3dbopi,  workspace->dDeltap_self[i]);
        rvec_ScaledAdd( temp, coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
        rvec_ScaledAdd( temp, coef.C2dbopi2, bo_ij->dBOp );
        rvec_ScaledAdd( temp, coef.C3dbopi2, workspace->dDeltap_self[i] );
        rvec4_Add( workspace->fCdDelta[i], temp );

        if (system->vflag_atom)
        {
          rvec_Scale(fi_tmp, -1.0, temp);
          rvec_ScaledSum( delij, 1., system->packed_atoms[i].x,-1., system->packed_atoms[j].x );
          v_tally_sys(system, i, fi_tmp, delij);
        }

        // forces on j
        rvec_Scale(           temp, -coef.C1dbo,    bo_ij->dBOp );
        rvec_ScaledAdd( temp,  coef.C3dbo,    workspace->dDeltap_self[j] );
        rvec_ScaledAdd( temp, -coef.C1dDelta, bo_ij->dBOp );
        rvec_ScaledAdd( temp,  coef.C3dDelta, workspace->dDeltap_self[j]);
        rvec_ScaledAdd( temp, -coef.C1dbopi,  bo_ij->dln_BOp_pi );
        rvec_ScaledAdd( temp, -coef.C2dbopi,  bo_ij->dBOp );
        rvec_ScaledAdd( temp,  coef.C4dbopi,  workspace->dDeltap_self[j]);
        rvec_ScaledAdd( temp, -coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
        rvec_ScaledAdd( temp, -coef.C2dbopi2, bo_ij->dBOp );
        rvec_ScaledAdd( temp,  coef.C4dbopi2, workspace->dDeltap_self[j]);
        rvec4_Add( workspace->fCdDelta[j], temp );

        if (system->vflag_atom)
        {
          rvec_Scale(fj_tmp, -1.0, temp);
          rvec_ScaledSum( delji, 1., system->packed_atoms[j].x,-1., system->packed_atoms[i].x );
          v_tally_sys(system, j, fj_tmp, delji);
        }

        // forces on k: i neighbor
        for( pk = Start_Index(i, bonds); pk < End_Index(i, bonds); ++pk ) 
        {
          nbr_k = &(bonds->select.bond_list[pk]);
          k = nbr_k->nbr;
          bo_t = &(bonds->bo_data_list[pk]);

          //rvec_Scale(     temp, -coef.C2dbo,    nbr_k->bo_data.dBOp);
          //rvec_ScaledAdd( temp, -coef.C2dDelta, nbr_k->bo_data.dBOp);
          //rvec_ScaledAdd( temp, -coef.C3dbopi,  nbr_k->bo_data.dBOp);
          //rvec_ScaledAdd( temp, -coef.C3dbopi2, nbr_k->bo_data.dBOp);
          rvec_Scale(     temp, -coef.C2dbo,    bo_t->dBOp);
          rvec_ScaledAdd( temp, -coef.C2dDelta, bo_t->dBOp);
          rvec_ScaledAdd( temp, -coef.C3dbopi,  bo_t->dBOp);
          rvec_ScaledAdd( temp, -coef.C3dbopi2, bo_t->dBOp);

          rvec4_Add( workspace->fCdDelta[k], temp );

          if( system->vflag_atom ) 
          {
            rvec_Scale(fk_tmp, -1.0, temp);
            rvec_ScaledSum(delki,1.,system->packed_atoms[k].x,-1.,system->packed_atoms[i].x);
            //system->pair_ptr->v_tally(k,fk_tmp,delki);
            v_tally_sys(system, k, fk_tmp, delki);
            rvec_ScaledSum(delkj,1.,system->packed_atoms[k].x,-1.,system->packed_atoms[j].x);
            v_tally_sys(system, k, fk_tmp, delkj);
            //system->pair_ptr->v_tally(k,fk_tmp,delkj);
          }
        }

        // forces on k: j neighbor
        for( pk = Start_Index(j, bonds); pk < End_Index(j, bonds); ++pk ) 
        {
          nbr_k = &(bonds->select.bond_list[pk]);
          k = nbr_k->nbr;
          bo_t = &(bonds->bo_data_list[pk]);

          rvec_Scale(     temp, -coef.C3dbo,    bo_t->dBOp );
          rvec_ScaledAdd( temp, -coef.C3dDelta, bo_t->dBOp);
          rvec_ScaledAdd( temp, -coef.C4dbopi,  bo_t->dBOp);
          rvec_ScaledAdd( temp, -coef.C4dbopi2, bo_t->dBOp);
          rvec4_Add( workspace->fCdDelta[k], temp );

          if( system->vflag_atom ) 
          {
            rvec_Scale(fk_tmp, -1.0, temp);
            rvec_ScaledSum(delki,1.,system->packed_atoms[k].x,-1.,system->packed_atoms[i].x);
            v_tally_sys(system, k, fk_tmp, delki);
            rvec_ScaledSum(delkj,1.,system->packed_atoms[k].x,-1.,system->packed_atoms[j].x);
            v_tally_sys(system, k, fk_tmp, delkj);
          }
        }
      }
    }
  }
}


///*****no if branch*****/
//void Add_All_dBond_to_Forces_C( reax_system_c *system, control_params *control,
//                                simulation_data *data, storage *workspace,
//                                reax_list **lists, mpi_datatypes *mpi_data)
//{
//  printf("in Add dBond Forces C\n");
//  int i, pj, pj_off;
//  reax_list *bonds = (*lists) + BONDS;
//  double *Cdbo_list = bonds->Cdbo_list;
//  double *Cdbopi_list = bonds->Cdbopi_list;
//  double *Cdbopi2_list = bonds->Cdbopi2_list;
//  int start_i, end_i, len_i;
//  bond_data *nbr_j, *nbr_k;
//  int pk, k, j, pk_off;
//  rvec temp;
//  double buf_c;
//  int buf_idx[MAX_BONDS];
//  rvec buf_dBOp[MAX_BONDS];
//  bond_order_data *bo_t;
//
//  for( i = 0; i < system->N; ++i )
//  {
//    start_i = Start_Index(i, bonds);
//    end_i   = End_Index(i, bonds);
//    len_i   = end_i - start_i;
//    buf_c = 0.0;
//    printf("i = %d\n", i);
//    for( pj = start_i; pj < end_i; ++pj )
//    {
//      //printf("i = %d, pj = %d\n", i, pj);
//      nbr_j = &(bonds->select.bond_list[pj]);
//      j = nbr_j->nbr;
//      pj_off = pj - start_i;
//
//      bond_order_data *bo_ij, *bo_ji;
//      dbond_coefficients coef;
//
//      /* Virial Tallying variables */
//      rvec fi_tmp, fj_tmp, fk_tmp, delij, delji, delki, delkj;
//
//      /* Initializations */
//      //bo_ij = &(nbr_j->bo_data);
//      bo_ij = &(bonds->bo_data_list[pj]);
//      //bo_ji = &(bonds->select.bond_list[ nbr_j->sym_index ].bo_data);
//      bo_ji = &(bonds->bo_data_list[ nbr_j->sym_index]);
//
//      //printf("-------00--------\n");
//
//      //coef.C1dbo = bo_ij->C1dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
//      //coef.C2dbo = bo_ij->C2dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
//      //coef.C3dbo = bo_ij->C3dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
//      coef.C1dbo = bo_ij->C1dbo * (Cdbo_list[pj] + Cdbo_list[nbr_j->sym_index]);
//      coef.C2dbo = bo_ij->C2dbo * (Cdbo_list[pj] + Cdbo_list[nbr_j->sym_index]);
//      coef.C3dbo = bo_ij->C3dbo * (Cdbo_list[pj] + Cdbo_list[nbr_j->sym_index]);
//
//      //printf("-------01--------\n");
//
//      //coef.C1dbopi = bo_ij->C1dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
//      //coef.C2dbopi = bo_ij->C2dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
//      //coef.C3dbopi = bo_ij->C3dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
//      //coef.C4dbopi = bo_ij->C4dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
//      coef.C1dbopi = bo_ij->C1dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);
//      coef.C2dbopi = bo_ij->C2dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);
//      coef.C3dbopi = bo_ij->C3dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);
//      coef.C4dbopi = bo_ij->C4dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);
//
//      //printf("-------02--------\n");
//
//      //coef.C1dbopi2 = bo_ij->C1dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
//      //coef.C2dbopi2 = bo_ij->C2dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
//      //coef.C3dbopi2 = bo_ij->C3dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
//      //coef.C4dbopi2 = bo_ij->C4dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
//      
//      //coef.C1dbopi2 = 0;//bo_ij->C1dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//      //coef.C2dbopi2 = 0;//bo_ij->C2dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//      //coef.C3dbopi2 = 0;//bo_ij->C3dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//      //coef.C4dbopi2 = 0;//bo_ij->C4dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//
//
//  
//    
//      if(i == 10822)
//      {
//      printf("C1: %d, %d, %d, %d\n", isNumber(bo_ij->C1dbopi2), isNumber(Cdbopi2_list[pj]), isNumber(Cdbopi2_list[nbr_j->sym_index]), isNumber(coef.C1dbopi2));
//      coef.C1dbopi2 = bo_ij->C1dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//      printf("C2: %d, %d, %d, %d\n", isNumber(bo_ij->C2dbopi2), isNumber(Cdbopi2_list[pj]), isNumber(Cdbopi2_list[nbr_j->sym_index]), isNumber(coef.C2dbopi2));
//      coef.C2dbopi2 = bo_ij->C2dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//      //printf("C3: %d, %d, %d, %d\n", isNumber(bo_ij->C3dbopi2), isNumber(Cdbopi2_list[pj]), isNumber(Cdbopi2_list[nbr_j->sym_index]), isNumber(coef.C3dbopi2));
//      //printf("C3: %d, %d, %d, %d\n", isnan(bo_ij->C3dbopi2), isnan(Cdbopi2_list[pj]), isnan(Cdbopi2_list[nbr_j->sym_index]), isnan(coef.C3dbopi2));
//      //printf("C31:  %d\n", isnan(bo_ij->C3dbopi2));
//      //printf("C32:  %d\n", isnan(Cdbopi2_list[pj]));
//      //printf("C33:  %d\n", isnan(Cdbopi2_list[nbr_j->sym_index]));
//      //printf("C34:  %d\n", isnan(coef.C3dbopi2));
//      
//      printf("C31:  %d\n", isinf(bo_ij->C3dbopi2));
//      printf("C32:  %d\n", isinf(Cdbopi2_list[pj]));
//      printf("C33:  %d\n", isinf(Cdbopi2_list[nbr_j->sym_index]));
//      printf("C34:  %d\n", isinf(coef.C3dbopi2));
//
//
//      coef.C3dbopi2 = bo_ij->C3dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//      printf("C4: %d, %d, %d, %d\n", isNumber(bo_ij->C4dbopi2), isNumber(Cdbopi2_list[pj]), isNumber(Cdbopi2_list[nbr_j->sym_index]), isNumber(coef.C4dbopi2));
//
//      coef.C4dbopi2 = bo_ij->C4dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//      }
//      else
//      {
//          coef.C1dbopi2 = bo_ij->C1dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//          coef.C2dbopi2 = bo_ij->C2dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//          coef.C3dbopi2 = bo_ij->C3dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//          coef.C4dbopi2 = bo_ij->C4dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
//      }
//
//
//      coef.C1dDelta = bo_ij->C1dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);
//      coef.C2dDelta = bo_ij->C2dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);
//      coef.C3dDelta = bo_ij->C3dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);
//
//
//
//      // forces on i
//      rvec_Scale(           temp, coef.C1dbo,    bo_ij->dBOp );
//      rvec_ScaledAdd( temp, coef.C2dbo,    workspace->dDeltap_self[i] );
//      rvec_ScaledAdd( temp, coef.C1dDelta, bo_ij->dBOp );
//      rvec_ScaledAdd( temp, coef.C2dDelta, workspace->dDeltap_self[i] );
//      rvec_ScaledAdd( temp, coef.C1dbopi,  bo_ij->dln_BOp_pi );
//      rvec_ScaledAdd( temp, coef.C2dbopi,  bo_ij->dBOp );
//      rvec_ScaledAdd( temp, coef.C3dbopi,  workspace->dDeltap_self[i]);
//      rvec_ScaledAdd( temp, coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
//      rvec_ScaledAdd( temp, coef.C2dbopi2, bo_ij->dBOp );
//      rvec_ScaledAdd( temp, coef.C3dbopi2, workspace->dDeltap_self[i] );
//      rvec4_Add( workspace->fCdDelta[i], temp );
//      
//      //printf("-------2--------\n");
//      buf_c += (-coef.C2dbo - coef.C2dDelta - coef.C3dbopi - coef.C3dbopi2);
//      //rvec_Copy(buf_dBOp[pj_off], nbr_j->bo_data.dBOp);
//      rvec_Copy(buf_dBOp[pj_off], bonds->bo_data_list[pj].dBOp);
//      buf_idx[pj_off] = j;
//
//      //printf("-------3--------\n");
//    }//for-pj
//   
//    //printf("-------4--------\n");
//    // forces on k: i neighbor
//    for( pj = start_i; pj < end_i; ++pj) 
//    {
//      pj_off = pj - start_i;
//      rvec_Scale(     temp, buf_c, buf_dBOp[pj_off]);
//      rvec4_Add( workspace->fCdDelta[buf_idx[pj_off]], temp);
//    }
//
//    //printf("-------5--------\n");
//  }//for-i
//  printf("End Add dBond Forces C\n");
//}



#ifdef MPE
#include <mpi.h>
#include <athread.h>
//#define LWPF_UNITS U(ADD_DBOND)
//#include "lwpf2.h"
extern SLAVE_FUN(Init_Forces_noQEq_Full_C_para)(init_forces_pack_t *);
extern SLAVE_FUN(Add_All_dBond_to_Forces_C_Para)(add_dbond_pack_t *);
extern SLAVE_FUN(Init_Forces_noQEq_HB_Full_C_para)(init_forces_pack_t *);

int r = 0;

//swcache_lock_t *locks_dbond = NULL;
//swcache_lock_t *locks_frc = NULL;

/*****no if branch*****/
void Add_All_dBond_to_Forces_C( reax_system_c *system, control_params *control,
                                simulation_data *data, storage *workspace,
                                reax_list **lists, mpi_datatypes *mpi_data)
{
  //printf("in Add dBond Forces C\n");
  int i, pj, pj_off;
  reax_list *bonds = (*lists) + BONDS;
  double *Cdbo_list = bonds->Cdbo_list;
  double *Cdbopi_list = bonds->Cdbopi_list;
  double *Cdbopi2_list = bonds->Cdbopi2_list;
  int start_i, end_i, len_i;
  bond_data *nbr_j, *nbr_k;
  int pk, k, j, pk_off;
  rvec temp;
  double buf_c;
  int buf_idx[MAX_BONDS];
  rvec buf_dBOp[MAX_BONDS];
  bond_order_data *bo_t;

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //perf_config_t conf;
  //conf.pcrc = PCRC_ALL;
  //conf.pcr0 = PC0_CYCLE;
  //conf.pcr1 = PC1_CYCLE;
  //conf.pcr2 = PC2_N_GLD;
  ////conf.pcr2 = PC2_N_GF_AND_A;
  ////conf.pcr2 = PC2_N_DMA_REQ;
  //lwpf_init(&conf);
 
  #define ADD_DBOND
  #ifdef  ADD_DBOND
  
  add_dbond_pack_t pm;
  pm.system       = system;
  pm.control      = control;
  pm.data         = data;
  pm.workspace    = workspace;
  pm.lists        = lists;
  pm.fCdDelta     = workspace->fCdDelta;
  pm.dDeltap_self = workspace->dDeltap_self;
  pm.Cdbo_list    = bonds->Cdbo_list;
  pm.Cdbopi_list  = bonds->Cdbopi_list;
  pm.Cdbopi2_list = bonds->Cdbopi2_list;
  pm.N            = system->N;
  pm.index        = bonds->index;
  pm.end_index    = bonds->end_index;
  pm.bond_list    = bonds->select.bond_list;
  pm.bo_data_list = bonds->bo_data_list;
  
  int nall = system->total_cap;
  //if(locks_dbond == NULL)
  //  locks_dbond = swcache_u_prepare_locks(nall, 3);
  //pm.locks_dbond  = locks_dbond;

  //if(locks_frc == NULL)
  //  locks_frc = swcache_u_prepare_locks(nall, 3);
  //pm.locks_frc  = locks_frc;
  pm.locks_frc  = bonds->locks_frc;


  if(athread_idle() == 0)
    athread_init();
  athread_spawn(Add_All_dBond_to_Forces_C_Para, &pm);
  athread_join();
  //if(locks_dbond != NULL) 
  //{
  //  free(locks_dbond); 
  //  locks_dbond = NULL;
  //}

  //if(locks_frc != NULL) 
  //{
  //  free(locks_frc); 
  //  locks_frc = NULL;
  //}
  //if(myrank == 0)
  //{
  //  lwpf_report_summary(stdout, &conf);
  //}


  #else
  for( i = 0; i < system->N; ++i )
  {
    start_i = Start_Index(i, bonds);
    end_i   = End_Index(i, bonds);
    len_i   = end_i - start_i;
    buf_c = 0.0;
    for( pj = start_i; pj < end_i; ++pj )
    {
      //printf("i = %d, pj = %d\n", i, pj);
      nbr_j = &(bonds->select.bond_list[pj]);
      j = nbr_j->nbr;
      pj_off = pj - start_i;

      bond_order_data *bo_ij, *bo_ji;
      dbond_coefficients coef;

      /* Virial Tallying variables */
      rvec fi_tmp, fj_tmp, fk_tmp, delij, delji, delki, delkj;

      /* Initializations */
      //bo_ij = &(nbr_j->bo_data);
      bo_ij = &(bonds->bo_data_list[pj]);
      //bo_ji = &(bonds->select.bond_list[ nbr_j->sym_index ].bo_data);
      bo_ji = &(bonds->bo_data_list[ nbr_j->sym_index]);

      //coef.C1dbo = bo_ij->C1dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
      //coef.C2dbo = bo_ij->C2dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
      //coef.C3dbo = bo_ij->C3dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
      coef.C1dbo = bo_ij->C1dbo * (Cdbo_list[pj] + Cdbo_list[nbr_j->sym_index]);
      coef.C2dbo = bo_ij->C2dbo * (Cdbo_list[pj] + Cdbo_list[nbr_j->sym_index]);
      coef.C3dbo = bo_ij->C3dbo * (Cdbo_list[pj] + Cdbo_list[nbr_j->sym_index]);

      //coef.C1dbopi = bo_ij->C1dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
      //coef.C2dbopi = bo_ij->C2dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
      //coef.C3dbopi = bo_ij->C3dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
      //coef.C4dbopi = bo_ij->C4dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
      coef.C1dbopi = bo_ij->C1dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);
      coef.C2dbopi = bo_ij->C2dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);
      coef.C3dbopi = bo_ij->C3dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);
      coef.C4dbopi = bo_ij->C4dbopi * (Cdbopi_list[pj] + Cdbopi_list[nbr_j->sym_index]);

      //coef.C1dbopi2 = bo_ij->C1dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
      //coef.C2dbopi2 = bo_ij->C2dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
      //coef.C3dbopi2 = bo_ij->C3dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
      //coef.C4dbopi2 = bo_ij->C4dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
      
      //coef.C1dbopi2 = 0;//bo_ij->C1dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
      //coef.C2dbopi2 = 0;//bo_ij->C2dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
      //coef.C3dbopi2 = 0;//bo_ij->C3dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
      //coef.C4dbopi2 = 0;//bo_ij->C4dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);

      coef.C1dbopi2 = bo_ij->C1dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
      coef.C2dbopi2 = bo_ij->C2dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
      coef.C3dbopi2 = bo_ij->C3dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);
      coef.C4dbopi2 = bo_ij->C4dbopi2 * (Cdbopi2_list[pj] + Cdbopi2_list[nbr_j->sym_index]);

      coef.C1dDelta = bo_ij->C1dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);
      coef.C2dDelta = bo_ij->C2dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);
      coef.C3dDelta = bo_ij->C3dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);

      // forces on i
      rvec_Scale(           temp, coef.C1dbo,    bo_ij->dBOp );
      rvec_ScaledAdd( temp, coef.C2dbo,    workspace->dDeltap_self[i] );
      rvec_ScaledAdd( temp, coef.C1dDelta, bo_ij->dBOp );
      rvec_ScaledAdd( temp, coef.C2dDelta, workspace->dDeltap_self[i] );
      rvec_ScaledAdd( temp, coef.C1dbopi,  bo_ij->dln_BOp_pi );
      rvec_ScaledAdd( temp, coef.C2dbopi,  bo_ij->dBOp );
      rvec_ScaledAdd( temp, coef.C3dbopi,  workspace->dDeltap_self[i]);
      rvec_ScaledAdd( temp, coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
      rvec_ScaledAdd( temp, coef.C2dbopi2, bo_ij->dBOp );
      rvec_ScaledAdd( temp, coef.C3dbopi2, workspace->dDeltap_self[i] );
      rvec4_Add( workspace->fCdDelta[i], temp );
      
      buf_c += (-coef.C2dbo - coef.C2dDelta - coef.C3dbopi - coef.C3dbopi2);
      //rvec_Copy(buf_dBOp[pj_off], nbr_j->bo_data.dBOp);
      rvec_Copy(buf_dBOp[pj_off], bonds->bo_data_list[pj].dBOp);
      buf_idx[pj_off] = j;
    }//for-pj
   
    // forces on k: i neighbor
    for( pj = start_i; pj < end_i; ++pj) 
    {
      pj_off = pj - start_i;
      rvec_Scale(     temp, buf_c, buf_dBOp[pj_off]);
      rvec4_Add( workspace->fCdDelta[buf_idx[pj_off]], temp);
    }

  }//for-i
  #endif
}



void Init_Forces_noQEq_HB_Full_C(init_forces_pack_t * param);
void Init_Forces_noQEq_Full_C(init_forces_pack_t * param)
{
  if(athread_idle() == 0)
    athread_init();
  
  reax_system_c *system       = param->system;
  control_params *control     = param->control;
  simulation_data *data       = param->data;
  storage *workspace          = param->workspace;
  reax_list **lists           = param->lists;

  reax_list *far_nbrs = *lists + FAR_NBRS_FULL;
  reax_list *bonds    = *lists + BONDS;
  reax_list *hbonds   = *lists + HBONDS;
  param->bond_list            = bonds->select.bond_list;
  param->bindex               = bonds->index;
  param->bend_index           = bonds->end_index;
  param->far_nbr_list_full    = far_nbrs->select.far_nbr_list_full;
  param->findex               = far_nbrs->index;
  param->fend_index           = far_nbrs->end_index;
  param->hbond_list           = hbonds->select.hbond_list;
  param->hindex               = hbonds->index;
  param->hend_index           = hbonds->end_index;


  param->n                    = system->n;
  param->N                    = system->N;
  param->num_atom_types       = system->reax_param.num_atom_types;
  param->bo_cut               = control->bo_cut;
  param->bond_cut             = control->bond_cut;
  param->hbond_cut            = control->hbond_cut;
  param->num_intrs            = bonds->num_intrs;
  param->bond_mark            = workspace->bond_mark;
  param->dDeltap_self         = workspace->dDeltap_self;
  //param->total_bond_order     = workspace->total_bond_order;
  param->num_bonds            = &(workspace->realloc.num_bonds);
  param->num_hbonds           = &(workspace->realloc.num_hbonds);
  param->bo_dboc              = workspace->bo_dboc;
  param->Cdbo_list            = bonds->Cdbo_list;
  param->Cdbopi_list          = bonds->Cdbopi_list;
  param->Cdbopi2_list         = bonds->Cdbopi2_list;
  param->BO_list              = bonds->BO_list;
  param->BOpi_list            = bonds->BOpi_list;
  param->bo_data_list         = bonds->bo_data_list;
 
  int ntypes = system->reax_param.num_atom_types;
  int i, j;
  init_forces_tbp pack_tbp[ntypes][ntypes];
  init_forces_sbp pack_sbp[ntypes];
  two_body_parameters **tbpr = system->reax_param.tbp;
  single_body_parameters *sbpr = system->reax_param.sbp;

  for (i = 0; i < ntypes; i++)//ntypes=4;
  {
    for (j = 0; j < ntypes; j++)
    {
      pack_tbp[i][j].p_bo1  = tbpr[i][j].p_bo1;
      pack_tbp[i][j].p_bo2  = tbpr[i][j].p_bo2;
      pack_tbp[i][j].p_bo3  = tbpr[i][j].p_bo3;
      pack_tbp[i][j].p_bo4  = tbpr[i][j].p_bo4;
      pack_tbp[i][j].p_bo5  = tbpr[i][j].p_bo5;
      pack_tbp[i][j].p_bo6  = tbpr[i][j].p_bo6;
      pack_tbp[i][j].r_s    = tbpr[i][j].r_s;
      pack_tbp[i][j].r_p    = tbpr[i][j].r_p;
      pack_tbp[i][j].r_pp   = tbpr[i][j].r_pp;
    }
    pack_sbp[i].r_s         = sbpr[i].r_s; 
    pack_sbp[i].r_pi        = sbpr[i].r_pi; 
    pack_sbp[i].r_pi_pi     = sbpr[i].r_pi_pi; 
    pack_sbp[i].p_hbond     = sbpr[i].p_hbond; 
    pack_sbp[i].padding     = 0; 
  }

  param->packed_atoms         = system->packed_atoms;
  param->pack_sbp             = pack_sbp;
  param->pack_tbp             = pack_tbp;
  
  //GPTLstart("Forces 1");
  athread_spawn(Init_Forces_noQEq_Full_C_para, param);
  Init_Forces_noQEq_HB_Full_C(param);
  athread_join();
  //GPTLstop("Forces 1");

  ////GPTLstart("Forces 2");
  ////Init_Forces_noQEq_HB_Full_C(param);
  //athread_spawn(Init_Forces_noQEq_HB_Full_C_para, param);
  //athread_join();
  ////GPTLstop("Forces 2");
}  
 

void Init_Forces_noQEq_HB_Full_C_org(init_forces_pack_t * param) 
{
  reax_system_c *system       = param->system;
  control_params *control     = param->control;
  simulation_data *data       = param->data;
  storage *workspace          = param->workspace;
  reax_list **lists           = param->lists;

  int i, j, pj;
  int start_i, end_i;
  int type_i, type_j;
  int num_hbonds;
  int ihb, jhb, ihb_top;
  int local, flag, renbr;
  double cutoff;
  reax_list *far_nbrs, *hbonds;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  far_neighbor_data_full *nbr_pj;
  reax_atom *atom_i, *atom_j;

  far_nbrs = *lists + FAR_NBRS_FULL;
  hbonds = *lists + HBONDS;

  num_hbonds = 0;
  renbr = (data->step-data->prev_steps) % control->reneighbor == 0;
  if (control->hbond_cut <= 0)
    return;
    
  //printf("system->n = %d\n", system->n);

  for(i = 0; i < system->n; ++i) 
  {
    //printf("i = %d\n", i);

    atom_i = &(system->my_atoms[i]);
    type_i  = system->packed_atoms[i].type;
    if (type_i < 0) continue;
    start_i = far_nbrs->index[i];
    end_i   = far_nbrs->end_index[i];
    sbp_i = &(system->reax_param.sbp[type_i]);

    local = 1;
    cutoff = control->hbond_cut;

    ihb = -1;
    ihb_top = -1;
    ihb = sbp_i->p_hbond;

    if(ihb == 1)
    {
      //printf("atom_i->Hindex=%d\n", atom_i->Hindex);
      //ihb_top = hbonds->end_index[atom_i->Hindex];
      ihb_top = hbonds->end_index[system->Hindex[i]];
    }
    else continue;;
    
      
    for(pj = start_i; pj < end_i; ++pj) 
    {
      //printf("pj = %d\n", pj);

      nbr_pj = &(far_nbrs->select.far_nbr_list_full[pj]);
      j = nbr_pj->nbr;
      atom_j = &(system->my_atoms[j]);

      if(nbr_pj->d <= cutoff)  flag = 1;
      else  flag = 0;

      if(flag) 
      {
        type_j = system->packed_atoms[j].type;
        if (type_j < 0) continue;
        sbp_j = &(system->reax_param.sbp[type_j]);
        twbp = &(system->reax_param.tbp[type_i][type_j]);

        if(local) 
        {
          if((ihb==1 || ihb==2)) 
          {
            jhb = sbp_j->p_hbond;
            if(ihb == 1 && jhb == 2) 
            {
              hbonds->select.hbond_list[ihb_top].nbr = j;
              hbonds->select.hbond_list[ihb_top].scl = 1;
              hbonds->select.hbond_list[ihb_top].ptr = nbr_pj;
              ++ihb_top;
              ++num_hbonds;
            }//if
          }//if
        }//if-local
      }//if-flag
    }//for-pj
    if(ihb == 1 )
    {
      //hbonds->end_index[atom_i->Hindex] = ihb_top;
      hbonds->end_index[system->Hindex[i]] = ihb_top;
    }
  }//for-i
  workspace->realloc.num_hbonds = num_hbonds;
}

void Init_Forces_noQEq_HB_Full_C(init_forces_pack_t * param) 
{
  reax_system_c *system       = param->system;
  control_params *control     = param->control;
  simulation_data *data       = param->data;
  storage *workspace          = param->workspace;
  reax_list **lists           = param->lists;

  int i, j, pj;
  int start_i, end_i;
  int type_i, type_j;
  int num_hbonds;
  int ihb, jhb, ihb_top;
  int local, flag, renbr;
  double cutoff;
  reax_list *far_nbrs, *hbonds;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  far_neighbor_data_full *nbr_pj;
  reax_atom *atom_i, *atom_j;

  far_nbrs = *lists + FAR_NBRS_FULL;
  hbonds = *lists + HBONDS;

  num_hbonds = 0;
  //renbr = (data->step-data->prev_steps) % control->reneighbor == 0;
  if (control->hbond_cut <= 0)
    return;
    
  for(i = 0; i < system->n; ++i) 
  {
    //atom_i = &(system->my_atoms[i]);
    type_i  = system->packed_atoms[i].type;
    if (type_i < 0) continue;
    start_i = far_nbrs->index[i];
    end_i   = far_nbrs->end_index[i];
    sbp_i = &(system->reax_param.sbp[type_i]);

    cutoff = control->hbond_cut;
    ihb = -1;
    ihb_top = -1;
    ihb = sbp_i->p_hbond;
      
    if(ihb == 1)
    {
      ihb_top = hbonds->end_index[system->Hindex[i]];
      for(pj = start_i; pj < end_i; ++pj) 
      {
        nbr_pj = &(far_nbrs->select.far_nbr_list_full[pj]);
        //atom_j = &(system->my_atoms[j]);

        if(nbr_pj->d <= cutoff)
        {
          j = nbr_pj->nbr;
          type_j = system->packed_atoms[j].type;
          if (type_j < 0) continue;
          sbp_j = &(system->reax_param.sbp[type_j]);
          //twbp = &(system->reax_param.tbp[type_i][type_j]);
          if((ihb==1 || ihb==2)) 
          {
            jhb = sbp_j->p_hbond;
            if(ihb == 1 && jhb == 2) 
            {
              hbonds->select.hbond_list[ihb_top].nbr = j;
              hbonds->select.hbond_list[ihb_top].scl = 1;
              hbonds->select.hbond_list[ihb_top].ptr = nbr_pj;
              ++ihb_top;
              ++num_hbonds;
            }//if
          }//if
        }//if-cutoff
      }//for-pj
      hbonds->end_index[system->Hindex[i]] = ihb_top;
    }
  }//for-i
  workspace->realloc.num_hbonds = num_hbonds;
}


void Validate_Lists_C(init_forces_pack_t * param)
{
  reax_system_c *system    = param->system;
  storage *workspace     = param->workspace;
  reax_list **lists      = param->lists;
  simulation_data * data = param->data;
  int step = data->step;
  int n    = system->n;
  int N    = system->N;
  int numH = system->numH;
  MPI_Comm comm = param->comm;

  int i, comp, Hindex;
  reax_list *bonds, *hbonds;
  double saferzone = system->saferzone;
  /* bond list */
  if( N > 0 ) 
  {
    bonds = *lists + BONDS;
    for( i = 0; i < N; ++i ) 
    {
      //system->my_atoms[i].num_bonds = 
      system->num_bonds[i] = 
                    MAX((bonds->end_index[i] - bonds->index[i])*2, MIN_BONDS);
      if( i < N-1 ) comp = bonds->index[i+1];
      else comp = bonds->num_intrs;
      if( bonds->end_index[i] > comp ) 
      {
        fprintf( stderr, "step%d-bondchk failed: i=%d end(i)=%d str(i+1)=%d\n",
                 step, i, bonds->end_index[i], comp );
        MPI_Abort( comm, INSUFFICIENT_MEMORY );
      }
    }//for
  }//if

  /* hbonds list */
  if( numH > 0 )
  {
    hbonds = *lists + HBONDS;
    for( i = 0; i < N; ++i ) 
    {
      //Hindex = system->my_atoms[i].Hindex;
      Hindex = system->Hindex[i];

      if( Hindex > -1 ) 
      {
        //system->my_atoms[i].num_hbonds =
        system->num_hbonds[i] =
          (int)(MAX( (hbonds->end_index[Hindex] - hbonds->index[Hindex]) *saferzone, 
                MIN_HBONDS ));

        if( Hindex < numH-1 ) comp = hbonds->index[Hindex+1];
        else  comp = hbonds->num_intrs;

        if( hbonds->end_index[Hindex] > comp ) 
        {
          fprintf(stderr,"step%d-hbondchk failed: H=%d end(H)=%d str(H+1)=%d\n",
                  step, Hindex, hbonds->end_index[Hindex], comp );
          MPI_Abort( comm, INSUFFICIENT_MEMORY );
        }//if
      }//if
    }//for
  
  }//if-numH
}

void sort_bonds_C(init_forces_pack_t * param)
{
  reax_list **lists   = param->lists;
  reax_list *bonds    = *lists + BONDS;
  reax_system_c *system = param->system;
  storage * workspace = param->workspace;

  int i, j;
  for (i = 0; i < system->N; i ++)
  {
    int start_i = Start_Index( i, bonds );
    int end_i = End_Index( i, bonds);
    int pj;
    for( pj = start_i; pj < end_i; ++pj ) 
    {
      bond_data *jbond = bonds->select.bond_list + pj;
      j = jbond->nbr;

      if(workspace->bond_mark[i] > workspace->bond_mark[j] + 1) 
        workspace->bond_mark[i] = workspace->bond_mark[j] + 1;

      //int start_j = Start_Index( j, bonds);
      //int end_j = End_Index( j, bonds );
      //int pk;//
      //int flag = 0;
      //for ( pk = start_j; pk < end_j; ++pk ) 
      //{
      //  bond_data *kbond = bonds->select.bond_list + pk;
      //  int k = kbond->nbr;
      //  if (k == i)
      //  {
      //    jbond->sym_index = pk;
      //    flag = 1;
      //  }
      //}//for
      //if (!flag)
      //  printf("sym not found %d %d\n", i, j);
    }//for-pj
  }//for-i
}

#endif


#ifdef CPE
#include "STUBS/mpi.h"
#include <slave.h>
#include <dma.h>
#include <math.h>
#include "dma_macros.h"
#include "poly_math.h"
#define ISTEP 64
#define SNSTEP 64
#define MAXBOND 35
//#define LWPF_KERNELS _K(ALL) K(BEFORE) K(GET)  K(PJ) K(BOP)   
//#define LWPF_UNIT U(FORCES)
//#include "lwpf.h"

//#define LWPF_UNIT U(ADD_DBOND)
//#define LWPF_KERNELS K(ALL) K(UPD) K(UPD2) K(FLSH)//  K(CMP) K(CMP_TOR)  
//#define DMA_FAST
//#include "lwpf2.h"


void Init_Forces_noQEq_Full_C_para(init_forces_pack_t * param)
{
  //lwpf_start(ALL);
  //lwpf_start(BEFORE);
  dma_init();
  init_forces_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(init_forces_pack_t));
  dma_syn();
  reax_system_c *system       = l_pm.system;
  control_params *control     = l_pm.control;
  storage *workspace          = l_pm.workspace;
  reax_list **lists           = l_pm.lists;
  //double *Cdbo_list           = l_pm.Cdbo_list;
  //double *Cdbopi_list         = l_pm.Cdbopi_list;
  //double *Cdbopi2_list        = l_pm.Cdbopi2_list;
  bond_order_data *bo_data_list = l_pm.bo_data_list;


  int i, j, pj;
  int start_i, end_i;
  int type_i, type_j;
  int btop_i, bend_i, num_bonds;
  int ihb, jhb, ihb_top, jhb_top;
  int local, flag;
  double cutoff;
  far_neighbor_data_full *nbr_pj;
  reax_atom              *atom_i, *atom_j;

  double bo_cut   = l_pm.bo_cut;
  double bond_cut = l_pm.bond_cut;
  int nlocal      = l_pm.n;
  int natoms      = l_pm.N;
  int num_intrs   = l_pm.num_intrs;
  int ntypes      = l_pm.num_atom_types;

  int ist, isz, ied;
  int ii, ioff;
  int fidx_st[ISTEP], fidx_ed[ISTEP], bidx_st[ISTEP], bidx_ed[ISTEP];
  atom_pack_t packed_atoms[ISTEP];
  far_neighbor_data_full pj_nbr[SNSTEP];
  bond_data       *l_ibond, l_ibond_stor;//
  bond_data       l_blists[MAXBOND];
  double          l_Cdbo[MAXBOND];
  double          *l_iCdbo;
  double          l_Cdbopi[MAXBOND];
  double          *l_iCdbopi;
  double          l_Cdbopi2[MAXBOND];
  double          *l_iCdbopi2;

  double          l_BO[MAXBOND];
  double          *l_iBO;

  rvec2           l_BOpi[MAXBOND];
  double          *l_iBOpi;
  
  bond_order_data l_bo_data[MAXBOND];
  bond_order_data *l_ibo_data;


  int             l_bdmk[ISTEP];
  double          l_norders[ISTEP][2];
  rvec            l_dDeltap[ISTEP];
  //l_ibond = &l_ibond_stor;
  
  init_forces_tbp l_tbp[ntypes][ntypes];
  init_forces_sbp l_sbp[ntypes];
  init_forces_tbp *twbp; 
  init_forces_sbp *sbp_i, *sbp_j;

  pe_get(l_pm.pack_sbp, l_sbp, sizeof(init_forces_sbp)*ntypes);
  pe_get(l_pm.pack_tbp, l_tbp, sizeof(init_forces_tbp)*ntypes*ntypes);
  dma_syn();

  num_bonds = 0;
  btop_i = 0;
  
  intv8 reg_bonds[2];
  for(i = 0; i < 2; i++)
    reg_bonds[i] = 0;
  int *p_nbonds = (int*)(void*)reg_bonds;

  int jst, jed, jsz;
  volatile dma_desc get_desc = 0;
  volatile int get_reply = 0;
  dma_set_mode(&get_desc, PE_MODE);
  dma_set_op(&get_desc, DMA_GET);
  dma_set_reply(&get_desc, &get_reply);

  //lwpf_stop(BEFORE);
  for(ist = _MYID * ISTEP; ist < natoms; ist += ISTEP * 64)
  {
    ied = ist + ISTEP;
    if(ied > natoms)
      ied = natoms;
    isz = ied - ist;
  
    //lwpf_start(GET);
    pe_get(l_pm.packed_atoms+ist,     packed_atoms, sizeof(atom_pack_t)*isz);
    pe_get(l_pm.findex+ist,           fidx_st,      sizeof(int)*isz);
    pe_get(l_pm.fend_index+ist,       fidx_ed,      sizeof(int)*isz);
    pe_get(l_pm.bindex+ist,           bidx_st,      sizeof(int)*isz);
    pe_get(l_pm.bend_index+ist,       bidx_ed,      sizeof(int)*isz);
    pe_get(l_pm.bond_mark+ist,        l_bdmk,       sizeof(int)*isz);
    pe_get(l_pm.dDeltap_self+ist,     l_dDeltap,    sizeof(rvec)*isz);
    //pe_get(l_pm.total_bond_order+ist, l_norders,    sizeof(double)*isz);
    pe_get(l_pm.bo_dboc[ist], l_norders[0],    sizeof(double)*isz*2);
    dma_syn();

    //lwpf_stop(GET);
    for(ii = ist; ii < ied; ++ii )
    {
      i = ii;
      ioff = i - ist;
    
      if(i < nlocal) l_bdmk[ioff] = 0;
      else l_bdmk[ioff] = 1000;
      type_i = packed_atoms[ioff].type;
      if (type_i < 0) continue;
      start_i = fidx_st[ioff];
      end_i   = fidx_ed[ioff];
      btop_i  = bidx_ed[ioff];
      sbp_i   = &(l_sbp[type_i]);
      cutoff  = bond_cut;
      l_ibond = l_blists;
      l_iCdbo     = l_Cdbo;
      l_iCdbopi   = l_Cdbopi;
      l_iCdbopi2  = l_Cdbopi2;
      l_iBO       = l_BO;
      l_iBOpi     = l_BOpi[0];
      l_ibo_data  = l_bo_data;
      for(jst = start_i; jst < end_i; jst += SNSTEP) 
      {
        jsz = SNSTEP;
        if(jst + SNSTEP > end_i)
          jsz = end_i - jst;
      
        pe_get(l_pm.far_nbr_list_full+jst, pj_nbr, sizeof(far_neighbor_data_full)*jsz);
        dma_syn();
        //lwpf_start(PJ);
        for(pj = 0; pj < jsz; ++pj) 
        {
          nbr_pj = &(pj_nbr[pj]);
          j = nbr_pj->nbr;
          //if (j >= natoms || j < 0) continue;
          type_j = nbr_pj->type;
          if(nbr_pj->d <= cutoff) flag = 1;
          else flag = 0;

          if(flag) 
          {
            if (type_j < 0) continue;
            sbp_j = &(l_sbp[type_j]);
            twbp  = &(l_tbp[type_i][type_j]);
            int flag_bopsingle = 0;
          
            //lwpf_start(BOP);
            double r2, C12, C34, C56;
            double Cln_BOp_s, Cln_BOp_pi, Cln_BOp_pi2;
            double BO, BO_s, BO_pi, BO_pi2;
            bond_data *ibond, *pibond;
            bond_order_data *bo_ij;
            r2 = SQR(nbr_pj->d);
            if( sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0 ) 
            {
              C12  = twbp->p_bo1 * p_powd( nbr_pj->d / twbp->r_s, twbp->p_bo2);
              //BO_s = (1.0 + bo_cut) * exp(C12);
              if(C12 < -300)  
                BO_s = 0;
              else if(C12 > 700)
                BO_s = 1.0/0;
              else
                BO_s = (1.0 + bo_cut) * p_expd(C12);
            }
            else BO_s = C12 = 0.0;

            if(sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0) 
            {
              C34   = twbp->p_bo3 * p_powd(nbr_pj->d / twbp->r_p, twbp->p_bo4);
              //BO_pi = exp(C34);
              if(C34 < -300)  
                BO_pi = 0;
              else if(C34 > 700)
                BO_pi = 1.0/0;
              else
                BO_pi = p_expd(C34);
            }
            else BO_pi = C34 = 0.0;
  
            if(sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0) 
            {
              C56    = twbp->p_bo5 * p_powd(nbr_pj->d / twbp->r_pp, twbp->p_bo6);
              //BO_pi2 = exp(C56);
              if(C56 < -300)  
                BO_pi2 = 0;
              else if(C56 > 700)
                BO_pi2 = 1.0/0;
              else
                BO_pi2 = p_expd(C56);
            }
            else BO_pi2 = C56 = 0.0;
          
            BO = BO_s + BO_pi + BO_pi2;
          
            if(BO >= bo_cut) 
            {
              //if (btop_i - bidx_st[ioff] >= MAXBOND)
              //{
              //  call_printf("%d %d %d\n", _MYID, i, btop_i - bidx_st[ioff]);
              //}
              l_ibond = l_blists + ((btop_i - bidx_st[ioff]) % MAXBOND);
              l_iCdbo     = l_Cdbo    + ((btop_i - bidx_st[ioff]) % MAXBOND);
              l_iCdbopi   = l_Cdbopi  + ((btop_i - bidx_st[ioff]) % MAXBOND);
              l_iCdbopi2  = l_Cdbopi2 + ((btop_i - bidx_st[ioff]) % MAXBOND);
              l_iBO       = l_BO      + ((btop_i - bidx_st[ioff]) % MAXBOND);
              l_iBOpi     = &(l_BOpi[(btop_i - bidx_st[ioff]) % MAXBOND][0]);
              l_ibo_data  = l_bo_data + ((btop_i - bidx_st[ioff]) % MAXBOND);

              l_ibond->nbr = nbr_pj->nbr;
              l_ibond->d   = nbr_pj->d;
              rvec_Copy(l_ibond->dvec, nbr_pj->dvec);
              ivec_Copy(l_ibond->rel_box, nbr_pj->rel_box);
              l_ibond->dbond_index = btop_i;

              //bo_ij = &(l_ibond->bo_data);
              //bo_ij->BO     = BO;
              *l_iBO        = BO;

              //bo_ij->BO_s   = BO_s;
              l_ibo_data->BO_s = BO_s;

              //bo_ij->BO_pi  = BO_pi;
              //bo_ij->BO_pi2 = BO_pi2;
              l_iBOpi[0]  = BO_pi;
              l_iBOpi[1]  = BO_pi2;



              Cln_BOp_s   = twbp->p_bo2 * C12 / r2;
              Cln_BOp_pi  = twbp->p_bo4 * C34 / r2;
              Cln_BOp_pi2 = twbp->p_bo6 * C56 / r2;

              //rvec_Scale(bo_ij->dln_BOp_s,  -bo_ij->BO_s  *Cln_BOp_s,  l_ibond->dvec);
              rvec_Scale(l_ibo_data->dln_BOp_s,  -l_ibo_data->BO_s  *Cln_BOp_s,  l_ibond->dvec);

              //rvec_Scale(bo_ij->dln_BOp_pi, -bo_ij->BO_pi *Cln_BOp_pi, l_ibond->dvec);
              //rvec_Scale(bo_ij->dln_BOp_pi2,-bo_ij->BO_pi2*Cln_BOp_pi2,l_ibond->dvec);
              //rvec_Scale(bo_ij->dln_BOp_pi, -l_iBOpi[0]*Cln_BOp_pi, l_ibond->dvec);
              //rvec_Scale(bo_ij->dln_BOp_pi2,-l_iBOpi[1]*Cln_BOp_pi2,l_ibond->dvec);
              rvec_Scale(l_ibo_data->dln_BOp_pi, -l_iBOpi[0]*Cln_BOp_pi, l_ibond->dvec);
              rvec_Scale(l_ibo_data->dln_BOp_pi2,-l_iBOpi[1]*Cln_BOp_pi2,l_ibond->dvec);


              //rvec_Scale(bo_ij->dBOp, -(bo_ij->BO_s  *Cln_BOp_s + bo_ij->BO_pi *Cln_BOp_pi + bo_ij->BO_pi2*Cln_BOp_pi2), l_ibond->dvec);
              //rvec_Scale(bo_ij->dBOp, -(bo_ij->BO_s  *Cln_BOp_s + l_iBOpi[0] *Cln_BOp_pi + l_iBOpi[1]*Cln_BOp_pi2), l_ibond->dvec);
              rvec_Scale(l_ibo_data->dBOp, -(l_ibo_data->BO_s  *Cln_BOp_s + l_iBOpi[0] *Cln_BOp_pi + l_iBOpi[1]*Cln_BOp_pi2), l_ibond->dvec);

              //rvec_Add(l_dDeltap[ioff], bo_ij->dBOp);
              rvec_Add(l_dDeltap[ioff], l_ibo_data->dBOp);
              //bo_ij->BO_s -= bo_cut;
              l_ibo_data->BO_s -= bo_cut;

              //bo_ij->BO   -= bo_cut;
              *l_iBO      -= bo_cut;

              //l_norders[ioff] += bo_ij->BO; 
              //l_norders[ioff][0] += bo_ij->BO; 
              l_norders[ioff][0] += *l_iBO; 
              
              //bo_ij->Cdbo = 0.0;
              *l_iCdbo = 0.0;
              //bo_ij->Cdbopi = 0.0;
              *l_iCdbopi = 0.0;
              
              //bo_ij->Cdbopi2 = 0.0;
              *l_iCdbopi2 = 0.0;

              flag_bopsingle  = 1;
              num_bonds += 1;
              btop_i    += 1;
            }//if
          
            //lwpf_stop(BOP);
          }//if
        }//for-pj
        //lwpf_stop(PJ);
      }//for-jst
      bidx_ed[ioff] = btop_i;

      //sort
      int i_start = bidx_st[ioff];
      int i_stop  = bidx_ed[ioff];
      if (i_start < i_stop)
      {
        my_bubble_sort(l_blists, 
                       l_Cdbo,  l_Cdbopi, l_Cdbopi2, l_BO, l_BOpi, l_bo_data, 
                       i_stop - i_start);
        pe_put(l_pm.bond_list + i_start,    l_blists,   sizeof(bond_data)*(i_stop-i_start));
        pe_put(l_pm.Cdbo_list + i_start,    l_Cdbo,     sizeof(double)*(i_stop-i_start));
        pe_put(l_pm.Cdbopi_list + i_start,  l_Cdbopi,   sizeof(double)*(i_stop-i_start));
        pe_put(l_pm.Cdbopi2_list + i_start, l_Cdbopi2,  sizeof(double)*(i_stop-i_start));
        pe_put(l_pm.BO_list + i_start,      l_BO,       sizeof(double)*(i_stop-i_start));
        pe_put(l_pm.BOpi_list + i_start,    l_BOpi,     sizeof(double)*(i_stop-i_start)*2);
        pe_put(l_pm.bo_data_list + i_start, l_bo_data, sizeof(bond_order_data)*(i_stop-i_start));
        dma_syn();
      }
    }//for-ii
  
    pe_put(l_pm.bend_index+ist,       bidx_ed,   sizeof(int)*isz);
    pe_put(l_pm.bond_mark+ist,        l_bdmk,    sizeof(int)*isz);
    pe_put(l_pm.dDeltap_self+ist,     l_dDeltap, sizeof(rvec)*isz);
    pe_put(l_pm.bo_dboc[ist], l_norders[0], sizeof(double)*isz*2);
    dma_syn();
  }//for-ist
  *p_nbonds = num_bonds;
  reg_reduce_inplace_intv8(reg_bonds, 1);

  if(_MYID == 0)
  {
    pe_put(l_pm.num_bonds, p_nbonds, sizeof(int));
    dma_syn();
  }//if-MYID
  //lwpf_stop(ALL);
}


#define DT_CACHE_FC fCdDelta,rvec4,5,3,6,FADDD
#define KSTEP 32
//read cache;
#define READ_C_H    9
#define READ_C_S    4
#define READ_C_LSZ  (1 << READ_C_S)
#define READ_C_LCNT (1 << (READ_C_H - READ_C_S))
#define READ_C_MM   (READ_C_LSZ - 1)
#define READ_C_LM   (READ_C_LCNT - 1)

void read_cache_cdbo(int i, double Cdbo_cache[][READ_C_LSZ], 
                    double Cdbopi_cache[][READ_C_LSZ],
                    double Cdbopi2_cache[][READ_C_LSZ],
                    int *read_ctag,
                    double *Cdbo_list, 
                    double *Cdbopi_list, 
                    double *Cdbopi2_list,
                    double *cdbo, double *cdbopi, double *cdbopi2)
{
  dma_init();
  if (read_ctag[(i >> READ_C_S) & READ_C_LM] != i >> READ_C_S)
  {
    pe_get(Cdbo_list + (i & ~READ_C_MM),    Cdbo_cache[(i >> READ_C_S) & READ_C_LM],    sizeof(double) * READ_C_LSZ);
    pe_get(Cdbopi_list + (i & ~READ_C_MM),  Cdbopi_cache[(i >> READ_C_S) & READ_C_LM],  sizeof(double) * READ_C_LSZ);
    pe_get(Cdbopi2_list + (i & ~READ_C_MM), Cdbopi2_cache[(i >> READ_C_S) & READ_C_LM], sizeof(double) * READ_C_LSZ);
    dma_syn();
    read_ctag[(i >> READ_C_S) & READ_C_LM] = i >> READ_C_S;
  }
  *cdbo     = Cdbo_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
  *cdbopi   = Cdbopi_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
  *cdbopi2  = Cdbopi2_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
}

void read_cache_fCdDelta(int i, rvec4 fCdDelta_cache[][READ_C_LSZ], 
                    int *fc_tag, rvec4 *fCdDelta, double *fc)
{
  dma_init();
  if (fc_tag[(i >> READ_C_S) & READ_C_LM] != i >> READ_C_S)
  {
    pe_get(fCdDelta + (i & ~READ_C_MM), 
          fCdDelta_cache[(i >> READ_C_S) & READ_C_LM], 
          sizeof(rvec4) * READ_C_LSZ);
    dma_syn();
    fc_tag[(i >> READ_C_S) & READ_C_LM] = i >> READ_C_S;
  }
  *fc = fCdDelta_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM][3];
}

/*****no if branch*****/
void Add_All_dBond_to_Forces_C_Para(add_dbond_pack_t *param) 
{
  //lwpf_enter(ADD_DBOND);
  //lwpf_start(ALL);

  dma_init();
  add_dbond_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(add_dbond_pack_t));
  dma_syn();

  reax_system_c *system           = l_pm.system;
  control_params *control         = l_pm.control;
  storage *workspace              = l_pm.workspace;
  reax_list **lists               = l_pm.lists;
  double *Cdbo_list               = l_pm.Cdbo_list;
  double *Cdbopi_list             = l_pm.Cdbopi_list;
  double *Cdbopi2_list            = l_pm.Cdbopi2_list;
  bond_order_data *bo_data_list   = l_pm.bo_data_list;
  rvec4 *fCdDelta                 = l_pm.fCdDelta;
  int i, pj, pj_off;
  //reax_list *bonds = (*lists) + BONDS;
  //double *Cdbo_list = bonds->Cdbo_list;
  //double *Cdbopi_list = bonds->Cdbopi_list;
  //double *Cdbopi2_list = bonds->Cdbopi2_list;
  int start_i, end_i, len_i;
  bond_data *nbr_j, *nbr_k;
  int pk, k, j, pk_off;
  //rvec temp;
  rvec4 temp;
  temp[0] = temp[1] = temp[2] = temp[3] = 0.0;
  double buf_c;
  int buf_idx[MAX_BONDS];
  rvec buf_dBOp[MAX_BONDS];
  bond_order_data *bo_t;
  
  int index[KSTEP], end_index[KSTEP];
  rvec dDeltap_self_i[KSTEP];
  rvec4 fCdDelta_i[KSTEP];

  bond_data bond_list_j[MAX_LEN];
  bond_order_data bo_data_j[MAX_LEN];
  double Cdbo_list_j[MAX_LEN];
  double Cdbopi_list_j[MAX_LEN];
  double Cdbopi2_list_j[MAX_LEN];

  //read_cache;
  int read_ctag[READ_C_LCNT];
  double Cdbo_cache[READ_C_LCNT][READ_C_LSZ];
  double Cdbopi_cache[READ_C_LCNT][READ_C_LSZ];
  double Cdbopi2_cache[READ_C_LCNT][READ_C_LSZ];
  double cdbo, cdbopi, cdbopi2;

  int fc_tag[READ_C_LCNT];
  rvec4 fCdDelta_cache[READ_C_LCNT][READ_C_LSZ];
  double fc;
  for(i = 0; i < READ_C_LCNT; i++)
  {
    read_ctag[i]  = -1;
    fc_tag[i]     = -1;
  }


  //swcache_lock_t *locks_dbond = l_pm.locks_dbond;
  //SWCACHE_INIT_U(DT_CACHE_FC, locks_dbond);
  swcache_lock_t *locks_frc = l_pm.locks_frc;
  SWCACHE_INIT_U(DT_CACHE_FC, locks_frc);

  int ist, isz, ied, ioff;
  int N = l_pm.N;

  for(ist = _MYID * KSTEP; ist < N; ist+=KSTEP * 64)
  {
    ied = ist + KSTEP;
    if(ied > N)
      ied = N;
    isz = ied - ist;

    pe_get(l_pm.index+ist,        index,          sizeof(int)*isz);
    pe_get(l_pm.end_index+ist,    end_index,      sizeof(int)*isz);
    pe_get(l_pm.dDeltap_self+ist, dDeltap_self_i, sizeof(rvec)*isz);
    pe_get(l_pm.fCdDelta+ist,     fCdDelta_i,     sizeof(rvec4)*isz);
    dma_syn();

    //for( i = 0; i < system->N; ++i )
    for(i = ist; i < ied; i++)
    {
      ioff = i - ist;
      //start_i = Start_Index(i, bonds);
      //end_i   = End_Index(i, bonds);
      //start_i = l_pm.index[i];//Start_Index(i, bonds);
      //end_i   = l_pm.end_index[i];//End_Index(i, bonds);
      
      start_i = index[ioff];//Start_Index(i, bonds);
      end_i   = end_index[ioff];//End_Index(i, bonds);
      len_i   = end_i - start_i;
      buf_c = 0.0;
      if(len_i <= 0) continue;
      pe_get(l_pm.Cdbo_list+start_i,    Cdbo_list_j,    sizeof(double)*len_i);
      pe_get(l_pm.Cdbopi_list+start_i,  Cdbopi_list_j,  sizeof(double)*len_i);
      pe_get(l_pm.Cdbopi2_list+start_i, Cdbopi2_list_j, sizeof(double)*len_i);
      pe_get(l_pm.bond_list+start_i,    bond_list_j,    sizeof(bond_data)*len_i);
      pe_get(l_pm.bo_data_list+start_i, bo_data_j,      sizeof(bond_order_data)*len_i);
      dma_syn();

      for( pj = start_i; pj < end_i; ++pj )
      {
        pj_off = pj - start_i;
        //nbr_j = &(bonds->select.bond_list[pj]);
        nbr_j = &(bond_list_j[pj_off]);
        j = nbr_j->nbr;

        bond_order_data *bo_ij, *bo_ji;
        dbond_coefficients coef;

        /* Virial Tallying variables */
        rvec fi_tmp, fj_tmp, fk_tmp, delij, delji, delki, delkj;

        /* Initializations */
        //bo_ij = &(bonds->bo_data_list[pj]);
        bo_ij = &(bo_data_j[pj_off]);
        //bo_ji = &(bonds->bo_data_list[ nbr_j->sym_index]);

        read_cache_cdbo(nbr_j->sym_index, Cdbo_cache, Cdbopi_cache, Cdbopi2_cache,
                        read_ctag, Cdbo_list, Cdbopi_list, Cdbopi2_list,
                        &cdbo, &cdbopi, &cdbopi2);
        read_cache_fCdDelta(j, fCdDelta_cache, fc_tag, fCdDelta, &fc);

        coef.C1dbo = bo_ij->C1dbo * (Cdbo_list_j[pj_off] + cdbo);//Cdbo_list[nbr_j->sym_index]);
        coef.C2dbo = bo_ij->C2dbo * (Cdbo_list_j[pj_off] + cdbo);//Cdbo_list[nbr_j->sym_index]);
        coef.C3dbo = bo_ij->C3dbo * (Cdbo_list_j[pj_off] + cdbo);//Cdbo_list[nbr_j->sym_index]);

        coef.C1dbopi = bo_ij->C1dbopi * (Cdbopi_list_j[pj_off] + cdbopi);//Cdbopi_list[nbr_j->sym_index]);
        coef.C2dbopi = bo_ij->C2dbopi * (Cdbopi_list_j[pj_off] + cdbopi);//Cdbopi_list[nbr_j->sym_index]);
        coef.C3dbopi = bo_ij->C3dbopi * (Cdbopi_list_j[pj_off] + cdbopi);//Cdbopi_list[nbr_j->sym_index]);
        coef.C4dbopi = bo_ij->C4dbopi * (Cdbopi_list_j[pj_off] + cdbopi);//Cdbopi_list[nbr_j->sym_index]);

        coef.C1dbopi2 = bo_ij->C1dbopi2 * (Cdbopi2_list_j[pj_off] + cdbopi2);//Cdbopi2_list[nbr_j->sym_index]);
        coef.C2dbopi2 = bo_ij->C2dbopi2 * (Cdbopi2_list_j[pj_off] + cdbopi2);//Cdbopi2_list[nbr_j->sym_index]);
        coef.C3dbopi2 = bo_ij->C3dbopi2 * (Cdbopi2_list_j[pj_off] + cdbopi2);//Cdbopi2_list[nbr_j->sym_index]);
        coef.C4dbopi2 = bo_ij->C4dbopi2 * (Cdbopi2_list_j[pj_off] + cdbopi2);//Cdbopi2_list[nbr_j->sym_index]);

        //coef.C1dDelta = bo_ij->C1dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);
        //coef.C2dDelta = bo_ij->C2dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);
        //coef.C3dDelta = bo_ij->C3dbo * (workspace->fCdDelta[i][3]+workspace->fCdDelta[j][3]);


        coef.C1dDelta = bo_ij->C1dbo * (fCdDelta_i[ioff][3]+fc);//workspace->fCdDelta[j][3]);
        coef.C2dDelta = bo_ij->C2dbo * (fCdDelta_i[ioff][3]+fc);//workspace->fCdDelta[j][3]);
        coef.C3dDelta = bo_ij->C3dbo * (fCdDelta_i[ioff][3]+fc);//workspace->fCdDelta[j][3]);

        // forces on i
        rvec4_Scale(           temp, coef.C1dbo,    bo_ij->dBOp );
        //rvec_ScaledAdd( temp, coef.C2dbo,    workspace->dDeltap_self[i] );
        rvec4_ScaledAdd( temp, coef.C2dbo,    dDeltap_self_i[ioff]);
        rvec4_ScaledAdd( temp, coef.C1dDelta, bo_ij->dBOp );
      
        //rvec_ScaledAdd( temp, coef.C2dDelta, workspace->dDeltap_self[i] );
        rvec4_ScaledAdd( temp, coef.C2dDelta, dDeltap_self_i[ioff]);
        rvec4_ScaledAdd( temp, coef.C1dbopi,  bo_ij->dln_BOp_pi );
        rvec4_ScaledAdd( temp, coef.C2dbopi,  bo_ij->dBOp );
        //rvec_ScaledAdd( temp, coef.C3dbopi,  workspace->dDeltap_self[i]);
        rvec4_ScaledAdd( temp, coef.C3dbopi,  dDeltap_self_i[ioff]);
        rvec4_ScaledAdd( temp, coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
        rvec4_ScaledAdd( temp, coef.C2dbopi2, bo_ij->dBOp );
        //rvec_ScaledAdd( temp, coef.C3dbopi2, workspace->dDeltap_self[i] );
        rvec4_ScaledAdd( temp, coef.C3dbopi2, dDeltap_self_i[ioff]);
       
        //lwpf_start(UPD);
        //rvec4_Add( workspace->fCdDelta[i], temp );
        SWCACHE_UPDATE(DT_CACHE_FC, i, (&temp));
        //lwpf_stop(UPD);
        
        buf_c += (-coef.C2dbo - coef.C2dDelta - coef.C3dbopi - coef.C3dbopi2);
        //rvec4_Copy(buf_dBOp[pj_off], bonds->bo_data_list[pj].dBOp);
        rvec4_Copy(buf_dBOp[pj_off], bo_data_j[pj_off].dBOp);
        buf_idx[pj_off] = j;
      }//for-pj
     
      // forces on k: i neighbor
      for( pj = start_i; pj < end_i; ++pj) 
      {
        pj_off = pj - start_i;
        rvec4_Scale(     temp, buf_c, buf_dBOp[pj_off]);
        //lwpf_start(UPD2);
        //rvec4_Add( workspace->fCdDelta[buf_idx[pj_off]], temp);
        SWCACHE_UPDATE(DT_CACHE_FC, buf_idx[pj_off], (&temp));
        //lwpf_stop(UPD2);
      }

    }//for-i
    //pe_put(l_pm.dDeltap_self+ist, dDeltap_self_j, sizeof(rvec)*isz);
    //dma_syn();

  }//ist
  //lwpf_start(FLSH);
  SWCACHE_FLUSH(DT_CACHE_FC);
  //lwpf_stop(FLSH);

  //lwpf_stop(ALL);
  //lwpf_exit(ADD_DBOND);
}

void Init_Forces_noQEq_HB_Full_C_para(init_forces_pack_t * param) 
{
  dma_init();
  init_forces_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(init_forces_pack_t));
  dma_syn();

  reax_system_c *system       = l_pm.system;
  control_params *control     = l_pm.control;
  simulation_data *data       = l_pm.data;
  storage *workspace          = l_pm.workspace;
  reax_list **lists           = l_pm.lists;
  int i, j, pj;
  int start_i, end_i;
  int type_i, type_j;
  int num_hbonds;
  int ihb, jhb, ihb_top;
  int local, flag, renbr;
  double cutoff;
  two_body_parameters *twbp;
  far_neighbor_data_full *nbr_pj;
  num_hbonds = 0;
  if (l_pm.hbond_cut <= 0) return;

  int ist, isz, ied, ioff;
  intv8 reg_hbonds[2];
  for(i = 0; i < 2; i++)
    reg_hbonds[i] = 0;
  int *p_nhbonds = (int*)(void*)reg_hbonds;
  int n = l_pm.n;
  int ntypes  = l_pm.num_atom_types;
  cutoff      = l_pm.hbond_cut;
  
  int fidx_st[ISTEP], fidx_ed[ISTEP], hidx_ed[ISTEP];
  int sys_hidx[ISTEP];
  atom_pack_t packed_atoms[ISTEP];
  far_neighbor_data_full pj_nbr[SNSTEP];
  hbond_data hblist[SNSTEP];

  init_forces_sbp l_sbp[ntypes];
  init_forces_sbp *sbp_i, *sbp_j;
  pe_get(l_pm.pack_sbp, l_sbp, sizeof(init_forces_sbp)*ntypes);
  dma_syn();

  volatile dma_desc get_desc = 0;
  volatile int get_reply = 0;
  dma_set_mode(&get_desc, PE_MODE);
  dma_set_op(&get_desc, DMA_GET);
  dma_set_reply(&get_desc, &get_reply);

  int jst, jed, jsz;
  int ihb_old, ihb_cnt;
  for(ist = _MYID*ISTEP; ist < n; ist+=ISTEP*64) 
  {
    ied = ist + ISTEP;
    if(ied > n)
      ied = n;
    isz = ied - ist;
    pe_get(l_pm.packed_atoms+ist, packed_atoms, sizeof(atom_pack_t)*isz);
    pe_get(l_pm.findex+ist,       fidx_st,      sizeof(int)*isz);
    pe_get(l_pm.fend_index+ist,   fidx_ed,      sizeof(int)*isz);
    pe_get(system->Hindex+ist,    sys_hidx,     sizeof(int)*isz);
    dma_syn();

    for(i = ist; i < ied; ++i) 
    {
      ioff = i - ist;
      type_i = packed_atoms[ioff].type;
      if (type_i < 0) continue;
      start_i = fidx_st[ioff];
      end_i   = fidx_ed[ioff];
      sbp_i   = &(l_sbp[type_i]);
      ihb = -1;
      ihb_top = -1;
      ihb = sbp_i->p_hbond;

      if(ihb == 1)
      {
        pe_get(l_pm.hend_index+sys_hidx[ioff], &ihb_top, sizeof(int));
        dma_syn();
        //ihb_top = l_pm.hend_index[sys_hidx[ioff]];
        ihb_old = ihb_top;
        ihb_cnt = 0;
        for(jst = start_i; jst < end_i; jst += SNSTEP) 
        {
          jsz = SNSTEP;
          if(jst + SNSTEP > end_i)
            jsz = end_i - jst;
          //if(jsz > 0)
          //{
          //  //pe_get(l_pm.far_nbr_list_full+jst, pj_nbr, 
          //  //        sizeof(far_neighbor_data_full)*jsz);
          //  //dma_syn();
          //  dma_set_size(&get_desc, sizeof(far_neighbor_data_full)*jsz);
          //  get_reply = 0;
          //  dma_rpl(get_desc, l_pm.far_nbr_list_full+jst, pj_nbr, get_reply);
          //  while(get_reply != 1);

          //}

          for(pj = 0; pj < jsz; ++pj) 
          {
            nbr_pj = &(l_pm.far_nbr_list_full[jst+pj]);
            //nbr_pj = &(pj_nbr[pj]);
            if(nbr_pj->d <= cutoff)
            {
              type_j = nbr_pj->type;
              j = nbr_pj->nbr;
              if (type_j < 0) continue;
              sbp_j = &(l_sbp[type_j]);
              if((ihb==1 || ihb==2)) 
              {
                jhb = sbp_j->p_hbond;
                if(ihb == 1 && jhb == 2) 
                {
                  hblist[ihb_cnt].nbr = j;
                  hblist[ihb_cnt].scl = 1;
                  hblist[ihb_cnt].ptr = nbr_pj;
                  ++ihb_cnt;
                }//if
              }//if
            }//if-cutoff
          }//for-pj 
        }//for-jst
        if(ihb_cnt > 0)
        {
          ihb_top += ihb_cnt;
          num_hbonds += ihb_cnt;
          pe_put(l_pm.hbond_list+ihb_old, hblist, 
            sizeof(hbond_data)*(ihb_cnt));
          dma_syn();
        }
        //l_pm.hend_index[sys_hidx[ioff]] = ihb_top;
        pe_put(l_pm.hend_index+sys_hidx[ioff], &ihb_top, sizeof(int));
        dma_syn();
      }//if-ihb
    }//for-i
  }//for-ist
  *p_nhbonds = num_hbonds;
  reg_reduce_inplace_intv8(reg_hbonds, 1);

  if(_MYID == 0)
  {
    pe_put(l_pm.num_hbonds, p_nhbonds, sizeof(int));
    dma_syn()
  }
}
#endif
