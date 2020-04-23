#include <math.h>
#include "reaxc_defs_sunway.h"
#include "reaxc_ctypes_sunway.h"
#include "reaxc_inlines_sw64.h"
#include "reaxc_bond_orders_sw64.h"
#include "stdio.h"
#include "stdlib.h"
#include "gptl.h"
#include "sunway.h"
#include "simd.h"
//#include "sleef_math.h"
#define ISTEP 32
#define MAX_LEN 20

#ifdef MPE
#include <mpi.h>
#include <athread.h>
//#define LWPF_UNITS U(BOFUNC)
//#include "lwpf2.h"

extern SLAVE_FUN(BO_Section1_C_Para)(bo_param_pack_t *);
extern SLAVE_FUN(BO_Section2_C_Para)(bo_param_pack_t *);
extern SLAVE_FUN(BO_Section3_C_Para)(bo_param_pack_t *);

void BO_Section1_C(bo_param_pack_t *param)
{
  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;

  int i, j, pj, type_i, type_j;
  int start_i, end_i, sym_index;
  double val_i, Deltap_i, Deltap_boc_i;
  double val_j, Deltap_j, Deltap_boc_j;
  double f1, f2, f3, f4, f5, f4f5, exp_f4, exp_f5;
  double exp_p1i,        exp_p2i, exp_p1j, exp_p2j;
  double temp, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
  double Cf45_ij, Cf45_ji, p_lp1; //u_ij, u_ji
  double A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji;
  double explp1, p_boc1, p_boc2;
  single_body_parameters *sbp_i, *sbp_j;
 
  #define SW_BO_SEC1
  #ifdef SW_BO_SEC1

  if(athread_idle() == 0)
    athread_init();
  athread_spawn(BO_Section1_C_Para, param);
  athread_join();

  #else

  for( i = 0; i < system->N; ++i ) 
  {
    type_i = system->packed_atoms[i].type;
    if (type_i < 0) continue;
    sbp_i = &(system->reax_param.sbp[type_i]);
    workspace->Deltap[i]           = workspace->bo_dboc[i][0] - sbp_i->valency;
    workspace->Deltap_boc[i]       = workspace->bo_dboc[i][0] - sbp_i->valency_val;
    workspace->bo_dboc[i][0] = 0;
  }
  #endif
}
void BO_Section2_C(bo_param_pack_t *param)
{
  //printf("Section2\n");

  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;

  int i, j, pj, type_i, type_j;
  int start_i, end_i, sym_index;
  double val_i, Deltap_i, Deltap_boc_i;
  double val_j, Deltap_j, Deltap_boc_j;
  double f1, f2, f3, f4, f5, f4f5, exp_f4, exp_f5;
  double exp_p1i,        exp_p2i, exp_p1j, exp_p2j;
  double temp, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
  double Cf45_ij, Cf45_ji, p_lp1; //u_ij, u_ji
  double A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji;
  double explp1, p_boc1, p_boc2;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  bond_order_data *bo_ij, *bo_ji;
  reax_list *bonds = (*lists) + BONDS;
  double *BO_list     = bonds->BO_list;
  rvec2 *BOpi_list    = bonds->BOpi_list;
  bond_order_data *bo_data_list = bonds->bo_data_list;

  p_boc1 = system->reax_param.gp.l[0];
  p_boc2 = system->reax_param.gp.l[1];

  #define SW_BO_SEC2
  #ifdef SW_BO_SEC2
 
  //bo_param_pack_t pm;
  //pm.system       = system;
  //pm.control      = control;
  //pm.data         = data;
  //pm.workspace    = workspace;
  //pm.lists        = lists;
  //pm.out_control  = out_control;
  //pm.packed_atoms = system->packed_atoms;
  //pm.bond_list    = bonds->select.bond_list;
  //pm.BO_list      = bonds->BO_list;
  //pm.BOpi_list    = bonds->BOpi_list;
  //pm.bo_data_list = bonds->bo_data_list;
  //pm.bond_list    = bonds->select.bond_list;
  //pm.index        = bonds->index;
  //pm.end_index    = bonds->end_index;
  //pm.bo_dboc      = workspace->bo_dboc;
  //pm.Deltap       = workspace->Deltap;
  //pm.Deltap_boc   = workspace->Deltap_boc;
  //pm.N            = system->N;
  //pm.ntypes       = system->reax_param.num_atom_types;
  //pm.p_boc1       = p_boc1;
  //pm.p_boc2       = p_boc2;


  //int s1, s2;
  //for(s1 = 0; s1 < pm.ntypes; s1++)
  //{
  //  pm.sbp[s1] = system->reax_param.sbp[s1].valency;
  //  for(s2 = 0; s2 < pm.ntypes; s2++)
  //  {
  //    pm.tbp[s1 * pm.ntypes + s2].p_boc3  = system->reax_param.tbp[s1][s2].p_boc3;
  //    pm.tbp[s1 * pm.ntypes + s2].p_boc4  = system->reax_param.tbp[s1][s2].p_boc4;
  //    pm.tbp[s1 * pm.ntypes + s2].p_boc5  = system->reax_param.tbp[s1][s2].p_boc5;
  //    pm.tbp[s1 * pm.ntypes + s2].v13cor  = system->reax_param.tbp[s1][s2].v13cor;
  //    pm.tbp[s1 * pm.ntypes + s2].ovc     = system->reax_param.tbp[s1][s2].ovc;
  //  }
  //}


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

  if(athread_idle() == 0)
    athread_init();
  athread_spawn(BO_Section2_C_Para, param);
  athread_join();
 
  //if(myrank == 0)
  //{
  //  lwpf_report_summary(stdout, &conf);
  //}

  #else
  for( i = 0; i < system->N; ++i ) 
  {
    type_i = system->packed_atoms[i].type;
    if (type_i < 0) continue;
    sbp_i = &(system->reax_param.sbp[type_i]);
    val_i = sbp_i->valency;
    Deltap_i     = workspace->Deltap[i];
    Deltap_boc_i = workspace->Deltap_boc[i];
    start_i = Start_Index(i, bonds);
    end_i   = End_Index(i, bonds);

    for( pj = start_i; pj < end_i; ++pj) 
    {
      j = bonds->select.bond_list[pj].nbr;
      type_j = system->packed_atoms[j].type;
      if (type_j < 0) continue;
      bo_ij = &(bo_data_list[pj]);
      // fprintf( stderr, "\tj:%d - ubo: %8.3f\n", j+1, bo_ij->BO );
     
      /*compute for sym_index*/
      int start_j = Start_Index( j, bonds);
      int end_j = End_Index( j, bonds );
      int pk;//
      for ( pk = start_j; pk < end_j; ++pk ) 
      {
        bond_data *kbond = bonds->select.bond_list + pk;
        int k = kbond->nbr;
        if (k == i)
        {
          bonds->select.bond_list[pj].sym_index = pk;
        }
      }//for
      /*end of sym_index*/

      //if( i < j || workspace->bond_mark[j] > 3 ) 
      {
        twbp = &( system->reax_param.tbp[type_i][type_j] );

        if( twbp->ovc < 0.001 && twbp->v13cor < 0.001 ) 
        {
          bo_ij->C1dbo = 1.000000;
          bo_ij->C2dbo = 0.000000;
          bo_ij->C3dbo = 0.000000;

          //bo_ij->C1dbopi = bo_ij->BO_pi;
          bo_ij->C1dbopi = BOpi_list[pj][0];
          bo_ij->C2dbopi = 0.000000;
          bo_ij->C3dbopi = 0.000000;
          bo_ij->C4dbopi = 0.000000;

          //bo_ij->C1dbopi2 = bo_ij->BO_pi2;
          bo_ij->C1dbopi2 = BOpi_list[pj][1];
          bo_ij->C2dbopi2 = 0.000000;
          bo_ij->C3dbopi2 = 0.000000;
          bo_ij->C4dbopi2 = 0.000000;

        }
        else 
        {
          val_j = system->reax_param.sbp[type_j].valency;
          Deltap_j = workspace->Deltap[j];
          Deltap_boc_j = workspace->Deltap_boc[j];

          /* on page 1 */
          if( twbp->ovc >= 0.001 ) 
          {
            /* Correction for overcoordination */
            exp_p1i = exp( -p_boc1 * Deltap_i );
            exp_p2i = exp( -p_boc2 * Deltap_i );
            exp_p1j = exp( -p_boc1 * Deltap_j );
            exp_p2j = exp( -p_boc2 * Deltap_j );

            f2 = exp_p1i + exp_p1j;
            f3 = -1.0 / p_boc2 * log( 0.5 * ( exp_p2i  + exp_p2j ) );
            f1 = 0.5 * ( ( val_i + f2 )/( val_i + f2 + f3 ) +
                         ( val_j + f2 )/( val_j + f2 + f3 ) );

            temp = f2 + f3;
            u1_ij = val_i + temp;
            u1_ji = val_j + temp;
            Cf1A_ij = 0.5 * f3 * (1.0 / SQR( u1_ij ) +
                                  1.0 / SQR( u1_ji ));
            Cf1B_ij = -0.5 * (( u1_ij - f3 ) / SQR( u1_ij ) +
                              ( u1_ji - f3 ) / SQR( u1_ji ));

            Cf1_ij = 0.50 * ( -p_boc1 * exp_p1i / u1_ij -
                              ((val_i+f2) / SQR(u1_ij)) *
                              ( -p_boc1 * exp_p1i +
                                exp_p2i / ( exp_p2i + exp_p2j ) ) +
                              -p_boc1 * exp_p1i / u1_ji -
                              ((val_j+f2) / SQR(u1_ji)) *
                              ( -p_boc1 * exp_p1i +
                                exp_p2i / ( exp_p2i + exp_p2j ) ));


            Cf1_ji = -Cf1A_ij * p_boc1 * exp_p1j +
              Cf1B_ij * exp_p2j / ( exp_p2i + exp_p2j );

          }
          else 
          {
            /* No overcoordination correction! */
            f1 = 1.0;
            Cf1_ij = Cf1_ji = 0.0;
          }

          if( twbp->v13cor >= 0.001 ) 
          {
            /* Correction for 1-3 bond orders */
            //exp_f4 =exp(-(twbp->p_boc4 * SQR( bo_ij->BO ) -
            exp_f4 =exp(-(twbp->p_boc4 * SQR(BO_list[pj]) -
                          Deltap_boc_i) * twbp->p_boc3 + twbp->p_boc5);
            //exp_f5 =exp(-(twbp->p_boc4 * SQR( bo_ij->BO ) -
            exp_f5 =exp(-(twbp->p_boc4 * SQR( BO_list[pj]) -
                          Deltap_boc_j) * twbp->p_boc3 + twbp->p_boc5);

            f4 = 1. / (1. + exp_f4);
            f5 = 1. / (1. + exp_f5);
            f4f5 = f4 * f5;

            /* Bond Order pages 8-9, derivative of f4 and f5 */
            Cf45_ij = -f4 * exp_f4;
            Cf45_ji = -f5 * exp_f5;
          }
          else 
          {
            f4 = f5 = f4f5 = 1.0;
            Cf45_ij = Cf45_ji = 0.0;
          }

          /* Bond Order page 10, derivative of total bond order */
          A0_ij = f1 * f4f5;
          //A1_ij = -2 * twbp->p_boc3 * twbp->p_boc4 * bo_ij->BO * (Cf45_ij + Cf45_ji);
          A1_ij = -2 * twbp->p_boc3 * twbp->p_boc4 * BO_list[pj] * (Cf45_ij + Cf45_ji);
          A2_ij = Cf1_ij / f1 + twbp->p_boc3 * Cf45_ij;
          A2_ji = Cf1_ji / f1 + twbp->p_boc3 * Cf45_ji;
          A3_ij = A2_ij + Cf1_ij / f1;
          A3_ji = A2_ji + Cf1_ji / f1;

          /* find corrected bond orders and their derivative coef */
          //bo_ij->BO    = bo_ij->BO    * A0_ij;
          BO_list[pj]    = BO_list[pj]  * A0_ij;
          //bo_ij->BO_pi = bo_ij->BO_pi * A0_ij *f1;
          //bo_ij->BO_pi2= bo_ij->BO_pi2* A0_ij *f1;
          BOpi_list[pj][0]= BOpi_list[pj][0]* A0_ij *f1;
          BOpi_list[pj][1]= BOpi_list[pj][1]* A0_ij *f1;

          //bo_ij->BO_s  = bo_ij->BO - ( bo_ij->BO_pi + bo_ij->BO_pi2 );
          //bo_ij->BO_s  = BO_list[pj] - ( bo_ij->BO_pi + bo_ij->BO_pi2 );
          bo_ij->BO_s  = BO_list[pj] - (BOpi_list[pj][0] + BOpi_list[pj][1]);

          //bo_ij->C1dbo = A0_ij + bo_ij->BO * A1_ij;
          //bo_ij->C2dbo = bo_ij->BO * A2_ij;
          //bo_ij->C3dbo = bo_ij->BO * A2_ji;
          bo_ij->C1dbo = A0_ij + BO_list[pj] * A1_ij;
          bo_ij->C2dbo = BO_list[pj] * A2_ij;
          bo_ij->C3dbo = BO_list[pj] * A2_ji;


          bo_ij->C1dbopi = f1*f1*f4*f5;
          //bo_ij->C2dbopi = bo_ij->BO_pi * A1_ij;
          //bo_ij->C3dbopi = bo_ij->BO_pi * A3_ij;
          //bo_ij->C4dbopi = bo_ij->BO_pi * A3_ji;
          bo_ij->C2dbopi = BOpi_list[pj][0] * A1_ij;
          bo_ij->C3dbopi = BOpi_list[pj][0] * A3_ij;
          bo_ij->C4dbopi = BOpi_list[pj][0] * A3_ji;


          bo_ij->C1dbopi2 = f1*f1*f4*f5;
          //bo_ij->C2dbopi2 = bo_ij->BO_pi2 * A1_ij;
          //bo_ij->C3dbopi2 = bo_ij->BO_pi2 * A3_ij;
          //bo_ij->C4dbopi2 = bo_ij->BO_pi2 * A3_ji;
          bo_ij->C2dbopi2 = BOpi_list[pj][1] * A1_ij;
          bo_ij->C3dbopi2 = BOpi_list[pj][1] * A3_ij;
          bo_ij->C4dbopi2 = BOpi_list[pj][1] * A3_ji;

        }

        /* neglect bonds that are < 1e-10 */
        //if( bo_ij->BO < 1e-10 ) bo_ij->BO = 0.0;
        if(BO_list[pj] < 1e-10 ) BO_list[pj] = 0.0;

        if( bo_ij->BO_s < 1e-10 )
          bo_ij->BO_s = 0.0;
        
        //if( bo_ij->BO_pi < 1e-10 )
        //  bo_ij->BO_pi = 0.0;
        //if( bo_ij->BO_pi2 < 1e-10 )
        //  bo_ij->BO_pi2 = 0.0;
        if(BOpi_list[pj][0] < 1e-10 )
          BOpi_list[pj][0] = 0.0;
        if(BOpi_list[pj][1] < 1e-10 )
          BOpi_list[pj][1] = 0.0;

        //workspace->total_bond_order[i] += bo_ij->BO; //now keeps total_BO
        //workspace->bo_dboc[i][0] += bo_ij->BO; //now keeps total_BO
        workspace->bo_dboc[i][0] += BO_list[pj]; //now keeps total_BO

      }
      //else 
      //{
      //  /* We only need to update bond orders from bo_ji
      //     everything else is set in uncorrected_bo calculations */
      //  sym_index     = bonds->select.bond_list[pj].sym_index;
      //  bo_ji         = &(bonds->select.bond_list[ sym_index ].bo_data);
      //  bo_ij->BO     = bo_ji->BO;
      //  bo_ij->BO_s   = bo_ji->BO_s;
      //  bo_ij->BO_pi  = bo_ji->BO_pi;
      //  bo_ij->BO_pi2 = bo_ji->BO_pi2;

      //  workspace->total_bond_order[i] += bo_ij->BO;// now keeps total_BO
      //}
    }
  }
  #endif
}

void BO_Section3_C(bo_param_pack_t *param)
{
  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;


  int i, j, pj, type_i, type_j;
  int start_i, end_i, sym_index;
  double val_i, Deltap_i, Deltap_boc_i;
  double val_j, Deltap_j, Deltap_boc_j;
  double f1, f2, f3, f4, f5, f4f5, exp_f4, exp_f5;
  double exp_p1i,        exp_p2i, exp_p1j, exp_p2j;
  double temp, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
  double Cf45_ij, Cf45_ji, p_lp1; //u_ij, u_ji
  double A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji;
  double explp1, p_boc1, p_boc2;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  bond_order_data *bo_ij, *bo_ji;
  reax_list *bonds = (*lists) + BONDS;
  double *BO_list     = bonds->BO_list;
  rvec2 *BOpi_list    = bonds->BOpi_list;
  bond_order_data *bo_data_list = bonds->bo_data_list;
  
  p_lp1 = system->reax_param.gp.l[15];

  #define SW_BO_SEC3
  #ifdef SW_BO_SEC3
  //bo_param_pack_t pm;
  //pm.system       = system;
  //pm.control      = control;
  //pm.data         = data;
  //pm.workspace    = workspace;
  //pm.lists        = lists;
  //pm.out_control  = out_control;
  //pm.packed_atoms = system->packed_atoms;
  //pm.ntypes       = system->reax_param.num_atom_types;
  //pm.bo_dboc      = workspace->bo_dboc;
  //pm.Delta        = workspace->Delta;
  //pm.Delta_e      = workspace->Delta_e;
  //pm.Delta_val      = workspace->Delta_val;
  //pm.vlpex        = workspace->vlpex;
  //pm.nlp          = workspace->nlp;
  //pm.nlp_temp          = workspace->nlp_temp;
  //pm.Clp          = workspace->Clp;
  //pm.Delta_lp     = workspace->Delta_lp;
  //pm.dDelta_lp    = workspace->dDelta_lp;
  //pm.Delta_lp_temp  = workspace->Delta_lp_temp;
  //pm.dDelta_lp_temp    = workspace->dDelta_lp_temp;

  //int ntypes = pm.ntypes; 
  //int s1;
  //for(s1 = 0; s1 < ntypes; s1++)
  //{
  //  pm.spm[s1].valency = system->reax_param.sbp[s1].valency;
  //  pm.spm[s1].valency_e = system->reax_param.sbp[s1].valency_e;
  //  pm.spm[s1].valency_boc = system->reax_param.sbp[s1].valency_boc;
  //  pm.spm[s1].valency_val = system->reax_param.sbp[s1].valency_val;
  //  pm.spm[s1].nlp_opt = system->reax_param.sbp[s1].nlp_opt;
  //  pm.spm[s1].mass = system->reax_param.sbp[s1].mass;
  //}
  if(athread_idle() == 0)
    athread_init();
  athread_spawn(BO_Section3_C_Para, param);
  athread_join();

  #else
  for( j = 0; j < system->N; ++j )
  {
    type_j = system->packed_atoms[j].type;
    sbp_j = &(system->reax_param.sbp[ type_j ]);

    workspace->Delta[j]     = workspace->bo_dboc[j][0] - sbp_j->valency;
    workspace->Delta_e[j]   = workspace->bo_dboc[j][0] - sbp_j->valency_e;
    workspace->bo_dboc[j][1] = workspace->bo_dboc[j][0] - sbp_j->valency_boc;
    workspace->Delta_val[j] = workspace->bo_dboc[j][0] - sbp_j->valency_val;

    workspace->vlpex[j] = workspace->Delta_e[j] - 2.0 * (int)(workspace->Delta_e[j]/2.0);
    explp1              = exp(-p_lp1 * SQR(2.0 + workspace->vlpex[j]));
    workspace->nlp[j]       = explp1 - (int)(workspace->Delta_e[j] / 2.0);
    workspace->Delta_lp[j]  = sbp_j->nlp_opt - workspace->nlp[j];
    workspace->Clp[j]       = 2.0 * p_lp1 * explp1 * (2.0 + workspace->vlpex[j]);
    workspace->dDelta_lp[j] = workspace->Clp[j];

    if( sbp_j->mass > 21.0 ) 
    {
      workspace->nlp_temp[j]        = 0.5 * (sbp_j->valency_e - sbp_j->valency);
      workspace->Delta_lp_temp[j]   = sbp_j->nlp_opt - workspace->nlp_temp[j];
      workspace->dDelta_lp_temp[j]  = 0.;
    }
    else 
    {
      workspace->nlp_temp[j]        = workspace->nlp[j];
      workspace->Delta_lp_temp[j]   = sbp_j->nlp_opt - workspace->nlp_temp[j];
      workspace->dDelta_lp_temp[j]  = workspace->Clp[j];
    }

  }//for-j
  #endif

}


void BO_C(param_pack_t *param)
{
  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;

  int i, j, pj, type_i, type_j;
  int start_i, end_i, sym_index;
  double val_i, Deltap_i, Deltap_boc_i;
  double val_j, Deltap_j, Deltap_boc_j;
  double f1, f2, f3, f4, f5, f4f5, exp_f4, exp_f5;
  double exp_p1i,        exp_p2i, exp_p1j, exp_p2j;
  double temp, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
  double Cf45_ij, Cf45_ji, p_lp1; //u_ij, u_ji
  double A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji;
  double explp1, p_boc1, p_boc2;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  bond_order_data *bo_ij, *bo_ji;
  reax_list *bonds = (*lists) + BONDS;
  double *BO_list     = bonds->BO_list;
  rvec2 *BOpi_list    = bonds->BOpi_list;
  bond_order_data *bo_data_list = bonds->bo_data_list;


  p_boc1 = system->reax_param.gp.l[0];
  p_boc2 = system->reax_param.gp.l[1];

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  bo_param_pack_t pm;
  pm.system       = system;
  pm.control      = control;
  pm.data         = data;
  pm.workspace    = workspace;
  pm.lists        = lists;
  pm.out_control  = out_control;
  pm.packed_atoms = system->packed_atoms;
  pm.bond_list    = bonds->select.bond_list;
  pm.BO_list      = bonds->BO_list;
  pm.BOpi_list    = bonds->BOpi_list;
  pm.bo_data_list = bonds->bo_data_list;
  pm.bond_list    = bonds->select.bond_list;
  pm.index        = bonds->index;
  pm.end_index    = bonds->end_index;
  pm.bo_dboc      = workspace->bo_dboc;
  pm.Deltap       = workspace->Deltap;
  pm.Deltap_boc   = workspace->Deltap_boc;
  pm.N            = system->N;
  pm.ntypes       = system->reax_param.num_atom_types;
  pm.p_boc1       = p_boc1;
  pm.p_boc2       = p_boc2;
  pm.Delta        = workspace->Delta;
  pm.Delta_e      = workspace->Delta_e;
  pm.Delta_val      = workspace->Delta_val;
  pm.vlpex        = workspace->vlpex;
  pm.nlp          = workspace->nlp;
  pm.nlp_temp          = workspace->nlp_temp;
  pm.Clp          = workspace->Clp;
  pm.Delta_lp     = workspace->Delta_lp;
  pm.dDelta_lp    = workspace->dDelta_lp;
  pm.Delta_lp_temp  = workspace->Delta_lp_temp;
  pm.dDelta_lp_temp    = workspace->dDelta_lp_temp;


  int s1, s2;
  for(s1 = 0; s1 < pm.ntypes; s1++)
  {
    pm.sbp[s1] = system->reax_param.sbp[s1].valency;
  
    pm.spm[s1].valency = system->reax_param.sbp[s1].valency;
    pm.spm[s1].valency_e = system->reax_param.sbp[s1].valency_e;
    pm.spm[s1].valency_boc = system->reax_param.sbp[s1].valency_boc;
    pm.spm[s1].valency_val = system->reax_param.sbp[s1].valency_val;
    pm.spm[s1].nlp_opt = system->reax_param.sbp[s1].nlp_opt;
    pm.spm[s1].mass = system->reax_param.sbp[s1].mass;

    for(s2 = 0; s2 < pm.ntypes; s2++)
    {
      pm.tbp[s1 * pm.ntypes + s2].p_boc3  = system->reax_param.tbp[s1][s2].p_boc3;
      pm.tbp[s1 * pm.ntypes + s2].p_boc4  = system->reax_param.tbp[s1][s2].p_boc4;
      pm.tbp[s1 * pm.ntypes + s2].p_boc5  = system->reax_param.tbp[s1][s2].p_boc5;
      pm.tbp[s1 * pm.ntypes + s2].v13cor  = system->reax_param.tbp[s1][s2].v13cor;
      pm.tbp[s1 * pm.ntypes + s2].ovc     = system->reax_param.tbp[s1][s2].ovc;
    }
  }

  #define SECTION1
  #ifdef  SECTION1
  BO_Section1_C(&pm);
  #else
 /* Calculate Deltaprime, Deltaprime_boc values */
  for( i = 0; i < system->N; ++i ) 
  {
    type_i = system->packed_atoms[i].type;
    if (type_i < 0) continue;
    sbp_i = &(system->reax_param.sbp[type_i]);
    workspace->Deltap[i]           = workspace->bo_dboc[i][0] - sbp_i->valency;
    workspace->Deltap_boc[i]       = workspace->bo_dboc[i][0] - sbp_i->valency_val;
    workspace->bo_dboc[i][0] = 0;
  }
  #endif

  #define SECTION2
  #ifdef SECTION2
  BO_Section2_C(&pm);
  #else
  /* Corrected Bond Order calculations */
  for( i = 0; i < system->N; ++i ) 
  {
    type_i = system->packed_atoms[i].type;
    if (type_i < 0) continue;
    sbp_i = &(system->reax_param.sbp[type_i]);
    val_i = sbp_i->valency;
    Deltap_i     = workspace->Deltap[i];
    Deltap_boc_i = workspace->Deltap_boc[i];
    start_i = Start_Index(i, bonds);
    end_i   = End_Index(i, bonds);

    for( pj = start_i; pj < end_i; ++pj) 
    {
      j = bonds->select.bond_list[pj].nbr;
      type_j = system->packed_atoms[j].type;
      if (type_j < 0) continue;
      //bo_ij = &( bonds->select.bond_list[pj].bo_data );
      bo_ij = &(bo_data_list[pj]);
      // fprintf( stderr, "\tj:%d - ubo: %8.3f\n", j+1, bo_ij->BO );

      //if( i < j || workspace->bond_mark[j] > 3 ) 
      {
        twbp = &( system->reax_param.tbp[type_i][type_j] );

        if( twbp->ovc < 0.001 && twbp->v13cor < 0.001 ) 
        {
          bo_ij->C1dbo = 1.000000;
          bo_ij->C2dbo = 0.000000;
          bo_ij->C3dbo = 0.000000;

          //bo_ij->C1dbopi = bo_ij->BO_pi;
          bo_ij->C1dbopi = BOpi_list[pj][0];
          bo_ij->C2dbopi = 0.000000;
          bo_ij->C3dbopi = 0.000000;
          bo_ij->C4dbopi = 0.000000;

          //bo_ij->C1dbopi2 = bo_ij->BO_pi2;
          bo_ij->C1dbopi2 = BOpi_list[pj][1];
          bo_ij->C2dbopi2 = 0.000000;
          bo_ij->C3dbopi2 = 0.000000;
          bo_ij->C4dbopi2 = 0.000000;

        }
        else 
        {
          val_j = system->reax_param.sbp[type_j].valency;
          Deltap_j = workspace->Deltap[j];
          Deltap_boc_j = workspace->Deltap_boc[j];

          /* on page 1 */
          if( twbp->ovc >= 0.001 ) 
          {
            /* Correction for overcoordination */
            exp_p1i = exp( -p_boc1 * Deltap_i );
            exp_p2i = exp( -p_boc2 * Deltap_i );
            exp_p1j = exp( -p_boc1 * Deltap_j );
            exp_p2j = exp( -p_boc2 * Deltap_j );

            f2 = exp_p1i + exp_p1j;
            f3 = -1.0 / p_boc2 * log( 0.5 * ( exp_p2i  + exp_p2j ) );
            f1 = 0.5 * ( ( val_i + f2 )/( val_i + f2 + f3 ) +
                         ( val_j + f2 )/( val_j + f2 + f3 ) );

            temp = f2 + f3;
            u1_ij = val_i + temp;
            u1_ji = val_j + temp;
            Cf1A_ij = 0.5 * f3 * (1.0 / SQR( u1_ij ) +
                                  1.0 / SQR( u1_ji ));
            Cf1B_ij = -0.5 * (( u1_ij - f3 ) / SQR( u1_ij ) +
                              ( u1_ji - f3 ) / SQR( u1_ji ));

            Cf1_ij = 0.50 * ( -p_boc1 * exp_p1i / u1_ij -
                              ((val_i+f2) / SQR(u1_ij)) *
                              ( -p_boc1 * exp_p1i +
                                exp_p2i / ( exp_p2i + exp_p2j ) ) +
                              -p_boc1 * exp_p1i / u1_ji -
                              ((val_j+f2) / SQR(u1_ji)) *
                              ( -p_boc1 * exp_p1i +
                                exp_p2i / ( exp_p2i + exp_p2j ) ));


            Cf1_ji = -Cf1A_ij * p_boc1 * exp_p1j +
              Cf1B_ij * exp_p2j / ( exp_p2i + exp_p2j );

          }
          else 
          {
            /* No overcoordination correction! */
            f1 = 1.0;
            Cf1_ij = Cf1_ji = 0.0;
          }

          if( twbp->v13cor >= 0.001 ) 
          {
            /* Correction for 1-3 bond orders */
            //exp_f4 =exp(-(twbp->p_boc4 * SQR( bo_ij->BO ) -
            exp_f4 =exp(-(twbp->p_boc4 * SQR(BO_list[pj]) -
                          Deltap_boc_i) * twbp->p_boc3 + twbp->p_boc5);
            //exp_f5 =exp(-(twbp->p_boc4 * SQR( bo_ij->BO ) -
            exp_f5 =exp(-(twbp->p_boc4 * SQR( BO_list[pj]) -
                          Deltap_boc_j) * twbp->p_boc3 + twbp->p_boc5);

            f4 = 1. / (1. + exp_f4);
            f5 = 1. / (1. + exp_f5);
            f4f5 = f4 * f5;

            /* Bond Order pages 8-9, derivative of f4 and f5 */
            Cf45_ij = -f4 * exp_f4;
            Cf45_ji = -f5 * exp_f5;
          }
          else 
          {
            f4 = f5 = f4f5 = 1.0;
            Cf45_ij = Cf45_ji = 0.0;
          }

          /* Bond Order page 10, derivative of total bond order */
          A0_ij = f1 * f4f5;
          //A1_ij = -2 * twbp->p_boc3 * twbp->p_boc4 * bo_ij->BO * (Cf45_ij + Cf45_ji);
          A1_ij = -2 * twbp->p_boc3 * twbp->p_boc4 * BO_list[pj] * (Cf45_ij + Cf45_ji);
          A2_ij = Cf1_ij / f1 + twbp->p_boc3 * Cf45_ij;
          A2_ji = Cf1_ji / f1 + twbp->p_boc3 * Cf45_ji;
          A3_ij = A2_ij + Cf1_ij / f1;
          A3_ji = A2_ji + Cf1_ji / f1;

          /* find corrected bond orders and their derivative coef */
          //bo_ij->BO    = bo_ij->BO    * A0_ij;
          BO_list[pj]    = BO_list[pj]  * A0_ij;
          //bo_ij->BO_pi = bo_ij->BO_pi * A0_ij *f1;
          //bo_ij->BO_pi2= bo_ij->BO_pi2* A0_ij *f1;
          BOpi_list[pj][0]= BOpi_list[pj][0]* A0_ij *f1;
          BOpi_list[pj][1]= BOpi_list[pj][1]* A0_ij *f1;

          //bo_ij->BO_s  = bo_ij->BO - ( bo_ij->BO_pi + bo_ij->BO_pi2 );
          //bo_ij->BO_s  = BO_list[pj] - ( bo_ij->BO_pi + bo_ij->BO_pi2 );
          bo_ij->BO_s  = BO_list[pj] - (BOpi_list[pj][0] + BOpi_list[pj][1]);

          //bo_ij->C1dbo = A0_ij + bo_ij->BO * A1_ij;
          //bo_ij->C2dbo = bo_ij->BO * A2_ij;
          //bo_ij->C3dbo = bo_ij->BO * A2_ji;
          bo_ij->C1dbo = A0_ij + BO_list[pj] * A1_ij;
          bo_ij->C2dbo = BO_list[pj] * A2_ij;
          bo_ij->C3dbo = BO_list[pj] * A2_ji;


          bo_ij->C1dbopi = f1*f1*f4*f5;
          //bo_ij->C2dbopi = bo_ij->BO_pi * A1_ij;
          //bo_ij->C3dbopi = bo_ij->BO_pi * A3_ij;
          //bo_ij->C4dbopi = bo_ij->BO_pi * A3_ji;
          bo_ij->C2dbopi = BOpi_list[pj][0] * A1_ij;
          bo_ij->C3dbopi = BOpi_list[pj][0] * A3_ij;
          bo_ij->C4dbopi = BOpi_list[pj][0] * A3_ji;


          bo_ij->C1dbopi2 = f1*f1*f4*f5;
          //bo_ij->C2dbopi2 = bo_ij->BO_pi2 * A1_ij;
          //bo_ij->C3dbopi2 = bo_ij->BO_pi2 * A3_ij;
          //bo_ij->C4dbopi2 = bo_ij->BO_pi2 * A3_ji;
          bo_ij->C2dbopi2 = BOpi_list[pj][1] * A1_ij;
          bo_ij->C3dbopi2 = BOpi_list[pj][1] * A3_ij;
          bo_ij->C4dbopi2 = BOpi_list[pj][1] * A3_ji;

        }

        /* neglect bonds that are < 1e-10 */
        //if( bo_ij->BO < 1e-10 ) bo_ij->BO = 0.0;
        if(BO_list[pj] < 1e-10 ) BO_list[pj] = 0.0;

        if( bo_ij->BO_s < 1e-10 )
          bo_ij->BO_s = 0.0;
        
        //if( bo_ij->BO_pi < 1e-10 )
        //  bo_ij->BO_pi = 0.0;
        //if( bo_ij->BO_pi2 < 1e-10 )
        //  bo_ij->BO_pi2 = 0.0;
        if(BOpi_list[pj][0] < 1e-10 )
          BOpi_list[pj][0] = 0.0;
        if(BOpi_list[pj][1] < 1e-10 )
          BOpi_list[pj][1] = 0.0;

        //workspace->total_bond_order[i] += bo_ij->BO; //now keeps total_BO
        //workspace->bo_dboc[i][0] += bo_ij->BO; //now keeps total_BO
        workspace->bo_dboc[i][0] += BO_list[pj]; //now keeps total_BO

      }
      //else 
      //{
      //  /* We only need to update bond orders from bo_ji
      //     everything else is set in uncorrected_bo calculations */
      //  sym_index     = bonds->select.bond_list[pj].sym_index;
      //  bo_ji         = &(bonds->select.bond_list[ sym_index ].bo_data);
      //  bo_ij->BO     = bo_ji->BO;
      //  bo_ij->BO_s   = bo_ji->BO_s;
      //  bo_ij->BO_pi  = bo_ji->BO_pi;
      //  bo_ij->BO_pi2 = bo_ji->BO_pi2;

      //  workspace->total_bond_order[i] += bo_ij->BO;// now keeps total_BO
      //}
    }
  }
  #endif


  #define SECTION3
  #ifdef SECTION3
  BO_Section3_C(&pm);
  #else

  p_lp1 = system->reax_param.gp.l[15];
  for( j = 0; j < system->N; ++j )
  {
    type_j = system->packed_atoms[j].type;
    sbp_j = &(system->reax_param.sbp[ type_j ]);

    workspace->Delta[j]     = workspace->bo_dboc[j][0] - sbp_j->valency;
    workspace->Delta_e[j]   = workspace->bo_dboc[j][0] - sbp_j->valency_e;
    workspace->bo_dboc[j][1] = workspace->bo_dboc[j][0] - sbp_j->valency_boc;
    workspace->Delta_val[j] = workspace->bo_dboc[j][0] - sbp_j->valency_val;

    workspace->vlpex[j] = workspace->Delta_e[j] - 2.0 * (int)(workspace->Delta_e[j]/2.0);
    explp1              = exp(-p_lp1 * SQR(2.0 + workspace->vlpex[j]));
    workspace->nlp[j]       = explp1 - (int)(workspace->Delta_e[j] / 2.0);
    workspace->Delta_lp[j]  = sbp_j->nlp_opt - workspace->nlp[j];
    workspace->Clp[j]       = 2.0 * p_lp1 * explp1 * (2.0 + workspace->vlpex[j]);
    workspace->dDelta_lp[j] = workspace->Clp[j];

    if( sbp_j->mass > 21.0 ) 
    {
      workspace->nlp_temp[j]        = 0.5 * (sbp_j->valency_e - sbp_j->valency);
      workspace->Delta_lp_temp[j]   = sbp_j->nlp_opt - workspace->nlp_temp[j];
      workspace->dDelta_lp_temp[j]  = 0.;
    }
    else 
    {
      workspace->nlp_temp[j]        = workspace->nlp[j];
      workspace->Delta_lp_temp[j]   = sbp_j->nlp_opt - workspace->nlp_temp[j];
      workspace->dDelta_lp_temp[j]  = workspace->Clp[j];
    }

  }//for-j
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
//#define LWPF_UNIT U(BOFUNC)
//#define LWPF_KERNELS K(ALL) K(BEFORE) K(CMP) K(LOG)   
//#include "lwpf2.h"

//read cache;
#define READ_C_H    9
#define READ_C_S    4
#define READ_C_LSZ  (1 << READ_C_S)
#define READ_C_LCNT (1 << (READ_C_H - READ_C_S))
#define READ_C_MM   (READ_C_LSZ - 1)
#define READ_C_LM   (READ_C_LCNT - 1)
#define KSTEP 128
void BO_Section1_C_Para(bo_param_pack_t *param)
{
  dma_init();
  bo_param_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(bo_param_pack_t));
  dma_syn();

  //reax_system_c *system         = param->system;
  //storage *workspace            = param->workspace;

  double Deltap[KSTEP];
  double Deltap_boc[KSTEP];
  rvec2 bo_dboc[KSTEP];
  atom_pack_t atom_i[KSTEP];

  int i, j, pj, type_i, type_j;
  //single_body_parameters *sbp_i, *sbp_j;
  bo_sbp_t *sbp_i;
  
  int ist, ied, isz, ioff;
  int N = l_pm.N;
  for( ist = _MYID*KSTEP; ist < N; ist += KSTEP*64) 
  {
    ied = ist + KSTEP;
    if(ied > N)
      ied = N;
    isz = ied - ist;

    pe_get(l_pm.Deltap+ist,       Deltap,     sizeof(double)*isz);
    pe_get(l_pm.Deltap_boc+ist,   Deltap_boc, sizeof(double)*isz);
    pe_get(l_pm.bo_dboc+ist,      bo_dboc,    sizeof(rvec2)*isz);
    pe_get(l_pm.packed_atoms+ist, atom_i,     sizeof(atom_pack_t)*isz);
    dma_syn();

    //for( i = 0; i < system->N; ++i ) 
    for(i = ist; i < ied; i++) 
    {
      ioff = i - ist;
      //type_i = system->packed_atoms[i].type;
      type_i = atom_i[ioff].type;//system->packed_atoms[i].type;
      if (type_i < 0) continue;
      //sbp_i = &(system->reax_param.sbp[type_i]);
      sbp_i = &l_pm.spm[type_i];
      Deltap[ioff]           = bo_dboc[ioff][0] - sbp_i->valency;
      Deltap_boc[ioff]       = bo_dboc[ioff][0] - sbp_i->valency_val;
      bo_dboc[ioff][0] = 0;
    }
    pe_put(l_pm.Deltap+ist,     Deltap,     sizeof(double)*isz);
    pe_put(l_pm.Deltap_boc+ist, Deltap_boc, sizeof(double)*isz);
    pe_put(l_pm.bo_dboc+ist,    bo_dboc,    sizeof(rvec2)*isz);
    dma_syn();
  }

}
void read_cache(int i,
                atom_pack_t atoms_cache[][READ_C_LSZ], int *read_ctag,
                atom_pack_t *packed_atoms, atom_pack_t *atom_i,
                double Deltap_cache[][READ_C_LSZ], double *Deltap, double *Deltap_i,
                double Deltap_boc_cache[][READ_C_LSZ],   double *Deltap_boc, double *Deltap_boc_i,
                int idx_cache[][READ_C_LSZ], int *index, int *i_st,
                int ed_idx_cache[][READ_C_LSZ], int *end_index, int *i_ed)
{
  dma_init();
  if (read_ctag[(i >> READ_C_S) & READ_C_LM] != i >> READ_C_S)
  {
    pe_get(packed_atoms + (i & ~READ_C_MM), 
          atoms_cache[(i >> READ_C_S) & READ_C_LM], 
          sizeof(atom_pack_t) * READ_C_LSZ);
    pe_get(Deltap + (i & ~READ_C_MM), Deltap_cache[(i >> READ_C_S) & READ_C_LM], 
          sizeof(double) * READ_C_LSZ);
    pe_get(Deltap_boc + (i & ~READ_C_MM), Deltap_boc_cache[(i >> READ_C_S) & READ_C_LM], 
          sizeof(double) * READ_C_LSZ);
    pe_get(index + (i & ~READ_C_MM), idx_cache[(i >> READ_C_S) & READ_C_LM], 
          sizeof(int) * READ_C_LSZ);
    pe_get(end_index + (i & ~READ_C_MM), ed_idx_cache[(i >> READ_C_S) & READ_C_LM], 
          sizeof(int) * READ_C_LSZ);

    dma_syn();
    read_ctag[(i >> READ_C_S) & READ_C_LM] = i >> READ_C_S;
  }
  *atom_i       = atoms_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
  *Deltap_i     = Deltap_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
  *Deltap_boc_i = Deltap_boc_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
  *i_st         = idx_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
  *i_ed         = ed_idx_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
}


void BO_Section2_C_Para(bo_param_pack_t *param)
{
  //lwpf_enter(BOFUNC);
  //lwpf_start(ALL);
  //lwpf_start(BEFORE);

  dma_init();
  bo_param_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(bo_param_pack_t));
  dma_syn();

  reax_system_c *system         = l_pm.system;
  control_params *control       = l_pm.control;
  simulation_data *data         = l_pm.data;
  storage *workspace            = l_pm.workspace;
  reax_list **lists             = l_pm.lists;
  output_controls *out_control  = l_pm.out_control;

  int i, j, pj, pj_off, type_i, type_j, pk, pk_off;
  int start_i, end_i, len_i, sym_index;
  double val_i, Deltap_i, Deltap_boc_i;
  double val_j, Deltap_j, Deltap_boc_j;
  double f1, f2, f3, f4, f5, f4f5, exp_f4, exp_f5;
  double exp_p1i,        exp_p2i, exp_p1j, exp_p2j;
  double temp, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
  double Cf45_ij, Cf45_ji, p_lp1; //u_ij, u_ji
  double A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji;
  double explp1, p_boc1, p_boc2;
  single_body_parameters *sbp_i, *sbp_j;
  //two_body_parameters *twbp;
  bo_tbp_t *twbp;
  bond_order_data *bo_ij, *bo_ji;
  //reax_list *bonds = (*lists) + BONDS;
  //double *BO_list     = bonds->BO_list;
  //rvec2 *BOpi_list    = bonds->BOpi_list;
  //bond_order_data *bo_data_list = bonds->bo_data_list;
  //p_boc1 = system->reax_param.gp.l[0];
  //p_boc2 = system->reax_param.gp.l[1];

  atom_pack_t *packed_atoms     = l_pm.packed_atoms;
  //double *Deltap_boc            = l_pm.Deltap_boc;
  //double *Deltap                = l_pm.Deltap;
  //rvec2 *bo_dboc                = l_pm.bo_dboc;
  double *BO_list               = l_pm.BO_list;
  rvec2 *BOpi_list              = l_pm.BOpi_list;
  bond_order_data *bo_data_list = l_pm.bo_data_list;
  bond_data *bond_list          = l_pm.bond_list;
  //int *index                    = l_pm.index;
  //int *end_index                = l_pm.end_index;
  int N                         = l_pm.N;
  int ntypes                    = l_pm.ntypes;
  p_boc1                        = l_pm.p_boc1;
  p_boc2                        = l_pm.p_boc2;

  
  atom_pack_t atom_i[ISTEP];
  int index[ISTEP], end_index[ISTEP];
  double Deltap_boc[ISTEP], Deltap[ISTEP];
  rvec2 bo_dboc[ISTEP];
  double BO_list_i[MAX_LEN];
  rvec2 BOpi_list_i[MAX_LEN];
  bond_order_data bo_data_i[MAX_LEN];
  bond_data bond_list_i[MAX_LEN];
  bond_data bond_list_j[MAX_LEN];

  //read_cache;
  int read_ctag[READ_C_LCNT];
  atom_pack_t atoms_cache[READ_C_LCNT][READ_C_LSZ];
  atom_pack_t atom_j;
  double Deltap_cache[READ_C_LCNT][READ_C_LSZ];
  double Deltap_boc_cache[READ_C_LCNT][READ_C_LSZ];
  double Deltap_cj, Deltap_boc_cj;
  int idx_cache[READ_C_LCNT][READ_C_LSZ];
  int ed_idx_cache[READ_C_LCNT][READ_C_LSZ];
  int start_j, end_j;
  for(i = 0; i < READ_C_LCNT; i++)
  {
    read_ctag[i] = -1;
  }
  //lwpf_stop(BEFORE);

  //lwpf_start(CMP);
  int ist, ied, isz, ioff;
  //if(_MYID == 10)
  //for(i = 0; i < system->N; i++) 
  for( ist = _MYID*ISTEP; ist < N; ist += ISTEP*64) 
  {
    ied = ist + ISTEP;
    if(ied > N)
      ied = N;
    isz = ied - ist;

    pe_get(l_pm.index+ist,        index,      sizeof(int)*isz);
    pe_get(l_pm.end_index+ist,    end_index,  sizeof(int)*isz);
    pe_get(l_pm.Deltap_boc+ist,   Deltap_boc, sizeof(double)*isz);
    pe_get(l_pm.Deltap+ist,       Deltap,     sizeof(double)*isz);
    pe_get(l_pm.bo_dboc+ist,      bo_dboc,    sizeof(rvec2)*isz);
    pe_get(l_pm.packed_atoms+ist, atom_i,     sizeof(atom_pack_t)*isz);
    dma_syn();

    for( i = ist; i < ied; ++i) 
    {
      ioff = i - ist;
      //type_i = system->packed_atoms[i].type;
      //type_i = packed_atoms[i].type;
      type_i = atom_i[ioff].type;
      if (type_i < 0) continue;
      //sbp_i = &(system->reax_param.sbp[type_i]);
      //val_i = sbp_i->valency;
      val_i = l_pm.sbp[type_i];
      //Deltap_i     = workspace->Deltap[i];
      //Deltap_boc_i = workspace->Deltap_boc[i];
      Deltap_i     = Deltap[ioff];
      Deltap_boc_i = Deltap_boc[ioff];

      //start_i = Start_Index(i, bonds);
      //end_i   = End_Index(i, bonds);
      //start_i = l_pm.index[i];//Start_Index(i, bonds);
      //end_i   = l_pm.end_index[i];//End_Index(i, bonds);
      start_i = index[ioff];//l_pm.index[i];//Start_Index(i, bonds);
      end_i   = end_index[ioff];//l_pm.end_index[i];//End_Index(i, bonds);
      len_i = end_i - start_i;
     
      if(len_i <=  0)  continue;
      pe_get(l_pm.BO_list+start_i,        BO_list_i,    sizeof(double)*len_i);
      pe_get(l_pm.BOpi_list+start_i,      BOpi_list_i,  sizeof(rvec2)*len_i);
      pe_get(l_pm.bo_data_list+start_i,   bo_data_i,    sizeof(bond_order_data)*len_i);
      pe_get(l_pm.bond_list+start_i, bond_list_i,  sizeof(bond_data)*len_i);
      dma_syn();

      for( pj = start_i; pj < end_i; ++pj) 
      {
        pj_off = pj - start_i;
        //j = bonds->select.bond_list[pj].nbr;
        j = bond_list_i[pj_off].nbr;
        read_cache(j, atoms_cache, read_ctag, packed_atoms, &atom_j,
                    Deltap_cache, l_pm.Deltap, &Deltap_cj,
                    Deltap_boc_cache, l_pm.Deltap_boc, &Deltap_boc_cj,
                    idx_cache, l_pm.index, &start_j,
                    ed_idx_cache,l_pm.end_index, &end_j);

        type_j = atom_j.type;
        //type_j = packed_atoms[j].type;
        if (type_j < 0) continue;
        bo_ij = &(bo_data_i[pj_off]);
        
        /*compute for sym_index*/
        //int start_j = l_pm.index[j];//Start_Index( j, bonds);
        //int end_j   = l_pm.end_index[j];//End_Index( j, bonds );
        int len_j = end_j - start_j;
        pe_get(l_pm.bond_list+start_j, bond_list_j,  sizeof(bond_data)*len_j);
        dma_syn();

        for ( pk = start_j; pk < end_j; ++pk ) 
        {
          pk_off = pk - start_j;
          int k = bond_list_j[pk_off].nbr;

          //bond_data *kbond = l_pm.bond_list + pk;
          //int k = kbond->nbr;
          if (k == i)
          {
            //l_pm.bond_list[pj].sym_index = pk;
            bond_list_i[pj_off].sym_index = pk;
          }
        }//for
        /*end of sym_index*/

        //twbp = &( system->reax_param.tbp[type_i][type_j] );
        twbp = &(l_pm.tbp[type_i*ntypes+type_j] );

        if( twbp->ovc < 0.001 && twbp->v13cor < 0.001 ) 
        {
          //bo_ij->C1dbo = 1.000000;
          //bo_ij->C2dbo = 0.000000;
          bo_ij->C3dbo = 0.000000;

          //bo_ij->C1dbopi = BOpi_list[pj][0];
          bo_ij->C1dbopi = BOpi_list_i[pj_off][0];
          bo_ij->C2dbopi = 0.000000;
          bo_ij->C3dbopi = 0.000000;
          bo_ij->C4dbopi = 0.000000;

          //bo_ij->C1dbopi2 = BOpi_list[pj][1];
          bo_ij->C1dbopi2 = BOpi_list_i[pj_off][1];
          bo_ij->C2dbopi2 = 0.000000;
          bo_ij->C3dbopi2 = 0.000000;
          bo_ij->C4dbopi2 = 0.000000;
        }
        else 
        {
          //val_j = system->reax_param.sbp[type_j].valency;
          val_j = l_pm.sbp[type_j];
          //Deltap_j = workspace->Deltap[j];
          //Deltap_boc_j = workspace->Deltap_boc[j];
          //Deltap_j      = l_pm.Deltap[j];
          //Deltap_boc_j  = l_pm.Deltap_boc[j];
          Deltap_j      = Deltap_cj;
          Deltap_boc_j  = Deltap_boc_cj;


          /* on page 1 */
          if( twbp->ovc >= 0.001 ) 
          {
            /* Correction for overcoordination */
            exp_p1i = p_expd( -p_boc1 * Deltap_i );
            exp_p2i = p_expd( -p_boc2 * Deltap_i );
            exp_p1j = p_expd( -p_boc1 * Deltap_j );
            exp_p2j = p_expd( -p_boc2 * Deltap_j );

            f2 = exp_p1i + exp_p1j;
            //lwpf_start(LOG);
            f3 = -1.0 / p_boc2 * log( 0.5 * ( exp_p2i  + exp_p2j ) );
            f1 = 0.5 * ( ( val_i + f2 )/( val_i + f2 + f3 ) +
                         ( val_j + f2 )/( val_j + f2 + f3 ) );

            //lwpf_stop(LOG);

            temp = f2 + f3;
            u1_ij = val_i + temp;
            u1_ji = val_j + temp;
            Cf1A_ij = 0.5 * f3 * (1.0 / SQR( u1_ij ) +
                                  1.0 / SQR( u1_ji ));
            Cf1B_ij = -0.5 * (( u1_ij - f3 ) / SQR( u1_ij ) +
                              ( u1_ji - f3 ) / SQR( u1_ji ));

            Cf1_ij = 0.50 * ( -p_boc1 * exp_p1i / u1_ij -
                              ((val_i+f2) / SQR(u1_ij)) *
                              ( -p_boc1 * exp_p1i +
                                exp_p2i / ( exp_p2i + exp_p2j ) ) +
                              -p_boc1 * exp_p1i / u1_ji -
                              ((val_j+f2) / SQR(u1_ji)) *
                              ( -p_boc1 * exp_p1i +
                                exp_p2i / ( exp_p2i + exp_p2j ) ));


            Cf1_ji = -Cf1A_ij * p_boc1 * exp_p1j +
              Cf1B_ij * exp_p2j / ( exp_p2i + exp_p2j );

          }
          else 
          {
            /* No overcoordination correction! */
            f1 = 1.0;
            Cf1_ij = Cf1_ji = 0.0;
          }

          if( twbp->v13cor >= 0.001 ) 
          {
            /* Correction for 1-3 bond orders */
            exp_f4 =p_expd(-(twbp->p_boc4 * SQR(BO_list_i[pj_off]) -
                       Deltap_boc_i) * twbp->p_boc3 + twbp->p_boc5);
            exp_f5 =p_expd(-(twbp->p_boc4 * SQR(BO_list_i[pj_off]) -
                          Deltap_boc_j) * twbp->p_boc3 + twbp->p_boc5);

            f4 = 1. / (1. + exp_f4);
            f5 = 1. / (1. + exp_f5);
            f4f5 = f4 * f5;

            /* Bond Order pages 8-9, derivative of f4 and f5 */
            Cf45_ij = -f4 * exp_f4;
            Cf45_ji = -f5 * exp_f5;
          }
          else 
          {
            f4 = f5 = f4f5 = 1.0;
            Cf45_ij = Cf45_ji = 0.0;
          }

          /* Bond Order page 10, derivative of total bond order */
          A0_ij = f1 * f4f5;
          A1_ij = -2 * twbp->p_boc3 * twbp->p_boc4 * BO_list_i[pj_off] * (Cf45_ij + Cf45_ji);
          A2_ij = Cf1_ij / f1 + twbp->p_boc3 * Cf45_ij;
          A2_ji = Cf1_ji / f1 + twbp->p_boc3 * Cf45_ji;
          A3_ij = A2_ij + Cf1_ij / f1;
          A3_ji = A2_ji + Cf1_ji / f1;

          /* find corrected bond orders and their derivative coef */
          BO_list_i[pj_off]    = BO_list_i[pj_off]  * A0_ij;
          BOpi_list_i[pj_off][0]= BOpi_list_i[pj_off][0]* A0_ij *f1;
          BOpi_list_i[pj_off][1]= BOpi_list_i[pj_off][1]* A0_ij *f1;

          bo_ij->BO_s  = BO_list_i[pj_off] - (BOpi_list_i[pj_off][0] + BOpi_list_i[pj_off][1]);
          bo_ij->C1dbo = A0_ij + BO_list_i[pj_off] * A1_ij;
          bo_ij->C2dbo = BO_list_i[pj_off] * A2_ij;
          bo_ij->C3dbo = BO_list_i[pj_off] * A2_ji;

          bo_ij->C1dbopi = f1*f1*f4*f5;
          bo_ij->C2dbopi = BOpi_list_i[pj_off][0] * A1_ij;
          bo_ij->C3dbopi = BOpi_list_i[pj_off][0] * A3_ij;
          bo_ij->C4dbopi = BOpi_list_i[pj_off][0] * A3_ji;

          bo_ij->C1dbopi2 = f1*f1*f4*f5;
          bo_ij->C2dbopi2 = BOpi_list_i[pj_off][1] * A1_ij;
          bo_ij->C3dbopi2 = BOpi_list_i[pj_off][1] * A3_ij;
          bo_ij->C4dbopi2 = BOpi_list_i[pj_off][1] * A3_ji;
        }

        /* neglect bonds that are < 1e-10 */
        if(BO_list_i[pj_off] < 1e-10 ) BO_list_i[pj_off] = 0.0;

        if( bo_ij->BO_s < 1e-10 )
          bo_ij->BO_s = 0.0;
        
        if(BOpi_list_i[pj_off][0] < 1e-10 )
          BOpi_list_i[pj_off][0] = 0.0;
        if(BOpi_list_i[pj_off][1] < 1e-10 )
          BOpi_list_i[pj_off][1] = 0.0;

        //workspace->bo_dboc[i][0] += BO_list[pj]; //now keeps total_BO
        //l_pm.bo_dboc[i][0] += BO_list[pj]; //now keeps total_BO
        bo_dboc[ioff][0] += BO_list_i[pj_off]; //now keeps total_BO
      }//for-pj
      pe_put(l_pm.bond_list+start_i, bond_list_i,  sizeof(bond_data)*len_i);///
      pe_put(l_pm.BO_list+start_i,      BO_list_i,    sizeof(double)*len_i);
      pe_put(l_pm.BOpi_list+start_i,    BOpi_list_i,  sizeof(rvec2)*len_i);
      pe_put(l_pm.bo_data_list+start_i, bo_data_i,    sizeof(bond_order_data)*len_i);
      dma_syn();

    }//for-i
    pe_put(l_pm.bo_dboc+ist, bo_dboc, sizeof(rvec2)*isz);
    dma_syn();
  }//for-ist
  //lwpf_stop(CMP);

  //lwpf_stop(ALL);
  //lwpf_exit(BOFUNC);
}

#define JSTEP 64
void BO_Section3_C_Para(bo_param_pack_t *param)
{
  dma_init();
  bo_param_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(bo_param_pack_t));
  dma_syn();

  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;


  int i, j, pj, type_i, type_j;
  int start_i, end_i, sym_index;
  double val_i, Deltap_i, Deltap_boc_i;
  double val_j, Deltap_j, Deltap_boc_j;
  double f1, f2, f3, f4, f5, f4f5, exp_f4, exp_f5;
  double exp_p1i,        exp_p2i, exp_p1j, exp_p2j;
  double temp, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
  double Cf45_ij, Cf45_ji, p_lp1; //u_ij, u_ji
  double A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji;
  double explp1, p_boc1, p_boc2;
  //single_body_parameters *sbp_i, *sbp_j;
  bo_sbp_t *sbp_j;
  two_body_parameters *twbp;
  bond_order_data *bo_ij, *bo_ji;
  reax_list *bonds = (*lists) + BONDS;
  double *BO_list     = bonds->BO_list;
  rvec2 *BOpi_list    = bonds->BOpi_list;
  bond_order_data *bo_data_list = bonds->bo_data_list;
  
  p_lp1 = system->reax_param.gp.l[15];
  int jst, jed, jsz, joff;
  int N = system->N;
  
  rvec2 bo_dboc_j[JSTEP];
  double Delta[JSTEP], Delta_e[JSTEP], Delta_val[JSTEP];
  double vlpex[JSTEP], nlp[JSTEP], Clp[JSTEP];
  double nlp_temp[JSTEP];
  double Delta_lp[JSTEP], dDelta_lp[JSTEP];
  double Delta_lp_temp[JSTEP], dDelta_lp_temp[JSTEP];
  atom_pack_t atom_j[JSTEP];

  //for( j = 0; j < system->N; ++j )
  for(jst = _MYID * JSTEP; jst < N; jst+=JSTEP*64)
  {
    jed = jst + JSTEP;
    if(jed > N)
      jed = N;
    jsz = jed - jst;
    pe_get(l_pm.bo_dboc+jst,        bo_dboc_j,      sizeof(rvec2) *jsz);
    pe_get(l_pm.Delta+jst,          Delta,          sizeof(double) *jsz);
    pe_get(l_pm.Delta_e+jst,        Delta_e,        sizeof(double) *jsz);
    pe_get(l_pm.Delta_val+jst,      Delta_val,      sizeof(double) *jsz);
    pe_get(l_pm.vlpex+jst,          vlpex,          sizeof(double) *jsz);
    pe_get(l_pm.nlp+jst,            nlp,            sizeof(double) *jsz);
    pe_get(l_pm.nlp_temp+jst,       nlp_temp,       sizeof(double) *jsz);
    pe_get(l_pm.Clp+jst,            Clp,            sizeof(double) *jsz);
    pe_get(l_pm.Delta_lp+jst,       Delta_lp,       sizeof(double) *jsz);
    pe_get(l_pm.dDelta_lp+jst,      dDelta_lp,      sizeof(double) *jsz);
    pe_get(l_pm.Delta_lp_temp+jst,  Delta_lp_temp,  sizeof(double) *jsz);
    pe_get(l_pm.dDelta_lp_temp+jst, dDelta_lp_temp, sizeof(double) *jsz);
    pe_get(l_pm.packed_atoms+jst,   atom_j,         sizeof(atom_pack_t)*jsz);

    dma_syn();

    for(j = jst; j < jed; j++)
    {
      joff = j - jst;
      //type_j = system->packed_atoms[j].type;
      type_j = atom_j[joff].type;//system->packed_atoms[j].type;
      //sbp_j = &(system->reax_param.sbp[ type_j ]);
      sbp_j = &(l_pm.spm[type_j]);

      Delta[joff]         = bo_dboc_j[joff][0] - sbp_j->valency;
      Delta_e[joff]       = bo_dboc_j[joff][0] - sbp_j->valency_e;
      bo_dboc_j[joff][1]  = bo_dboc_j[joff][0] - sbp_j->valency_boc;
      Delta_val[joff]     = bo_dboc_j[joff][0] - sbp_j->valency_val;

      vlpex[joff]     = Delta_e[joff] - 2.0 * (int)(Delta_e[joff]/2.0);
      explp1          = p_expd(-p_lp1 * SQR(2.0 + vlpex[joff]));
      nlp[joff]       = explp1 - (int)(Delta_e[joff] / 2.0);
      Delta_lp[joff]  = sbp_j->nlp_opt - nlp[joff];
      Clp[joff]       = 2.0 * p_lp1 * explp1 * (2.0 + vlpex[joff]);
      dDelta_lp[joff] = Clp[joff];

      if( sbp_j->mass > 21.0 ) 
      {
        nlp_temp[joff]        = 0.5 * (sbp_j->valency_e - sbp_j->valency);
        Delta_lp_temp[joff]   = sbp_j->nlp_opt - nlp_temp[joff];
        dDelta_lp_temp[joff]  = 0.;
      }
      else 
      {
        nlp_temp[joff]        = nlp[joff];
        Delta_lp_temp[joff]   = sbp_j->nlp_opt - nlp_temp[joff];
        dDelta_lp_temp[joff]  = Clp[joff];
      }
    }
    pe_put(l_pm.bo_dboc+jst,        bo_dboc_j,      sizeof(rvec2) *jsz);
    pe_put(l_pm.Delta+jst,          Delta,          sizeof(double) *jsz);
    pe_put(l_pm.Delta_e+jst,        Delta_e,        sizeof(double) *jsz);
    pe_put(l_pm.Delta_val+jst,      Delta_val,      sizeof(double) *jsz);
    pe_put(l_pm.vlpex+jst,          vlpex,          sizeof(double) *jsz);
    pe_put(l_pm.nlp+jst,            nlp,            sizeof(double) *jsz);
    pe_put(l_pm.Clp+jst,            Clp,            sizeof(double) *jsz);
    pe_put(l_pm.Delta_lp+jst,       Delta_lp,       sizeof(double) *jsz);
    pe_put(l_pm.dDelta_lp+jst,      dDelta_lp,      sizeof(double) *jsz);
    pe_put(l_pm.Delta_lp_temp+jst,  Delta_lp_temp,  sizeof(double) *jsz);
    pe_put(l_pm.dDelta_lp_temp+jst, dDelta_lp_temp, sizeof(double) *jsz);
    dma_syn();

  }//for-j

}
#endif
