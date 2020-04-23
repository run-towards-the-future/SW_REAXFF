#include <math.h>
#include "reaxc_defs_sunway.h"
#include "reaxc_ctypes_sunway.h"
#include "reaxc_inlines_sw64.h"
#include "reaxc_hydrogen_bonds_sw64.h"
#include "stdio.h"
#include "stdlib.h"
#include "gptl.h"
#include "sunway.h"
#include "simd.h"
#include "sleef_math.h"
#define JSTEP 32
#define BLEN 20
#define HLEN 64

#ifdef MPE
#include <mpi.h>
#include <athread.h>
//#define LWPF_UNITS U(HBOND)
//#include "lwpf2.h"

extern SLAVE_FUN(Hydrogen_Bonds_C_ev0_cpe)(hb_param_pack_t *);
extern SLAVE_FUN(Hydrogen_Bonds_C_ev1_cpe)(hb_param_pack_t *);

//swcache_lock_t *locks_frc = NULL;
void Hydrogen_Bonds_C(hb_param_pack_t *param)
{
  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;

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

  #define SW_HYDROGEN
  #ifdef  SW_HYDROGEN
  hb_param_pack_t pm;
  pm.system       = system;
  pm.control      = control;
  pm.data         = data;
  pm.workspace    = workspace;
  pm.lists        = lists;
  
  double  packed_eng[4];
  packed_eng[0] = packed_eng[1] = packed_eng[2] = packed_eng[3] = 0;
  int nall = bonds->num_intrs;
  //if(locks_frc == NULL)
  //  locks_frc = swcache_u_prepare_locks(nall, 3);

  //pm.locks_frc      = locks_frc;
  pm.locks_frc      = bonds->locks_frc;

  pm.packed_eng     = packed_eng;
  pm.n              = system->n;
  pm.bond_list      = bond_list;
  pm.hbond_list     = hbond_list;
  pm.BO_list        = BO_list;
  pm.Cdbo_list      = Cdbo_list;
  pm.packed_atoms   = system->packed_atoms;
  pm.Hindex         = system->Hindex;
  pm.fCdDelta       = workspace->fCdDelta;
  pm.bindex         = bonds->index;
  pm.end_bindex     = bonds->end_index;
  pm.hbindex        = hbonds->index;
  pm.end_hbindex    = hbonds->end_index;
  pm.ntypes         = system->reax_param.num_atom_types;

  int ntypes = pm.ntypes;
  int s1, s2, s3, hsum = 0;
  for(s1 = 0; s1 < ntypes; s1++)
  {
    pm.sbp[s1] = system->reax_param.sbp[s1].p_hbond;
    for(s2 = 0; s2 < ntypes; s2++)
    {
      for(s3 = 0; s3 < ntypes; s3++)
      {
        pm.hbp[hsum] = system->reax_param.hbp[s1][s2][s3];
        hsum++;
      }
    }
  }

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
  
  if(system->evflag)
    athread_spawn(Hydrogen_Bonds_C_ev1_cpe, &pm);
  else
    athread_spawn(Hydrogen_Bonds_C_ev0_cpe, &pm);

  athread_join();

  //if(myrank == 0)
  //{
  //  lwpf_report_summary(stdout, &conf);
  //}

  system->eng_vdwl += packed_eng[0];
  data->my_en.e_hb += packed_eng[1];
  //if(locks_frc != NULL)
  //{
  //  free(locks_frc);
  //  locks_frc = NULL;
  //}
  return;

  #else

  int jst, jed, jsz, joff;
  //for( j = 0; j < system->n; ++j )
  for( jst = 0; jst < system->n; jst+=JSTEP)
  {
    jed = jst + JSTEP;
    if(jed > system->n)
      jed = system->n;
    jsz = jed - jst;

    for(j = jst; j < jed; j++)
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

              // dcos terms
              rvec4_ScaledAdd( workspace->fCdDelta[i], +CEhb2, dcos_theta_di );
              rvec4_ScaledAdd( workspace->fCdDelta[j], +CEhb2, dcos_theta_dj );
              rvec4_ScaledAdd( workspace->fCdDelta[k], +CEhb2, dcos_theta_dk );
              // dr terms
              rvec4_ScaledAdd( workspace->fCdDelta[j], -CEhb3/r_jk, dvec_jk );
              rvec4_ScaledAdd( workspace->fCdDelta[k], +CEhb3/r_jk, dvec_jk );
              
              if(system->evflag)
                system->eng_vdwl += e_hb;
            }//if
          }//for-itr
        }//for-pk
      }//if
    }//for-j
  }//for-jst
  #endif
}

#endif




#ifdef CPE
#include "STUBS/mpi.h"
#include <dma.h>
#include "slave.h"

#include "poly_math.h"
#define DMA_FAST
#include "dma_macros.h"

//#define LWPF_UNIT U(HBOND)
//#define LWPF_KERNELS K(ALL) K(BEFORE)  K(CMP) K(BOND) K(IF) K(PK)
//#include "lwpf2.h"

//read cache;
#define READ_C_H    9
#define READ_C_S    4
#define READ_C_LSZ  (1 << READ_C_S)
#define READ_C_LCNT (1 << (READ_C_H - READ_C_S))
#define READ_C_MM   (READ_C_LSZ - 1)
#define READ_C_LM   (READ_C_LCNT - 1)
//index, end_index for atom k;
#define K_C_H    9
#define K_C_S    4
#define K_C_LSZ  (1 << K_C_S)
#define K_C_LCNT (1 << (K_C_H - K_C_S))
#define K_C_MM   (K_C_LSZ - 1)
#define K_C_LM   (K_C_LCNT - 1)

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
void read_cache_index(int k, 
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

#define DT_CACHE_F fCdDelta,rvec4,5,3,6,FADDD
#define EVFLAG 1
#include "reaxc_hydrogen_bonds_cpe.h"
#undef EVFLAG

#define EVFLAG 0
#include "reaxc_hydrogen_bonds_cpe.h"
#undef EVFLAG

#endif
