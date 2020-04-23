#include <math.h>
#include "reaxc_defs_sunway.h"
#include "reaxc_ctypes_sunway.h"
#include "reaxc_inlines_sw64.h"
#include "reaxc_bonds_sw64.h"
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
//#define LWPF_UNITS U(BONDS)
//#include "lwpf2.h"

extern SLAVE_FUN(Bonds_C_evflag_cpe)(bonds_pack_t *);
//swcache_lock_t *locks_frc = NULL;

void Bonds_C(bonds_pack_t *param) 
{
  return;

  reax_system_c *system         = param->system;
  control_params *control       = param->control;
  simulation_data *data         = param->data;
  storage *workspace            = param->workspace;
  reax_list **lists             = param->lists;
  output_controls *out_control  = param->out_control;

  int i, j, pj, natoms;
  int start_i, end_i;
  int type_i, type_j;
  double ebond, pow_BOs_be2, exp_be12, CEbo;
  double gp3, gp4, gp7, gp10, gp37;
  double exphu, exphua1, exphub1, exphuov, hulpov, estriph;
  double decobdbo, decobdboua, decobdboub;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  bond_order_data *bo_ij;
  reax_list *bonds;

  bonds = (*lists) + BONDS;
    
  
  gp3 = system->reax_param.gp.l[3];
  gp4 = system->reax_param.gp.l[4];
  gp7 = system->reax_param.gp.l[7];
  gp10 = system->reax_param.gp.l[10];
  gp37 = (int) system->reax_param.gp.l[37];
  natoms = system->n;

  //#define SW_BONDS
  #ifdef  SW_BONDS
  bonds_pack_t pm;
  pm.system       = system;
  pm.control      = control;
  pm.data         = data;
  pm.workspace    = workspace;
  pm.lists        = lists;
  double  packed_eng[4];
  packed_eng[0] = packed_eng[1] = packed_eng[2] = packed_eng[3] = 0;
  pm.packed_eng     = packed_eng;
  pm.gp3            = gp3;
  pm.gp4            = gp4;
  pm.gp7            = gp7;
  pm.gp10           = gp10;
  pm.gp37           = gp37;
  pm.n              = system->n;
  pm.evflag         = system->evflag;
 
  pm.bond_list      = bonds->select.bond_list;
  pm.bo_data_list   = bonds->bo_data_list;
  pm.Cdbo_list      = bonds->Cdbo_list;
  pm.Cdbopi_list    = bonds->Cdbopi_list;
  pm.Cdbopi2_list   = bonds->Cdbopi2_list;
  pm.BO_list        = bonds->BO_list;
  pm.BOpi_list      = bonds->BOpi_list;
  pm.index          = bonds->index;
  pm.end_index      = bonds->end_index;
  pm.packed_atoms   = system->packed_atoms;
  pm.bo_dboc        = workspace->bo_dboc;
  pm.Delta          = workspace->Delta;
  pm.fCdDelta       = workspace->fCdDelta;
  pm.ntypes         = system->reax_param.num_atom_types;
  
  int nall = bonds->num_intrs;
  //if(locks_frc == NULL)
  //  locks_frc = swcache_u_prepare_locks(nall, 3);
  //pm.locks_frc      = locks_frc;
  pm.locks_frc      = bonds->locks_frc;

  int s1, s2;
  for(s1 = 0; s1 < pm.ntypes; s1++)
  {
    pm.sbp[s1] = system->reax_param.sbp[s1].mass;
    for(s2 = 0; s2 < pm.ntypes; s2++)
    {
      pm.tbp[s1*pm.ntypes+s2].p_be1  = system->reax_param.tbp[s1][s2].p_be1; 
      pm.tbp[s1*pm.ntypes+s2].p_be2  = system->reax_param.tbp[s1][s2].p_be2; 
      pm.tbp[s1*pm.ntypes+s2].De_s   = system->reax_param.tbp[s1][s2].De_s; 
      pm.tbp[s1*pm.ntypes+s2].De_p   = system->reax_param.tbp[s1][s2].De_p; 
      pm.tbp[s1*pm.ntypes+s2].De_pp  = system->reax_param.tbp[s1][s2].De_pp; 
    }
  }

  //printf("in CPE\n");
  //perf_config_t conf;
  //conf.pcrc = PCRC_ALL;
  //conf.pcr0 = PC0_CYCLE;
  //conf.pcr1 = PC1_CYCLE;
  ////conf.pcr2 = PC2_N_GLD;
  ////conf.pcr2 = PC2_N_GF_AND_A;
  //conf.pcr2 = PC2_N_DMA_REQ;
  //lwpf_init(&conf);
  //int myrank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(athread_idle() == 0)
    athread_init();
  athread_spawn(Bonds_C_evflag_cpe, &pm);
  athread_join();

  system->eng_vdwl += packed_eng[0];
  data->my_en.e_bond += packed_eng[1];

  //if(myrank == 0)
  //{
  //  lwpf_report_summary(stdout, &conf);
  //}

  //if(locks_frc != NULL)
  //{
  //  free(locks_frc);
  //  locks_frc = NULL;
  //}
  return;
  
  #else
  double *Cdbo_list     = bonds->Cdbo_list;
  double *Cdbopi_list   = bonds->Cdbopi_list;
  double *Cdbopi2_list  = bonds->Cdbopi2_list;
  double *BO_list       = bonds->BO_list;
  rvec2 *BOpi_list      = bonds->BOpi_list;

  for( i = 0; i < natoms; ++i ) 
  {
    start_i = Start_Index(i, bonds);
    end_i = End_Index(i, bonds);

    for( pj = start_i; pj < end_i; ++pj ) 
    {
      j = bonds->select.bond_list[pj].nbr;

      if( system->packed_atoms[i].orig_id > system->packed_atoms[j].orig_id )
	      continue;
      if( system->packed_atoms[i].orig_id == system->packed_atoms[j].orig_id ) 
      {
        if (system->packed_atoms[j].x[2] <  system->packed_atoms[i].x[2]) continue;
      	if (system->packed_atoms[j].x[2] == system->packed_atoms[i].x[2] &&
      	    system->packed_atoms[j].x[1] <  system->packed_atoms[i].x[1]) continue;
        if (system->packed_atoms[j].x[2] == system->packed_atoms[i].x[2] &&
      	    system->packed_atoms[j].x[1] == system->packed_atoms[i].x[1] &&
      	    system->packed_atoms[j].x[0] <  system->packed_atoms[i].x[0]) continue;
      }

      /* set the pointers */
      type_i = system->packed_atoms[i].type;
      type_j = system->packed_atoms[j].type;
      sbp_i = &( system->reax_param.sbp[type_i] );
      sbp_j = &( system->reax_param.sbp[type_j] );
      twbp = &( system->reax_param.tbp[type_i][type_j] );
      //bo_ij = &( bonds->select.bond_list[pj].bo_data );
      bo_ij = &( bonds->bo_data_list[pj]);

      /* calculate the constants */
      if (bo_ij->BO_s == 0.0) pow_BOs_be2 = 0.0;
      else pow_BOs_be2 = pow( bo_ij->BO_s, twbp->p_be2 );
      exp_be12 = exp( twbp->p_be1 * ( 1.0 - pow_BOs_be2 ) );
      CEbo = -twbp->De_s * exp_be12 *
	        ( 1.0 - twbp->p_be1 * twbp->p_be2 * pow_BOs_be2 );

      /* calculate the Bond Energy */
      data->my_en.e_bond += ebond =
	      -twbp->De_s * bo_ij->BO_s * exp_be12
	      -twbp->De_p * BOpi_list[pj][0]
	      -twbp->De_pp * BOpi_list[pj][1];
        //-twbp->De_p * bo_ij->BO_pi
	      //-twbp->De_pp * bo_ij->BO_pi2;


      /* tally into per-atom energy */
      if( system->evflag)
      {
        system->eng_vdwl += ebond;
      }

      /* calculate derivatives of Bond Orders */
      //bo_ij->Cdbo += CEbo;
      Cdbo_list[pj] += CEbo;

      //bo_ij->Cdbopi -= (CEbo + twbp->De_p);
      Cdbopi_list[pj] -= (CEbo + twbp->De_p);

      //bo_ij->Cdbopi2 -= (CEbo + twbp->De_pp);
      Cdbopi2_list[pj] -= (CEbo + twbp->De_pp);

      /* Stabilisation terminal triple bond */
      //if( bo_ij->BO >= 1.00 ) 
      if(BO_list[pj] >= 1.00 ) 
      {
	      if( gp37 == 2 ||
	      (sbp_i->mass == 12.0000 && sbp_j->mass == 15.9990) ||
	      (sbp_j->mass == 12.0000 && sbp_i->mass == 15.9990) ) 
        {
	        //exphu = exp( -gp7 * SQR(bo_ij->BO - 2.50) );
          //exphua1 = exp(-gp3 * (workspace->bo_dboc[i][0]-bo_ij->BO));
	        //exphub1 = exp(-gp3 * (workspace->bo_dboc[j][0]-bo_ij->BO));

	        exphu = exp( -gp7 * SQR(BO_list[pj] - 2.50) );
	        exphua1 = exp(-gp3 * (workspace->bo_dboc[i][0]-BO_list[pj]));
	        exphub1 = exp(-gp3 * (workspace->bo_dboc[j][0]-BO_list[pj]));
	        exphuov = exp(gp4 * (workspace->Delta[i] + workspace->Delta[j]));
	        hulpov = 1.0 / (1.0 + 25.0 * exphuov);

	        estriph = gp10 * exphu * hulpov * (exphua1 + exphub1);
	        data->my_en.e_bond += estriph;

	        decobdbo = gp10 * exphu * hulpov * (exphua1 + exphub1) *
	          ( gp3 - 2.0 * gp7 * (BO_list[pj]-2.50) );
	          //( gp3 - 2.0 * gp7 * (bo_ij->BO-2.50) );

	        decobdboua = -gp10 * exphu * hulpov *
	          (gp3*exphua1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));
	        decobdboub = -gp10 * exphu * hulpov *
	          (gp3*exphub1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));

	        /* tally into per-atom energy */
	        if(system->evflag)
          {
            system->eng_vdwl += estriph;
          }

	        //bo_ij->Cdbo += decobdbo;
	        Cdbo_list[pj] += decobdbo;
	        workspace->fCdDelta[i][3] += decobdboua;
	        workspace->fCdDelta[j][3] += decobdboub;
	      }//if
      }//if
    }//for-pj
  }//for-i
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
//#define LWPF_UNIT U(BONDS)
//#define LWPF_KERNELS K(ALL) K(FJ) K(FI) K(FL)
//#include "lwpf2.h"

//read cache;
#define READ_C_H    8
#define READ_C_S    3
#define READ_C_LSZ  (1 << READ_C_S)
#define READ_C_LCNT (1 << (READ_C_H - READ_C_S))
#define READ_C_MM   (READ_C_LSZ - 1)
#define READ_C_LM   (READ_C_LCNT - 1)

void read_cache(int i,
                atom_pack_t atoms_cache[][READ_C_LSZ], 
                rvec2 bo_dboc_cache[][READ_C_LSZ],
                double Delta_cache[][READ_C_LSZ],
                int *read_ctag,
                atom_pack_t *packed_atoms,
                rvec2 *bo_dboc,
                double *Delta,
                atom_pack_t *atom_i,
                rvec2 bo_dboc_i,
                double *Delta_i)
{
  dma_init();
  if (read_ctag[(i >> READ_C_S) & READ_C_LM] != i >> READ_C_S)
  {
    pe_get(packed_atoms + (i & ~READ_C_MM), atoms_cache[(i >> READ_C_S) & READ_C_LM], sizeof(atom_pack_t) * READ_C_LSZ);
    pe_get(bo_dboc + (i & ~READ_C_MM), bo_dboc_cache[(i >> READ_C_S) & READ_C_LM], sizeof(rvec2) * READ_C_LSZ);
    pe_get(Delta + (i & ~READ_C_MM), Delta_cache[(i >> READ_C_S) & READ_C_LM], sizeof(double) * READ_C_LSZ);
    dma_syn();
    read_ctag[(i >> READ_C_S) & READ_C_LM] = i >> READ_C_S;
  }
  *atom_i = atoms_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
  *Delta_i = Delta_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM];
  (bo_dboc_i)[0] = bo_dboc_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM][0];
  (bo_dboc_i)[1] = bo_dboc_cache[(i >> READ_C_S) & READ_C_LM][i & READ_C_MM][1];
}

#define DT_CACHE_F fCdDelta,rvec4,5,3,6,FADDD
void Bonds_C_evflag_cpe(bonds_pack_t *param) 
{
  //lwpf_enter(BONDS);
  //lwpf_start(ALL);

  dma_init();
  bonds_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(bonds_pack_t));
  dma_syn();

  reax_system_c *system         = l_pm.system;
  control_params *control       = l_pm.control;
  simulation_data *data         = l_pm.data;
  storage *workspace            = l_pm.workspace;
  reax_list **lists             = l_pm.lists;
  output_controls *out_control  = l_pm.out_control;

  int i, j, pj, pj_off, natoms;
  int start_i, end_i, len_i;
  int type_i, type_j;
  double ebond, pow_BOs_be2, exp_be12, CEbo;
  double gp3, gp4, gp7, gp10, gp37;
  double exphu, exphua1, exphub1, exphuov, hulpov, estriph;
  double decobdbo, decobdboua, decobdboub;
  //single_body_parameters *sbp_i, *sbp_j;
  //two_body_parameters *twbp;
  double sbp_i_mass, sbp_j_mass;
  bonds_tbp_t *twbp;

  bond_order_data *bo_ij;
  reax_list *bonds;

  //bonds = (*lists) + BONDS;
  double *Cdbo_list     = l_pm.Cdbo_list;
  double *Cdbopi_list   = l_pm.Cdbopi_list;
  double *Cdbopi2_list  = l_pm.Cdbopi2_list;
  double *BO_list       = l_pm.BO_list;
  rvec2 *BOpi_list      = l_pm.BOpi_list;
  rvec4 *fCdDelta       = l_pm.fCdDelta;

  gp3     = l_pm.gp3; //system->reax_param.gp.l[3];
  gp4     = l_pm.gp4; //system->reax_param.gp.l[4];
  gp7     = l_pm.gp7; //system->reax_param.gp.l[7];
  gp10    = l_pm.gp10;//system->reax_param.gp.l[10];
  gp37    = l_pm.gp37;//(int) system->reax_param.gp.l[37];
  natoms  = l_pm.n;//system->n;

  int index[ISTEP], end_index[ISTEP];
  atom_pack_t atom_i[ISTEP];
  double Delta_i[ISTEP];
  rvec2 bo_dboc_i[ISTEP];

  bond_data blist_j[MAX_LEN];
  bond_order_data bo_data_j[MAX_LEN];
  double BO_list_j[MAX_LEN];
  rvec2 BOpi_list_j[MAX_LEN];
  double Cdbo_list_j[MAX_LEN], Cdbopi_list_j[MAX_LEN], Cdbopi2_list_j[MAX_LEN];

  //read_cache;
  int read_ctag[READ_C_LCNT];
  atom_pack_t atoms_cache[READ_C_LCNT][READ_C_LSZ];
  atom_pack_t atom_j;
  double Delta_cache[READ_C_LCNT][READ_C_LSZ];
  rvec2 bo_dboc_cache[READ_C_LCNT][READ_C_LSZ];
  double Delta_j;
  rvec2 bo_dboc_j;

  for(i = 0; i < READ_C_LCNT; i++)
  {
    read_ctag[i] = -1;
  }

  swcache_lock_t *locks_frc = l_pm.locks_frc;
  SWCACHE_INIT_U(DT_CACHE_F, locks_frc);

  doublev4 eng_virial[2];
  eng_virial[0] = eng_virial[1] = 0;
  double *eng_vdwl = (double*)(void*)eng_virial;
  double *eng_ebond = eng_vdwl+1;

  int ist, ied, isz, ioff;

  //for( i = 0; i < natoms; ++i ) 
  for(ist = _MYID * ISTEP; ist < natoms; ist+=64*ISTEP) 
  {
    ied = ist + ISTEP;
    if(ied > natoms)
      ied = natoms;
    isz = ied - ist;

    pe_get(l_pm.index+ist,        index,      sizeof(int)*isz);
    pe_get(l_pm.end_index+ist,    end_index,  sizeof(int)*isz);
    pe_get(l_pm.packed_atoms+ist, atom_i,     sizeof(atom_pack_t)*isz);
    pe_get(l_pm.Delta+ist,        Delta_i,    sizeof(double)*isz);
    pe_get(l_pm.bo_dboc+ist,      bo_dboc_i,  sizeof(rvec2)*isz);
    dma_syn();

    for(i = ist; i < ied; i++)
    {
      ioff = i - ist;
      start_i = index[ioff];//l_pm.index[i];//Start_Index(i, bonds);
      end_i   = end_index[ioff];//l_pm.end_index[i];//End_Index(i, bonds);
      len_i = end_i - start_i;
      if(len_i <= 0) continue;
      pe_get(l_pm.bond_list+start_i,    blist_j,        sizeof(bond_data)*len_i);
      pe_get(l_pm.bo_data_list+start_i, bo_data_j,      sizeof(bond_order_data)*len_i);
      pe_get(l_pm.BO_list+start_i,      BO_list_j,      sizeof(double)*len_i);
      pe_get(l_pm.BOpi_list+start_i,    BOpi_list_j,    sizeof(rvec2)*len_i);
      pe_get(l_pm.Cdbo_list+start_i,    Cdbo_list_j,    sizeof(double)*len_i);
      pe_get(l_pm.Cdbopi_list+start_i,  Cdbopi_list_j,  sizeof(double)*len_i);
      pe_get(l_pm.Cdbopi2_list+start_i, Cdbopi2_list_j, sizeof(double)*len_i);
      dma_syn();

      rvec4 fci_tmp;
      fci_tmp[0] = fci_tmp[1] = fci_tmp[2] = fci_tmp[3] = 0.0;

      for( pj = start_i; pj < end_i; ++pj ) 
      {
        pj_off = pj-start_i;
        //j = bonds->select.bond_list[pj].nbr;
        j = blist_j[pj_off].nbr;
     
        read_cache(j, atoms_cache, bo_dboc_cache, Delta_cache,read_ctag,
                  l_pm.packed_atoms, l_pm.bo_dboc, l_pm.Delta,
                  &atom_j, bo_dboc_j, &Delta_j);
        int flag = 1;
        if(atom_i[ioff].orig_id > atom_j.orig_id) flag = 0;
        if(atom_i[ioff].orig_id == atom_j.orig_id) 
        {
          if (atom_j.x[2] <  atom_i[ioff].x[2]) flag = 0;
        	if (atom_j.x[2] == atom_i[ioff].x[2] &&
        	    atom_j.x[1] <  atom_i[ioff].x[1]) flag = 0;
          if (atom_j.x[2] == atom_i[ioff].x[2] &&
        	    atom_j.x[1] == atom_i[ioff].x[1] &&
        	    atom_j.x[0] <  atom_i[ioff].x[0]) flag = 0;
        }
     
        if(flag)
        {
          rvec4 fcj_tmp;
          fcj_tmp[0] = fcj_tmp[1] = fcj_tmp[2] = fcj_tmp[3] = 0.0;

          /* set the pointers */
          type_i = atom_i[ioff].type;
          type_j = atom_j.type;
          sbp_i_mass = l_pm.sbp[type_i];//system->reax_param.sbp[type_i];
          sbp_j_mass = l_pm.sbp[type_j];//system->reax_param.sbp[type_j];
          //twbp = &( system->reax_param.tbp[type_i][type_j] );
          twbp = &(l_pm.tbp[type_i*l_pm.ntypes+type_j] );
          //bo_ij = &( bonds->bo_data_list[pj]);
          bo_ij = &(bo_data_j[pj_off]);

          /* calculate the constants */
          if (bo_ij->BO_s == 0.0) pow_BOs_be2 = 0.0;
          else pow_BOs_be2 = p_powd( bo_ij->BO_s, twbp->p_be2 );
          exp_be12 = p_expd( twbp->p_be1 * ( 1.0 - pow_BOs_be2 ) );
          CEbo = -twbp->De_s * exp_be12 *
	            ( 1.0 - twbp->p_be1 * twbp->p_be2 * pow_BOs_be2 );

          /* calculate the Bond Energy */
          //data->my_en.e_bond += ebond =
          *eng_ebond += ebond =
	          -twbp->De_s * bo_ij->BO_s * exp_be12
	          -twbp->De_p * BOpi_list_j[pj_off][0]
	          -twbp->De_pp * BOpi_list_j[pj_off][1];
            //-twbp->De_p * BOpi_list[pj][0]
	          //-twbp->De_pp * BOpi_list[pj][1];


          /* tally into per-atom energy */
          if(l_pm.evflag)
          {
            //system->eng_vdwl += ebond;
            *eng_vdwl += ebond;
          }

          /* calculate derivatives of Bond Orders */
          Cdbo_list_j[pj_off]     += CEbo;
          Cdbopi_list_j[pj_off]   -= (CEbo + twbp->De_p);
          Cdbopi2_list_j[pj_off]  -= (CEbo + twbp->De_pp);

          /* Stabilisation terminal triple bond */
          //if(BO_list[pj] >= 1.00 ) 
          if(BO_list_j[pj_off] >= 1.00 ) 
          {
	          if( gp37 == 2 ||
	          (sbp_i_mass == 12.0000 && sbp_j_mass == 15.9990) ||
	          (sbp_j_mass == 12.0000 && sbp_i_mass == 15.9990) ) 
            {
	            //exphu   = p_expd( -gp7 * SQR(BO_list[pj] - 2.50) );
	            exphu   = p_expd( -gp7 * SQR(BO_list_j[pj_off] - 2.50) );
	            //exphua1 = p_expd(-gp3 * (workspace->bo_dboc[i][0]-BO_list_j[pj_off]));
	            exphua1 = p_expd(-gp3 * (bo_dboc_i[ioff][0]-BO_list_j[pj_off]));
	            exphub1 = p_expd(-gp3 * (bo_dboc_j[0]-BO_list_j[pj_off]));
	            //exphuov = p_expd(gp4 * (workspace->Delta[i] + workspace->Delta[j]));
	            exphuov = p_expd(gp4 * (Delta_i[ioff] + Delta_j));
	            hulpov = 1.0 / (1.0 + 25.0 * exphuov);

	            estriph = gp10 * exphu * hulpov * (exphua1 + exphub1);
	            //data->my_en.e_bond += estriph;
	            *eng_ebond += estriph;

	            decobdbo = gp10 * exphu * hulpov * (exphua1 + exphub1) *
	                      (gp3 - 2.0 * gp7 * (BO_list_j[pj_off]-2.50));

	            decobdboua = -gp10 * exphu * hulpov *
	                          (gp3*exphua1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));
	            decobdboub = -gp10 * exphu * hulpov *
	                          (gp3*exphub1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));

	            /* tally into per-atom energy */
	            if(l_pm.evflag)
              {
                //system->eng_vdwl += estriph;
                *eng_vdwl += estriph;
              }

	            Cdbo_list_j[pj_off] += decobdbo;
	            //workspace->fCdDelta[i][3] += decobdboua;
              fci_tmp[3] += decobdboua;
              fcj_tmp[3] = decobdboub;
	            //workspace->fCdDelta[j][3] += fcj_tmp[3];
              //lwpf_start(FJ);
              SWCACHE_UPDATE(DT_CACHE_F, j, (&fcj_tmp));
              //lwpf_stop(FJ);
	          }//if
          }//if
        }//if-flag
      }//for-pj
      pe_put(l_pm.Cdbo_list+start_i,    Cdbo_list_j,    sizeof(double)*len_i);
      pe_put(l_pm.Cdbopi_list+start_i,  Cdbopi_list_j,  sizeof(double)*len_i);
      pe_put(l_pm.Cdbopi2_list+start_i, Cdbopi2_list_j, sizeof(double)*len_i);
      dma_syn();
      //workspace->fCdDelta[i][3] += fci_tmp[3];
      //lwpf_start(FI);
      SWCACHE_UPDATE(DT_CACHE_F, i, (&fci_tmp));
      //lwpf_stop(FI);
    }//for-i
  }//for-ist
  
  //lwpf_start(FL);
  SWCACHE_FLUSH(DT_CACHE_F);
  //lwpf_stop(FL);
  reg_reduce_inplace_doublev4(eng_virial, 1);
  if(_MYID == 0)
  {
    pe_put(l_pm.packed_eng, eng_virial, sizeof(double) * 4);
    dma_syn();
  }
  //lwpf_stop(ALL);
  //lwpf_exit(BONDS);

}
#endif
