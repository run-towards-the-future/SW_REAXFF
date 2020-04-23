#ifndef REAXC_BONDS_ORDERS_SW64_H
#define REAXC_BONDS_ORDERS_SW64_H
#include "reaxc_ctypes_sunway.h"
#include "reaxc_nsdef_sunway.h"
#include "swcache.h"
#ifdef __cplusplus
extern "C"{
#else
#endif

  #define NT 4
  #define NT2 16

  typedef struct merge_sbp_t
  {
    double valency, p_ovun2, p_ovun5, mass;
  }merge_sbp_t;

  typedef struct merge_tbp_t
  {
    double p_ovun1;
    double p_be1, p_be2;
    double De_s, De_p, De_pp;
  }merge_tbp_t;

  typedef struct merge_bonds_eng_t
  {
    NSDEF(reax_system_c) *system;
    NSDEF(control_params) *control;
    NSDEF(simulation_data) *data;
    NSDEF(storage) *workspace;
    NSDEF(reax_list) **lists;
    NSDEF(output_controls) *out_control;
    
    NSDEF(rvec2) *BOpi_list, *bo_dboc;
    NSDEF(rvec4) *fCdDelta;
    NSDEF(atom_pack_t) *packed_atoms;
    NSDEF(bond_data) *bond_list;
    NSDEF(bond_order_data) *bo_data_list;
    double *BO_list;
    double *Cdbo_list, *Cdbopi_list, *Cdbopi2_list;
    double *Delta, *dDelta_lp;
    double *Delta_lp_temp;
    double p_lp3;
    double p_ovun3, p_ovun4, p_ovun6, p_ovun7, p_ovun8;
    double gp3, gp4, gp7, gp10, gp37;

    int *index, *end_index;
    int n, ntypes, evflag;
    int enobondsflag;
    double *packed_eng;
    NSDEF(swcache_lock_t) *locks_frc;
   
    merge_sbp_t sbp[NT];
    merge_tbp_t tbp[NT2];
  }merge_bonds_eng_t;
#ifdef __cplusplus
}
#endif
#endif
