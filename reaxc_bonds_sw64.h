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

  //typedef struct bo_tbp_t
  //{
  //  double p_boc3, p_boc4, p_boc5;
  //  double v13cor, ovc;
  //}bo_tbp_t;

  typedef struct bonds_tbp_t
  {
    /* Bond Energy parameters */
    double p_be1, p_be2;
    double De_s, De_p, De_pp;
  }bonds_tbp_t;

  typedef struct bonds_pack_t
  {
    NSDEF(reax_system_c) *system;
    NSDEF(control_params) *control;
    NSDEF(simulation_data) *data;
    NSDEF(storage) *workspace;
    NSDEF(reax_list) **lists;
    NSDEF(output_controls) *out_control;
    
    NSDEF(atom_pack_t) *packed_atoms;
    NSDEF(rvec2) *BOpi_list, *bo_dboc;
    NSDEF(bond_data) *bond_list;
    NSDEF(bond_order_data) *bo_data_list;
    NSDEF(rvec4)*fCdDelta;
    double *Cdbo_list, *Cdbopi_list, *Cdbopi2_list;
    double *BO_list;
  
    double *Delta;
    int *index, *end_index;
    int n, ntypes, evflag;
    double gp3, gp4, gp7, gp10, gp37;
   
    double sbp[NT];
    bonds_tbp_t tbp[NT2];
    double *packed_eng;
    NSDEF(swcache_lock_t) *locks_frc;
  }bonds_pack_t;
#ifdef __cplusplus
}
#endif
#endif
