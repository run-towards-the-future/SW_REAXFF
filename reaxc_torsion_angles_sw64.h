#ifndef REAXC_TORSION_ANGLES_SW64_H
#define REAXC_TORSION_ANGLES_SW64_H
#include "reaxc_ctypes_sunway.h"
#include "reaxc_nsdef_sunway.h"
#include "swcache.h"
#ifdef __cplusplus
extern "C"{
#else
#endif

/////torsion+valence angles/////
  
  #define NT 4
  #define NT3 64
  typedef struct merge_sbp
  {
    double p_val3, p_val5;
  }merge_sbp_t;

  typedef struct merge_param_pack_t
  {
    NSDEF(reax_system_c) *system;
    NSDEF(control_params) *control;
    NSDEF(simulation_data) *data;
    NSDEF(storage) *workspace;
    NSDEF(reax_list) **lists;
    NSDEF(output_controls) *out_control;
    
    NSDEF(rvec2) *BOpi_list, *bo_dboc;
    NSDEF(rvec4) *fCdDelta;
    NSDEF(bond_data) *bond_list;
    NSDEF(atom_pack_t) *packed_atoms;
    NSDEF(three_body_parameters)  *thbp_prm;
    NSDEF(four_body_parameters)   *fbp_prm;

    double *Cdbo_list, *Cdbopi_list, *Cdbopi2_list, *BO_list;
    double *vlpex, *nlp, *dDelta_lp, *Delta, *Delta_val;

    int *fbp_cnt, fbp_tot, *thbp_cnt, thbp_tot;
    int *index, *end_index;
    int n, ntypes;
    double *packed_eng, *gpl;
    double thb_cut, thb_cutsq;
    double p_tor2, p_tor3, p_tor4, p_cot2;
    double p_val6, p_val8, p_val9, p_val10;
    double p_pen2, p_pen3, p_pen4, p_coa2, p_coa3, p_coa4;
    
    NSDEF(swcache_lock_t) *locks;//Cdbo_list;
    NSDEF(swcache_lock_t) *locks_frc;//fCdDelta;
    merge_sbp_t   sbp[NT];
  }merge_param_pack_t;

  

#ifdef __cplusplus
}
#endif
#endif

