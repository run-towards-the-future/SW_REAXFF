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

  typedef struct bo_tbp_t
  {
    double p_boc3, p_boc4, p_boc5;
    double v13cor, ovc;
  }bo_tbp_t;

  typedef struct bo_sbp_t
  {
    double valency, valency_e, valency_boc, valency_val;
    double nlp_opt, mass;
  }bo_sbp_t;

  typedef struct bo_param_pack_t
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
    double *BO_list;
    double *Deltap, *Deltap_boc;

    double *Delta_e, *Delta_val, *Delta;
    double *Delta_lp, *dDelta_lp;
    double *Delta_lp_temp, *dDelta_lp_temp;
    double *nlp, *vlpex, *Clp, *nlp_temp;

    int *index, *end_index;
    int N, ntypes;
    double p_boc1, p_boc2;
   
    double sbp[NT];
    bo_sbp_t spm[NT];
    bo_tbp_t tbp[NT2];
  }bo_param_pack_t;




#ifdef __cplusplus
}
#endif
#endif
