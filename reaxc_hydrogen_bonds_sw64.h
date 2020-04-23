#ifndef REAXC_HYDROGEN_BONDS_SW64_H
#define REAXC_HYDROGEN_BONDS_SW64_H
#include "reaxc_ctypes_sunway.h"
#include "reaxc_nsdef_sunway.h"
#include "swcache.h"
#ifdef __cplusplus
extern "C"{
#else
#endif

/////hydrogen bonds/////
  
  #define NT 4
  #define NT3 64

  typedef struct hydrogen_param_pack_t
  {
    NSDEF(reax_system_c) *system;
    NSDEF(control_params) *control;
    NSDEF(simulation_data) *data;
    NSDEF(storage) *workspace;
    NSDEF(reax_list) **lists;
    NSDEF(output_controls) *out_control;
    
    NSDEF(rvec4) *fCdDelta;
    NSDEF(bond_data) *bond_list;
    NSDEF(hbond_data) *hbond_list;
    NSDEF(atom_pack_t) *packed_atoms;

    double *Cdbo_list, *BO_list;
    int *bindex, *end_bindex;
    int *hbindex, *end_hbindex;
    int *Hindex;
    int n, ntypes;
    double *packed_eng;
   
    int sbp[NT];
    NSDEF(hbond_parameters) hbp[NT3];

    NSDEF(swcache_lock_t) *locks_frc;//fCdDelta;
  }hb_param_pack_t;
#ifdef __cplusplus
}
#endif
#endif
