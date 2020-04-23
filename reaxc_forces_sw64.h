#ifndef REAXC_BONDS_ORDERS_SW64_H
#define REAXC_BONDS_ORDERS_SW64_H
#include "reaxc_ctypes_sunway.h"
#include "reaxc_nsdef_sunway.h"
#include "swcache.h"
#ifdef __cplusplus
extern "C"{
#else
#endif
  typedef struct add_dbond_pack_t
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
    rvec4 *fCdDelta;
    rvec *dDeltap_self;
    double *Cdbo_list, *Cdbopi_list, *Cdbopi2_list;

    int *index, *end_index;
    int N;
    //swcache_lock_t *locks_dbond;//fCdDelta;
    swcache_lock_t *locks_frc;//fCdDelta;
  }add_dbond_pack_t;

#ifdef __cplusplus
}
#endif
#endif
