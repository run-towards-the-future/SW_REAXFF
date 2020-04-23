#ifndef NPAIR_FULL_BIN_GHOST_SW64_H_
#define NPAIR_FULL_BIN_GHOST_SW64_H_
#ifdef __cplusplus
extern "C"{
#endif
  typedef struct nbin_pack_t
  {
    int nbinx,nbiny,nbinz;
    int mbinx,mbiny,mbinz;
    int mbinxlo,mbinylo,mbinzlo;
    double bininvx,bininvy,bininvz;
  }nbin_pack_t;

  typedef struct npair_full_bin_ghost_param_t
  {
    int *binpackhead, *binpacknn;
    bin_pack_atom_t *binpack;
    double *cutneighsq;
    int nstencil, *stencil, *stencilxyz;
    int nlocal, nghost, ntotal, ntypes, mbins, nall;
    int maxchunk;
    int mbinx, mbiny, mbinz;
    double bboxlo[3], bboxhi[3];
    nbin_pack_t nbin_in;
  } npair_full_bin_ghost_param_t;
  void npair_full_bin_ghost_sunway_build_packed(npair_full_bin_ghost_param_t *);
#ifdef __cplusplus
}
#endif
#endif
