#include "sunway.h"
#include <math.h>
#define BNDL_SIZE 1
#define MIN(A,B) ((A)<(B)?(A):(B))
int coord2bin(double *x, int *_ix, int *_iy, int *_iz, 
              nbin_pack_t *nbin_in_, double *bboxlo, double *bboxhi)
{
  int ix, iy, iz;
  int nbinx,nbiny,nbinz;
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;
  double bininvx,bininvy,bininvz;
  //double *bboxlo,*bboxhi;
  nbin_pack_t nbin_in = *nbin_in_;
  nbinx = nbin_in.nbinx;
  nbiny = nbin_in.nbiny;
  nbinz = nbin_in.nbinz;
  mbinx = nbin_in.mbinx;
  mbiny = nbin_in.mbiny;
  mbinz = nbin_in.mbinz;
  mbinxlo = nbin_in.mbinxlo;
  mbinylo = nbin_in.mbinylo;
  mbinzlo = nbin_in.mbinzlo;
  bininvx = nbin_in.bininvx;
  bininvy = nbin_in.bininvy;
  bininvz = nbin_in.bininvz;

  if (x[0] >= bboxhi[0])
    ix = (int)((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = (int)((x[0]-bboxlo[0])*bininvx);
    ix = MIN(ix,nbinx-1);
  } else
    ix = (int)((x[0]-bboxlo[0])*bininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = (int)((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = (int)((x[1]-bboxlo[1])*bininvy);
    iy = MIN(iy,nbiny-1);
  } else
    iy = (int)((x[1]-bboxlo[1])*bininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = (int)((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = (int)((x[2]-bboxlo[2])*bininvz);
    iz = MIN(iz,nbinz-1);
  } else
    iz = (int)((x[2]-bboxlo[2])*bininvz) - 1;

  ix -= mbinxlo;
  iy -= mbinylo;
  iz -= mbinzlo;
  *_ix = ix;
  *_iy = iy;
  *_iz = iz;
  return iz*mbiny*mbinx + iy*mbinx + ix;
}

#ifdef MPE
extern SLAVE_FUN(npair_full_bin_ghost_sunway_build_packed_para)(reaxc_neigh_param_t *pm);

void npair_full_bin_ghost_sunway_build_packed(reaxc_neigh_param_t *pm)
{
  if (athread_idle() == 0)
    athread_init();
  athread_spawn(npair_full_bin_ghost_sunway_build_packed_para, pm);
  athread_join();  
  
} 
#endif

#ifdef CPE
#define BIN_PAGESIZE 64
#define NEIGH_PAGESIZE 128
#define DMA_FAST
#include "dma_macros.h"

void npair_full_bin_ghost_sunway_build_packed_para(reaxc_neigh_param_t *pm)
{
  dma_init();
  reaxc_neigh_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(reaxc_neigh_param_t));
  dma_syn();

  int nstencil = l_pm.nstencil;
  int mbins = l_pm.mbins;
  int stencil[l_pm.nstencil];
  int stencilxyz[l_pm.nstencil][3];
  double cutneighsq[l_pm.ntypes + 1][l_pm.ntypes + 1];
  bin_pack_atom_t binpacki[BIN_PAGESIZE], binpackk[BIN_PAGESIZE];
  int binpackheadi[BNDL_SIZE + 1], binpackheadk[BNDL_SIZE + 1];
  int firstneigh[BIN_PAGESIZE][NEIGH_PAGESIZE];
  int numneigh[BIN_PAGESIZE], totalneigh[BIN_PAGESIZE];
  pe_get(l_pm.stencil, stencil, nstencil * sizeof(int));
  pe_get(l_pm.stencilxyz, stencilxyz, 3 * nstencil * sizeof(int));
  pe_get(l_pm.cutneighsq, cutneighsq, sizeof(double)*(l_pm.ntypes+1)*(l_pm.ntypes+1));
  dma_syn();
  int ibins;

  for (ibins = BNDL_SIZE * _MYID; ibins < mbins; ibins += BNDL_SIZE * 64)
  {
    int ibine = ibins + BNDL_SIZE;
    if (ibine > mbins)
      ibine = mbins;
    int ibin, ibinpacks, ibinpacke, ibinpackn;
    if (ibine - ibins + 1 > 0) 
    {
      pe_get(l_pm.binpackhead + ibins, binpackheadi, (ibine - ibins + 1) * sizeof(int));
      dma_syn();
    
      ibinpacks = binpackheadi[0];
      ibinpacke = binpackheadi[ibine - ibins];
      ibinpackn = ibinpacke - ibinpacks;
      if (ibinpackn > 0)
      {
        pe_get(l_pm.binpack + ibinpacks, binpacki, sizeof(bin_pack_atom_t) * ibinpackn);
        int i;
        for (i = 0; i < ibinpackn; i ++)
          numneigh[i] = totalneigh[i] = 0;
        dma_syn();
      }
    }//if-ibine
    int k;

    for (k = 0; k < l_pm.nstencil; k ++)
    {
      int kbins = ibins + stencil[k];
      int kbine = ibine + stencil[k];
      if (kbins < 0) kbins = 0;
      if (kbine > mbins) kbine = mbins;
      if (kbins >= kbine) continue;
      int kbinpacks, kbinpacke, kbinpackn;
      if (kbine - kbins + 1 > 0) 
      {
        pe_get(l_pm.binpackhead + kbins, binpackheadk, (kbine-kbins+1) * sizeof(int));
        dma_syn();
        kbinpacks = binpackheadk[0];
        kbinpacke = binpackheadk[kbine - kbins];
        kbinpackn = kbinpacke - kbinpacks;
        if (kbinpackn > 0){
          pe_get(l_pm.binpack + kbinpacks, binpackk, sizeof(bin_pack_atom_t)*kbinpackn);
          dma_syn();
        }
      }//if-kbine

      for (ibin = ibins; ibin < ibine; ibin ++)
      {
        int iboff = ibin - ibins;
        if (ibin + stencil[k] > 0 && ibin + stencil[k] < mbins)
        {
          int js = binpackheadk[ibin - kbins + stencil[k]];
          int je = binpackheadk[ibin - kbins + stencil[k] + 1];
          int j;
          int i;
          for (i = binpackheadi[iboff]; i < binpackheadi[iboff + 1]; i ++)
          {
          
              bin_pack_atom_t *iatom = binpacki + i - ibinpacks;
              if(iatom->id >= l_pm.nlocal)
              {
                int xbin, ybin, zbin;
                int tmp = coord2bin( iatom->x, &xbin, &ybin, &zbin,
                        &l_pm.nbin_in, l_pm.bboxlo, l_pm.bboxhi);

                int xbin2 = xbin + stencilxyz[k][0];
                int ybin2 = ybin + stencilxyz[k][1];
                int zbin2 = zbin + stencilxyz[k][2];
                if (xbin2 < 0 || xbin2 >= l_pm.mbinx ||
                    ybin2 < 0 || ybin2 >= l_pm.mbiny ||
                    zbin2 < 0 || zbin2 >= l_pm.mbinz) continue;
              }
              
            for (j = js; j < je; j ++)
            {
              bin_pack_atom_t *jatom = binpackk + j - kbinpacks;
              int jtype = jatom->type;
            
              if (iatom->id != jatom->id)
              {
                int itype = iatom->type;
                double delx = iatom->x[0] - jatom->x[0];
                double dely = iatom->x[1] - jatom->x[1];
                double delz = iatom->x[2] - jatom->x[2];
                double rsq = delx * delx + dely * dely + delz * delz;
                if (rsq <= cutneighsq[itype][jtype])
                {
                  int ioff = i - ibinpacks;
                  firstneigh[ioff][numneigh[ioff] ++] = jatom->id;
                  if (numneigh[ioff] >= NEIGH_PAGESIZE)
                  {
                    pe_put(iatom->firstneigh + totalneigh[ioff], firstneigh[ioff], 
                           sizeof(int) * NEIGH_PAGESIZE);
                    totalneigh[ioff] += numneigh[ioff];
                    numneigh[ioff] = 0;
                    dma_syn();
                  }
                  //l_pm.binpacknn[i] ++;
                  //iatom->firstneigh[l_pm.binpacknn[i] ++] = jatom->id;
                }//if
              }//if
            }//for-i
          }//for-j
        }//if
      }//for-ibin
    }//for-k

    if (ibinpackn > 0)
    {
      int i;
      for (i = 0; i < ibinpackn; i ++)
      {
        if (numneigh[i] > 0)
        {
          pe_put(binpacki[i].firstneigh + totalneigh[i], firstneigh[i], 
                  sizeof(int) * numneigh[i]);
          totalneigh[i] += numneigh[i];
          dma_syn();
        }
      }//for-i
      pe_put(l_pm.binpacknn + ibinpacks, totalneigh, sizeof(int) * ibinpackn);
      dma_syn();
    }//if
  }//for-ibin
}




#endif
