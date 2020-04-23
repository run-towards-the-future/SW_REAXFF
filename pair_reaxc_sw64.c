#include "sunway.h"
#include "gptl.h"
#include <math.h>
#include <stdlib.h>
#include <simd.h>
#include "pair_reaxc_sw64.h"

#ifdef MPE
#include <math.h>
#include <athread.h>
extern SLAVE_FUN(write_reax_atoms_and_pack_para)(write_atoms_param_t *);
extern SLAVE_FUN(write_reax_lists_c_para)(write_lists_param_t *);
void write_reax_atoms_and_pack(write_atoms_param_t *pm)
{
  if(athread_idle() == 0)
    athread_init();

  athread_spawn(write_reax_atoms_and_pack_para, pm);
  athread_join();

  //int i;
  //for (i = 0; i < pm->csys->N; i ++)
  //{
  //  pm->csys->my_atoms[i].orig_id    = pm->tag[i];
  //  pm->csys->my_atoms[i].type       = pm->map[pm->type[i]];
  //  pm->csys->my_atoms[i].x[0]       = pm->x[i][0];
  //  pm->csys->my_atoms[i].x[1]       = pm->x[i][1];
  //  pm->csys->my_atoms[i].x[2]       = pm->x[i][2];
  //  pm->csys->my_atoms[i].q          = pm->q[i];
  //  pm->csys->my_atoms[i].num_bonds  = pm->num_bonds[i];
  //  pm->csys->my_atoms[i].num_hbonds = pm->num_hbonds[i];

  //  pm->csys->packed_atoms[i].orig_id    = pm->tag[i];
  //  pm->csys->packed_atoms[i].type       = pm->map[pm->type[i]];
  //  pm->csys->packed_atoms[i].x[0]       = pm->x[i][0];
  //  pm->csys->packed_atoms[i].x[1]       = pm->x[i][1];
  //  pm->csys->packed_atoms[i].x[2]       = pm->x[i][2];
  //  pm->csys->packed_atoms[i].q          = pm->q[i];
  //}
}

void write_reax_lists_c(write_lists_param_t *pm)
{
  if(athread_idle() == 0)
    athread_init();

  athread_spawn(write_reax_lists_c_para, pm);
  athread_join();
  //int numall = pm->numall;
  //int *ilist = pm->ilist;
  //int *numneigh = pm->numneigh;
  //int **firstneigh = pm->firstneigh;
  //reax_list *far_nbrs = pm->far_nbrs;
  //far_neighbor_data_full *far_list = far_nbrs->select.far_nbr_list_full;
  //int maxfar = pm->maxfar;
  //atom_pack_t *patoms = pm->patoms;
  //double nonb_cutsq = pm->nonb_cutsq;
  //int ii;
  //for (ii = 0; ii < numall; ii ++)
  //{
  //  int i = ilist[ii];
  //  int *jlist = firstneigh[i];
  //  int jnum = numneigh[i];
  //  int jj;
  //  far_nbrs->index[i] = i * maxfar;
  //  int knum = 0;
  //  far_neighbor_data_full *klist = far_list + i * maxfar;
  //  for (jj = 0; jj < jnum; jj ++)
  //  {
  //    int j = jlist[jj];
  //    j &= NEIGHMASK;
  //    rvec dvec;
  //    dvec[0] = patoms[j].x[0] - patoms[i].x[0];
  //    dvec[1] = patoms[j].x[1] - patoms[i].x[1];
  //    dvec[2] = patoms[j].x[2] - patoms[i].x[2];
  //    double d_sqr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
  //    if (d_sqr <= nonb_cutsq)
  //    {
  //      far_neighbor_data_full *nbr = klist + knum;
  //      nbr->nbr = j;
  //      nbr->d = sqrt(d_sqr);
  //      nbr->dvec[0] = dvec[0];
  //      nbr->dvec[1] = dvec[1];
  //      nbr->dvec[2] = dvec[2];
  //      nbr->rel_box[0] = 0;
  //      nbr->rel_box[1] = 0;
  //      nbr->rel_box[2] = 0;
  //      nbr->type = patoms[j].type;
  //      nbr->orig_id = patoms[j].orig_id;
  //      nbr->q = patoms[j].q;
  //      knum ++;
  //    }
  //  }
  //  far_nbrs->end_index[i] = i * maxfar + knum;
  //}
}
#endif


#ifdef CPE
#include <dma.h>
#define ISTEP 32
#define SNSTEP 128 //128(err on sunway_big)
#define JCACHE_HBIT 10
#define JCACHE_SBIT 3
#define JCACHE_LINESIZE (1 << JCACHE_SBIT)
#define JCACHE_LINECNT (1 << (JCACHE_HBIT - JCACHE_SBIT))
#define JCACHE_MMASK (JCACHE_LINESIZE - 1)
#define JCACHE_LMASK (JCACHE_LINECNT - 1)
#include "poly_math.h"
#define DMA_FAST
#include "dma_macros.h"

void write_reax_atoms_and_pack_para(write_atoms_param_t *pm)
{  
  dma_init();
  write_atoms_param_t l_pm;
  pe_get(pm,&l_pm,sizeof(write_atoms_param_t));
  dma_syn();
  
  reax_system_c *csys = l_pm.csys;
  
  double x[ISTEP][3];
  int type[ISTEP];
  double q[ISTEP];
  rc_tagint tag[ISTEP];
  int num_bonds[ISTEP];
  int num_hbonds[ISTEP];
  int N = csys->N;
  int ntypes = l_pm.ntypes;
  int map[ntypes + 1];
  
  //reax_atom_part my_atoms;
  atom_pack_t    packed_atoms[ISTEP];
  int my_bonds[2];
  
  int ist,ied,isz;
  int i;
  
  pe_get(l_pm.map,map,(ntypes + 1) * sizeof(int));
  
  for(ist = _MYID * ISTEP; ist < N; ist += ISTEP * 64){

    ied = ist + ISTEP;
    if(ied > N)
      ied = N;
    isz = ied - ist;

    pe_get(l_pm.tag + ist,tag,isz * sizeof(rc_tagint));
    pe_get(l_pm.type + ist,type,isz * sizeof(int));
    pe_get(l_pm.x[ist],x,isz * sizeof(double) * 3);
    pe_get(l_pm.q + ist,q,isz * sizeof(double));
    pe_get(l_pm.num_bonds + ist,num_bonds,isz * sizeof(int));
    pe_get(l_pm.num_hbonds + ist,num_hbonds,isz * sizeof(int));
    dma_syn();

    for (i = ist; i < ied; i++)
    {
      int ii = i - ist; 
      
      //my_atoms.orig_id    = tag[ii];
      //my_atoms.type       = map[type[ii]];
      //my_atoms.x[0]       = x[ii][0];
      //my_atoms.x[1]       = x[ii][1];
      //my_atoms.x[2]       = x[ii][2];
      //my_atoms.q          = q[ii];
      //my_atoms.num_bonds  = num_bonds[ii];
      //my_atoms.num_hbonds = num_hbonds[ii];
      
      //my_bonds[0]  = num_bonds[ii];
      //my_bonds[1] = num_hbonds[ii];

      //pe_put(&(csys->my_atoms[i]),&my_atoms,sizeof(reax_atom_part));
      //pe_put(&(csys->my_atoms[i]),my_bonds,sizeof(int)*2);

      packed_atoms[ii].orig_id    = tag[ii];
      packed_atoms[ii].type       = map[type[ii]];
      packed_atoms[ii].x[0]       = x[ii][0];
      packed_atoms[ii].x[1]       = x[ii][1];
      packed_atoms[ii].x[2]       = x[ii][2];
      packed_atoms[ii].q          = q[ii]; 
      //pe_syn();
    }
    pe_put(&(csys->packed_atoms[ist]),packed_atoms, isz * sizeof(atom_pack_t));
    pe_put(&(csys->num_bonds[ist]),   num_bonds,    isz * sizeof(int));
    pe_put(&(csys->num_hbonds[ist]),  num_hbonds,   isz * sizeof(int));
    dma_syn();

  }
}


void write_reax_lists_c_para(write_lists_param_t *pm)
{
  dma_init();
  write_lists_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(write_lists_param_t));
  dma_syn();

  int numall          = l_pm.numall;
  int maxfar          = l_pm.maxfar;
  double nonb_cutsq   = l_pm.nonb_cutsq;
  double (*x)[3]      = l_pm.x;
  int *ilist          = l_pm.ilist;
  int *numneigh       = l_pm.numneigh;
  int **firstneigh    = l_pm.firstneigh;
  reax_list *far_nbrs = l_pm.far_nbrs;
  atom_pack_t *patoms = l_pm.patoms;
  //far_neighbor_data_full *far_list = far_nbrs->select.far_nbr_list_full;
  far_neighbor_data_full *far_list = l_pm.far_nbr_list_full;

  int ii, ioff;
  int ist, ied, isz;
  int *fn[ISTEP], nn[ISTEP], end_idx[ISTEP], idx[ISTEP];
  int jlist_buf[SNSTEP];
  far_neighbor_data_full fnk[SNSTEP];
  double xi[ISTEP][3];
  
  atom_pack_t j_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_tag[JCACHE_LINECNT];
  atom_pack_t atom_j;

  double *xjc;
  volatile dma_desc get_desc = 0;
  volatile get_reply = 0;
  dma_set_mode(&get_desc, PE_MODE);
  dma_set_size(&get_desc, sizeof(atom_pack_t)*JCACHE_LINESIZE);
  dma_set_op(&get_desc, DMA_GET);
  dma_set_reply(&get_desc, &get_reply);

  int i;
  for (i = 0; i < JCACHE_LINECNT; i ++)
    jcache_tag[i] = -1;

  double stmp;
  for(ist = _MYID * ISTEP; ist < numall; ist += ISTEP * 64)
  {
    ied = ist + ISTEP;
    if(ied > numall)
      ied = numall;
    isz = ied - ist;
    
    pe_get(x[ist],   xi, sizeof(double) * isz *3);
    pe_get(firstneigh + ist,fn, sizeof(int*) * isz);
    pe_get(numneigh + ist,  nn, sizeof(int) * isz);
    dma_syn();
     
    for (ii = ist; ii < ied; ii ++)
    {
      int i = ii;
      ioff = i - ist;
      int *jlist = fn[ioff];
      int jnum = nn[ioff];
      int jj;
      double xtmp, ytmp, ztmp;
      idx[ioff] = i * maxfar;      
      
      xtmp = xi[ioff][0];
      ytmp = xi[ioff][1];
      ztmp = xi[ioff][2];
      int knum = 0;
      far_neighbor_data_full *nbr = far_list + i * maxfar; 
      int jst, jed, jsz;
      for(jst = 0; jst < jnum; jst += SNSTEP)
      {
        jsz = SNSTEP;
        if(jsz + jst > jnum)
          jsz = jnum - jst;
        pe_get(jlist + jst, jlist_buf, sizeof(int)*jsz);
        dma_syn();
        int cnt_fnk = 0;
        for (jj = 0; jj < jsz; jj ++)
        {
          int j = jlist_buf[jj];
          j &= NEIGHMASK;
          int line = j >> JCACHE_SBIT & JCACHE_LMASK;
          int tag = j >> JCACHE_HBIT;
          if(jcache_tag[line] != tag)
          {
            int mem = j & ~JCACHE_MMASK;
            get_reply = 0;
            dma(get_desc, l_pm.patoms + mem, j_cache[line]);
            while(get_reply != 1);
            jcache_tag[line] = tag;
          }

          atom_j = j_cache[line][j & JCACHE_MMASK];
          rvec dvec;
          //dvec[0] = j_cache[line][j & JCACHE_MMASK].x[0] - xtmp;
          //dvec[1] = j_cache[line][j & JCACHE_MMASK].x[1] - ytmp;
          //dvec[2] = j_cache[line][j & JCACHE_MMASK].x[2] - ztmp;
          
          dvec[0] = atom_j.x[0] - xtmp;
          dvec[1] = atom_j.x[1] - ytmp;
          dvec[2] = atom_j.x[2] - ztmp;

          
          double d_sqr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
          if (d_sqr <= nonb_cutsq)
          {
            fnk[cnt_fnk].nbr = j;
            inv_sqrt(d_sqr, stmp);
            fnk[cnt_fnk].d = d_sqr * stmp;
            //fnk[cnt_fnk].d = sqrt(d_sqr);
            fnk[cnt_fnk].dvec[0] = dvec[0];
            fnk[cnt_fnk].dvec[1] = dvec[1];
            fnk[cnt_fnk].dvec[2] = dvec[2];
            fnk[cnt_fnk].rel_box[0] = 0;
            fnk[cnt_fnk].rel_box[1] = 0;
            fnk[cnt_fnk].rel_box[2] = 0;
            //fnk[cnt_fnk].type = j_cache[line][j & JCACHE_MMASK].type;
            //fnk[cnt_fnk].orig_id = j_cache[line][j & JCACHE_MMASK].orig_id;
            //fnk[cnt_fnk].q = j_cache[line][j & JCACHE_MMASK].q;
            
            fnk[cnt_fnk].type     = atom_j.type;
            fnk[cnt_fnk].orig_id  = atom_j.orig_id;
            fnk[cnt_fnk].q        = atom_j.q;

            cnt_fnk++;
          }//if
        }//for-jj
        if(!cnt_fnk) continue;
        //if(cnt_fnk > 0)
        //{
          pe_put(nbr + knum, fnk, sizeof(far_neighbor_data_full) * cnt_fnk);
          dma_syn();
        //}
        knum += cnt_fnk;
      }//for-jst

      end_idx[ioff] = i * maxfar + knum;
    }//for-ii

    //pe_put(far_nbrs->end_index+ist, end_idx, sizeof(int) * isz);
    //pe_put(far_nbrs->index+ist, idx, sizeof(int) * isz);
    pe_put(l_pm.end_index+ist, end_idx, sizeof(int) * isz);
    pe_put(l_pm.index+ist, idx, sizeof(int) * isz);

    dma_syn();
  }//for-ist
}
#endif
