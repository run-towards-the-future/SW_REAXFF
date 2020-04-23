#include "sunway.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "reaxc_ctypes_sunway.h"
#include "simd.h"
#define SQR(x) ((x)*(x))
#define EV_TO_KCAL_PER_MOL 14.4
static inline int llf(double dx, double dy, double dz)
{
  double SMALL = 0.0001;
  if (dz < -SMALL) return 1;
  if (dz < SMALL){
    if (dy < -SMALL) return 1;
    if (dy < SMALL && dx < -SMALL) return 1;
  }
  return 0;
}
static inline int urb(double dx, double dy, double dz)
{
  double SMALL = 0.0001;
  if (dz > SMALL) return 1;
  if (dz > -SMALL){
    if (dy > SMALL) return 1;
    if (dy > -SMALL && dx > SMALL) return 1;
  }
  return 0;
}

#ifdef MPE
#include <mpi.h>
#include <athread.h>
//#define LWPF_UNITS U(SPARSE)
//#include "lwpf2.h"

extern SLAVE_FUN(fix_qeq_reax_sunway_a2s)(fix_qeq_pack_t *);
extern SLAVE_FUN(sparse_matvec_C_para)(fix_qeq_pack_t *);
extern SLAVE_FUN(compute_H_Full_C_para1)(fix_qeq_pack_t *);
extern SLAVE_FUN(compute_H_Full_C_para2)(fix_qeq_pack_t *);
extern SLAVE_FUN(compute_H_Full_C_para)(fix_qeq_pack_t *);
int r = 0;
void sparse_matvec_C(fix_qeq_pack_t *param)
{
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
  athread_spawn(sparse_matvec_C_para, param);
  athread_join();
  
  //if(myrank == 0)
  //{
  //  lwpf_report_summary(stdout, &conf);
  //}

}

void sparse_matvec_C_spawn(fix_qeq_pack_t *param) {
  athread_spawn(sparse_matvec_C_para, param);
}

void sparse_matvec_C_join() {
  athread_join();
}

void compute_H_Full_C(fix_qeq_pack_t *param)
{
  if(athread_idle() == 0)
    athread_init();
  //if(r == 0)
  //{
  //  perf_config_t conf;
  //  conf.pcrc = PCRC_ALL;
  //  conf.pcr0 = PC0_CNT_INST;
  //  conf.pcr1 = PC1_CNT_JMP_FAIL;
  //  //conf.pcr2 = PC2_CNT_GLD;
  //  conf.pcr2 = PC2_CNT_GST;
  //  lwpf_init(&conf);
  //}
  //int rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //atom_in_t *atom_in = malloc(sizeof(atom_in_t) * (pm->nlocal + pm->nghost + 32));
  atom_in_t *atom_in = malloc(sizeof(atom_in_t) * (param->list.inum + 32));
  param->atom_in = (void*)((long)atom_in | 255) + 1;
  
  //GPTLstart("a2s");
  athread_spawn(fix_qeq_reax_sunway_a2s, param);
  athread_join();
  //GPTLstop("a2s");
  
  int i, j, ii, jj, flag,jnum;
  int *jlist;
  double dx, dy, dz, r_sqr;
  const double SMALL = 0.0001;
  double swbsq      = param->swbsq;
  int groupbit      = param->groupbit;
  int m_fill        = param->m_fill;
  double **shld     = param->shld;
  int *type         = param->atom.type;
  rc_tagint *tag    = param->atom.tag;
  double (*x)[3]    = param->atom.x;
  int *mask         = param->atom.mask;
  int nlocal        = param->atom.nlocal;
  int inum          = param->list.inum;  
  int *numneigh     = param->list.numneigh;
  int **firstneigh  = param->list.firstneigh;
  int *firstnbr     = param->H.firstnbr;
  int *numnbrs      = param->H.numnbrs;
  int *H_jlist      = param->H.jlist;
  double *val       = param->H.val;
  double *Tap       = param->Tap;
  int maxHlist      = param->maxHlist;

  r_sqr = 0;
  //nlocal = 20736, inum = 65908;
  
  //GPTLstart("para1");
  //athread_spawn(compute_H_Full_C_para1, param);
  //athread_join();
  //GPTLstop("para1");
  //
  //GPTLstart("sum");
  //m_fill = 0;
  ////for (ii = 0; ii < inum; ii ++)
  //for (ii = 0; ii < nlocal; ii ++)
  //{
  //  //firstnbr[ii] = m_fill;
  //  firstnbr[ii] = ii * 512;
  //  //m_fill += numnbrs[ii];
  //}
  //GPTLstop("sum");

  GPTLstart("para");
  //athread_spawn(compute_H_Full_C_para2, param);
  athread_spawn(compute_H_Full_C_para, param);
  athread_join();
  //for(ii = 0; ii < nlocal; ii++) 
  //{
  //  int foff = 0;
  //  i = ii;
  //  numnbrs[i] = 0;
  //  firstnbr[i] = maxHlist * i;
  //  if(mask[i] & groupbit) 
  //  {
  //    jlist = firstneigh[i];
  //    jnum = numneigh[i];
  //    for(jj = 0; jj < jnum; jj++) 
  //    {
  //      j = jlist[jj];
  //      dx = x[j][0] - x[i][0];
  //      dy = x[j][1] - x[i][1];
  //      dz = x[j][2] - x[i][2];
  //      r_sqr = SQR(dx) + SQR(dy) + SQR(dz);

  //      if(r_sqr <= swbsq)
  //      {
  //        numnbrs[i] ++;
  //        H_jlist[firstnbr[i] + foff] = j;
  //        
  //        double Taper, denom;
  //        double r = sqrt(r_sqr);
  //        double gamma = shld[type[i]][type[j]];
  //        Taper = Tap[7] * r + Tap[6];
  //        Taper = Taper * r + Tap[5];
  //        Taper = Taper * r + Tap[4];
  //        Taper = Taper * r + Tap[3];
  //        Taper = Taper * r + Tap[2];
  //        Taper = Taper * r + Tap[1];
  //        Taper = Taper * r + Tap[0];
  //        denom = r * r * r + gamma;
  //        denom = pow(denom,0.3333333333333);

  //        val[firstnbr[i] + foff] = Taper * EV_TO_KCAL_PER_MOL / denom;
  //        foff++;
  //      }//if
  //    }//for-jj
  //    //numnbrs[i] = foff;
  //  }//if
  //}//for-ii

  GPTLstop("para");
 
  //for(ii = 0; ii < inum; ii ++)
  //{
  //  i = ii;
  //  numnbrs[i] = 0;
  //  if (mask[i] & groupbit)
  //  {
  //    jlist = firstneigh[i];
  //    jnum = numneigh[i];
  //    firstnbr[i] = m_fill;
  //    for(jj = 0; jj < jnum; jj++) 
  //    {
  //      j = jlist[jj];
  //      dx = x[j][0] - x[i][0];
  //      dy = x[j][1] - x[i][1];
  //      dz = x[j][2] - x[i][2];
  //      r_sqr = SQR(dx) + SQR(dy) + SQR(dz);

  //      if(r_sqr <= swbsq) 
  //      {
  //        //if(i >= nlocal && j >= nlocal) continue;
  //        //if(i >= nlocal)
  //        //{
  //        //  if(tag[j] > tag[i]) continue;
  //        //  if(tag[j] == tag[i] && urb(dx, dy, dz)) continue;
  //        //}
  //        //if(j >= nlocal)
  //        //{
  //        //  if(tag[i] > tag[j]) continue;
  //        //  if(tag[i] == tag[j] && llf(dx, dy, dz)) continue;
  //        //}
  //        numnbrs[i] ++;
  //      }//if
  //    }//for-jj
  //  }//if
  //}//for-ii

  //m_fill = 0;
  //for (ii = 0; ii < inum; ii ++)
  //{
  //  firstnbr[ii] = m_fill;
  //  m_fill += numnbrs[ii];
  //}
  //
  //for(ii = 0; ii < inum; ii++) 
  //{
  //  int foff = 0;
  //  i = ii;
  //  if(mask[i] & groupbit) 
  //  {
  //    jlist = firstneigh[i];
  //    jnum = numneigh[i];
  //    for(jj = 0; jj < jnum; jj++) 
  //    {
  //      j = jlist[jj];
  //      dx = x[j][0] - x[i][0];
  //      dy = x[j][1] - x[i][1];
  //      dz = x[j][2] - x[i][2];
  //      r_sqr = SQR(dx) + SQR(dy) + SQR(dz);

  //      if(r_sqr <= swbsq)
  //      {
  //        //if(i >= nlocal && j >= nlocal) continue;
  //        //if(i >= nlocal)
  //        //{
  //        //  if(tag[j] > tag[i]) continue;
  //        //  if(tag[j] == tag[i] && urb(dx, dy, dz)) continue;
  //        //}
  //        //if(j >= nlocal)
  //        //{
  //        //  if(tag[i] > tag[j]) continue;
  //        //  if(tag[i] == tag[j] && llf(dx, dy, dz)) continue;
  //        //}
  //        H_jlist[firstnbr[i] + foff] = j;
  //        
  //        double Taper, denom;
  //        double r = sqrt(r_sqr);
  //        double gamma = shld[type[i]][type[j]];
  //        Taper = Tap[7] * r + Tap[6];
  //        Taper = Taper * r + Tap[5];
  //        Taper = Taper * r + Tap[4];
  //        Taper = Taper * r + Tap[3];
  //        Taper = Taper * r + Tap[2];
  //        Taper = Taper * r + Tap[1];
  //        Taper = Taper * r + Tap[0];
  //        denom = r * r * r + gamma;
  //        denom = pow(denom,0.3333333333333);

  //        val[firstnbr[i] + foff] = Taper * EV_TO_KCAL_PER_MOL / denom;
  //        foff++;
  //      }//if
  //    }//for-jj
  //    //numnbrs[i] = foff;
  //  }//if
  //}//for-ii


  //if(r == 10 && rank == 0)
  //{
  //  lwpf_finish(stdout);
  //}
  //r++;
  free(atom_in);
}
#endif


#ifdef CPE
#include "STUBS/mpi.h"
#include <slave.h>
#include <dma.h>
#include "poly_math.h"
#define DMA_FAST
#include "dma_macros.h"

#define ISTEP 64
#define SNSTEP 512
#define JCACHE_HBIT 12
#define JCACHE_SBIT 4
#define JCACHE_LINESIZE (1 << JCACHE_SBIT)
#define JCACHE_LINECNT (1 << (JCACHE_HBIT - JCACHE_SBIT))
#define JCACHE_MMASK (JCACHE_LINESIZE - 1)
#define JCACHE_LMASK (JCACHE_LINECNT - 1)
#define dma_rpl(desc, mem, ldm, reply) \
  asm("dma %0, %1, %2\n\t" : : "r"(desc), "r"(mem), "r"(ldm), "r"(&reply) : "memory");
//#define LWPF_UNIT U(SPARSE)
//#define LWPF_KERNELS K(ALL) K(BEFORE) K(GET) K(JLST) K(JGET) K(ITRJ) K(JBUF) K(JCACHE) K(JCMP)
//#include "lwpf2.h"
void fix_qeq_reax_sunway_a2s(fix_qeq_pack_t *pm)
{
  dma_init();
  fix_qeq_pack_t l_pm;
  pe_get(pm, &l_pm, sizeof(fix_qeq_pack_t));
  dma_syn();

  double x[ISTEP][3];
  rc_tagint tag[ISTEP];
  int type[ISTEP];

  atom_in_t atom_in[ISTEP];
  int ist, ied, isz, i;
  int inum = l_pm.list.inum;
  for(ist = ISTEP * _MYID; ist < inum; ist += ISTEP * 64)
  {
    int ied = ist + ISTEP;
    if(ied > inum)
      ied = inum;
    int isz = ied - ist;
    pe_get(l_pm.atom.x[ist],   x,    sizeof(double) * 3 * isz);
    pe_get(l_pm.atom.type+ist, type, sizeof(int) * isz);
    pe_get(l_pm.atom.tag+ist,  tag,  sizeof(rc_tagint) * isz);
    dma_syn();
    for(i = 0; i < isz; i++)
    {
      atom_in[i].x[0] = x[i][0];
      atom_in[i].x[1] = x[i][1];
      atom_in[i].x[2] = x[i][2];
      atom_in[i].tag  = tag[i];
      atom_in[i].type = type[i];
    }
    pe_put(l_pm.atom_in+ist, atom_in, sizeof(atom_in_t) * isz);
    dma_syn();
  }//for-ist
}

void sparse_matvec_C_para(fix_qeq_pack_t *param)
{
  //lwpf_enter(SPARSE);
  //lwpf_start(ALL);
  //lwpf_start(BEFORE);
  dma_init();
  int i, j, itr_j, ii, ioff;
  fix_qeq_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(fix_qeq_pack_t));
  dma_syn();

  double *x     = l_pm.x;
  double *b     = l_pm.b;
  double *eta   = l_pm.eta;
  int groupbit  = l_pm.groupbit;
  int nn        = l_pm.list.inum;
  int NN        = l_pm.list.nums;
  int *ilist    = l_pm.list.ilist;
  int ntypes    = l_pm.atom.ntypes;
  int *mask     = l_pm.atom.mask;
  int *type     = l_pm.atom.type;
  int *firstnbr = l_pm.H.firstnbr;
  int *numnbrs  = l_pm.H.numnbrs;
  int *jlist    = l_pm.H.jlist;
  double *val   = l_pm.H.val;
 
  int ist, isz, ied;
  int mi[ISTEP], ti[ISTEP], fi[ISTEP], nni[ISTEP];
  double xi[ISTEP];
  int jlist_buf[SNSTEP];
  double val_buf[SNSTEP], bi[ISTEP];
  double l_eta[ntypes+1];
  pe_get(eta,  l_eta,   sizeof(double)*(ntypes+1));
  dma_syn();

  int jst, jed, jsz;
  double j_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_tag[JCACHE_LINECNT];
  for (i = 0; i < JCACHE_LINECNT; i ++)
    jcache_tag[i] = -1;

  volatile dma_desc get_desc = 0;
  volatile get_reply = 0;
  dma_set_mode(&get_desc, PE_MODE);
  dma_set_size(&get_desc, sizeof(double)*JCACHE_LINESIZE);
  dma_set_op(&get_desc, DMA_GET);
  dma_set_reply(&get_desc, &get_reply);
  //lwpf_stop(BEFORE);

  for(ist = _MYID * ISTEP; ist < nn; ist += ISTEP * 64)
  { 
    ied = ist + ISTEP;
    if(ied > nn)
      ied = nn;
    isz = ied - ist;
    //lwpf_start(GET);
    pe_get(mask+ist,    mi,   sizeof(int)*isz);
    pe_get(type+ist,    ti,   sizeof(int)*isz);
    pe_get(x+ist,       xi,   sizeof(double)*isz);
    pe_get(firstnbr+ist,fi,   sizeof(int)*isz);
    pe_get(numnbrs+ist, nni,  sizeof(int)*isz);
    dma_syn();
    //lwpf_stop(GET);

    for (ii = ist; ii < ied; ++ii) 
    {
      i = ii;
      ioff = i - ist;
      if(mi[ioff] & groupbit)
      { 
        bi[ioff] = l_eta[ti[ioff]] * xi[ioff];
        //if(ii < nn)   bi[ioff] = l_eta[ti[ioff]] * xi[ioff];
        //else  bi[ioff] = 0;
      }

      int start_i = fi[ioff];
      int end_i   = fi[ioff]+nni[ioff];
      int len_i = end_i - start_i;
      //lwpf_start(JLST);
      for(jst = start_i; jst < end_i; jst += SNSTEP)
      {
        jsz = SNSTEP;
        if(jst + SNSTEP > end_i)
          jsz = end_i - jst;
        
        //lwpf_start(JGET);
        pe_get(jlist + jst, jlist_buf, sizeof(int)*jsz);
        pe_get(val + jst,   val_buf,   sizeof(double)*jsz);
        //pe_get(jlist+start_i, jlist_buf, sizeof(int)*len_i);
        //pe_get(val+start_i,   val_buf,   sizeof(double)*len_i);
        dma_syn();
        //lwpf_stop(JGET);

        //lwpf_start(ITRJ);
        for(itr_j = 0; itr_j < jsz; ++itr_j) 
        //for(itr_j = 0; itr_j < len_i; ++itr_j) 
        {
          //lwpf_start(JBUF);
          j = jlist_buf[itr_j];
          //if (j < 0 || j > NN) call_printf("%d\n", j);
          int line = j >> JCACHE_SBIT & JCACHE_LMASK;
          int tag = j >> JCACHE_HBIT;
          //lwpf_stop(JBUF);

          //lwpf_start(JCACHE);
          if(jcache_tag[line] != tag)
          {
            int mem = j & ~JCACHE_MMASK;
            get_reply = 0;
            dma_rpl(get_desc, x+mem, j_cache[line], get_reply);
            while(get_reply != 1);
            jcache_tag[line] = tag;
          }
          //lwpf_stop(JCACHE);
          
          //lwpf_start(JCMP);
          double xj = j_cache[line][j & JCACHE_MMASK];
          bi[ioff] += val_buf[itr_j] * xj;
          //lwpf_stop(JCMP);
        }//for-itr_j
        //lwpf_stop(ITRJ);
      }//jst
      //lwpf_stop(JLST);
    }//for-ii

    pe_put(b+ist, bi, sizeof(double)*isz);
    dma_syn();
  }//for-ist
  //lwpf_stop(ALL);
  //lwpf_exit(SPARSE);
}

#define ISTEP 64
#define SNSTEP 512
#define JCACHE_HBIT 10
#define JCACHE_SBIT 3
#define JCACHE_LINESIZE (1 << JCACHE_SBIT)
#define JCACHE_LINECNT (1 << (JCACHE_HBIT - JCACHE_SBIT))
#define JCACHE_MMASK (JCACHE_LINESIZE - 1)
#define JCACHE_LMASK (JCACHE_LINECNT - 1)
#define dma_rpl(desc, mem, ldm, reply) \
  asm("dma %0, %1, %2\n\t" : : "r"(desc), "r"(mem), "r"(ldm), "r"(&reply) : "memory");

void compute_H_Full_C_para1(fix_qeq_pack_t *param)
{
  dma_init();
  fix_qeq_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(fix_qeq_pack_t));
  dma_syn();

  int i, j, ii, jj, flag,jnum;
  int *jlist;
  double dx, dy, dz, r_sqr;
  const double SMALL = 0.0001;
  double swbsq      = l_pm.swbsq;
  int groupbit      = l_pm.groupbit;
  int m_fill        = l_pm.m_fill;
  double **shld     = l_pm.shld;
  int *type         = l_pm.atom.type;
  rc_tagint *tag    = l_pm.atom.tag;
  double (*x)[3]    = l_pm.atom.x;
  int *mask         = l_pm.atom.mask;
  int nlocal        = l_pm.atom.nlocal;
  int inum          = l_pm.list.inum;  
  int *numneigh     = l_pm.list.numneigh;
  int **firstneigh  = l_pm.list.firstneigh;
  int *firstnbr     = l_pm.H.firstnbr;
  int *numnbrs      = l_pm.H.numnbrs;
  int *H_jlist      = l_pm.H.jlist;
  double *val       = l_pm.H.val;
  double *Tap       = l_pm.Tap;

  m_fill = 0;
  r_sqr = 0;
  int ist, isz, ied, ioff;
  int mi[ISTEP];
  int nni[ISTEP], nbi[ISTEP], *fni[ISTEP];
  atom_in_t atom_in[ISTEP];
  rc_tagint l_tag[ISTEP];
  int jlist_buf[SNSTEP];
  
  int jst, jed, jsz;
  atom_in_t j_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_tag[JCACHE_LINECNT];
  for (i = 0; i < JCACHE_LINECNT; i ++)
    jcache_tag[i] = -1;

  volatile dma_desc get_desc = 0;
  volatile get_reply = 0;
  dma_set_mode(&get_desc, PE_MODE);
  dma_set_size(&get_desc, sizeof(atom_in_t)*JCACHE_LINESIZE);
  dma_set_op(&get_desc, DMA_GET);
  dma_set_reply(&get_desc, &get_reply);

  for(ist = _MYID * ISTEP; ist < nlocal; ist += ISTEP * 64)
  { 
    ied = ist + ISTEP;
    if(ied > inum)
      ied = inum;
    isz = ied - ist;

    pe_get(mask+ist,        mi,     sizeof(int)*isz);
    pe_get(numnbrs+ist,     nbi,    sizeof(int)*isz);
    pe_get(numneigh+ist,    nni,    sizeof(int)*isz);
    pe_get(firstneigh+ist,  fni,    sizeof(int*)*isz);
    pe_get(l_pm.atom_in+ist,atom_in,sizeof(atom_in_t)*isz);
    dma_syn();
    for(ii = ist; ii < ied; ii ++)
    {
      i = ii;
      ioff = i - ist;
      nbi[ioff] = 0;
      if (mi[ioff] & groupbit)
      {
        jlist = fni[ioff];
        jnum  = nni[ioff];
        for(jst = 0; jst < jnum; jst += SNSTEP)
        {
          jsz = SNSTEP;
          if(jst + SNSTEP > jnum)
            jsz = jnum - jst;
        
          pe_get(jlist+jst, jlist_buf, sizeof(int)*jsz);
          dma_syn();

          for(jj = 0; jj < jsz; jj++) 
          {
            j = jlist_buf[jj];
            int line = j >> JCACHE_SBIT & JCACHE_LMASK;
            int jtag = j >> JCACHE_HBIT;
            if(jcache_tag[line] != jtag)
            {
              int mem = j & ~JCACHE_MMASK;
              get_reply = 0;
              dma_rpl(get_desc, l_pm.atom_in+mem, j_cache[line], get_reply);
              while(get_reply != 1);
              jcache_tag[line] = jtag;
            }
            dx = j_cache[line][j & JCACHE_MMASK].x[0] - atom_in[ioff].x[0];
            dy = j_cache[line][j & JCACHE_MMASK].x[1] - atom_in[ioff].x[1];
            dz = j_cache[line][j & JCACHE_MMASK].x[2] - atom_in[ioff].x[2];
            r_sqr = SQR(dx) + SQR(dy) + SQR(dz);

            if(r_sqr <= swbsq) 
            {
              nbi[ioff] ++;
            }//if
          }//for-jj
        }//for-jst 
      }//if
    }//for-ii
    pe_put(numnbrs+ist, nbi, sizeof(int)*isz);
    dma_syn();
  }//for-ist
}

void compute_H_Full_C_para2(fix_qeq_pack_t *param)
{
  dma_init();
  fix_qeq_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(fix_qeq_pack_t));
  dma_syn();

  int i, j, ii, jj, flag,jnum;
  int *jlist;
  double dx, dy, dz, r_sqr;
  const double SMALL = 0.0001;
  double swbsq      = l_pm.swbsq;
  int groupbit      = l_pm.groupbit;
  int m_fill        = l_pm.m_fill;
  double **shld     = l_pm.shld;
  int ntypes        = l_pm.atom.ntypes;
  int *type         = l_pm.atom.type;
  rc_tagint *tag    = l_pm.atom.tag;
  double (*x)[3]    = l_pm.atom.x;
  int *mask         = l_pm.atom.mask;
  int nlocal        = l_pm.atom.nlocal;
  int inum          = l_pm.list.inum;  
  int *numneigh     = l_pm.list.numneigh;
  int **firstneigh  = l_pm.list.firstneigh;
  int *firstnbr     = l_pm.H.firstnbr;
  int *numnbrs      = l_pm.H.numnbrs;
  int *H_jlist      = l_pm.H.jlist;
  double *val       = l_pm.H.val;
  double *Tap       = l_pm.Tap;

  m_fill = 0; r_sqr = 0;
  int ist, isz, ied, ioff;
  int mi[ISTEP];
  int nni[ISTEP],*fni[ISTEP], fbi[ISTEP];
  int l_Hjlist[SNSTEP];
  double l_Hval[SNSTEP];
  atom_in_t atom_in[ISTEP];
  int jlist_buf[SNSTEP];

  double l_shld[ntypes+1][ntypes+1];
  double l_Tap[8];
  pe_get(l_pm.shld[0], l_shld, sizeof(double)*(ntypes+1)*(ntypes+1));
  pe_get(l_pm.Tap,     l_Tap,  sizeof(double)*8);
  dma_syn();

  int jst, jed, jsz;
  atom_in_t j_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_tag[JCACHE_LINECNT];
  for (i = 0; i < JCACHE_LINECNT; i ++)
    jcache_tag[i] = -1;

  volatile dma_desc get_desc = 0;
  volatile get_reply = 0;
  dma_set_mode(&get_desc, PE_MODE);
  dma_set_size(&get_desc, sizeof(atom_in_t)*JCACHE_LINESIZE);
  dma_set_op(&get_desc, DMA_GET);
  dma_set_reply(&get_desc, &get_reply);

  for(ist = _MYID * ISTEP; ist < nlocal; ist += ISTEP * 64)
  { 
    ied = ist + ISTEP;
    if(ied > inum)
      ied = inum;
    isz = ied - ist;
    pe_get(mask+ist,        mi,     sizeof(int)*isz);
    pe_get(numneigh+ist,    nni,    sizeof(int)*isz);
    pe_get(firstneigh+ist,  fni,    sizeof(int*)*isz);
    pe_get(firstnbr+ist,    fbi,    sizeof(int)*isz);
    pe_get(l_pm.atom_in+ist,atom_in,sizeof(atom_in_t)*isz);
    dma_syn();

    for(ii = ist; ii < ied; ii++) 
    {
      int foff = 0;//
      i = ii;
      ioff = i - ist;
      if(mi[ioff] & groupbit) 
      {
        jlist = fni[ioff];
        jnum  = nni[ioff];
        for(jst = 0; jst < jnum; jst += SNSTEP)
        {
          jsz = SNSTEP;
          if(jst + SNSTEP > jnum)
            jsz = jnum - jst;
          pe_get(jlist+jst, jlist_buf, sizeof(int)*jsz);
          //pe_get(H_jlist+fbi[ioff]+toff, l_Hjlist, sizeof(int)*)
          dma_syn();
          int tfoff = 0;//
          int ti = atom_in[ioff].type;
          for(jj = 0; jj < jsz; jj++) 
          {
            j = jlist_buf[jj];
            int line = j >> JCACHE_SBIT & JCACHE_LMASK;
            int jtag = j >> JCACHE_HBIT;
            if(jcache_tag[line] != jtag)
            {
              int mem = j & ~JCACHE_MMASK;
              get_reply = 0;
              dma_rpl(get_desc, l_pm.atom_in+mem, j_cache[line], get_reply);
              while(get_reply != 1);
              jcache_tag[line] = jtag;
            }
            dx = j_cache[line][j & JCACHE_MMASK].x[0] - atom_in[ioff].x[0];
            dy = j_cache[line][j & JCACHE_MMASK].x[1] - atom_in[ioff].x[1];
            dz = j_cache[line][j & JCACHE_MMASK].x[2] - atom_in[ioff].x[2];
            r_sqr = SQR(dx) + SQR(dy) + SQR(dz);
            int tj = j_cache[line][j & JCACHE_MMASK].type;
            if(r_sqr <= swbsq)
            {                             
              double Taper, denom;
              double r = sqrt(r_sqr);
              /* int ti, tj; */
              /* ti = atom_in[ioff].type; */
              /* tj = j_cache[line][j & JCACHE_MMASK].type; */
              //if ((ti < 0 || ti > 10 || tj < 0 || tj > 10)) call_printf("%d %d %d\n", j, ti, tj);
              double gamma = l_shld[ti][tj];
              Taper = l_Tap[7] * r + l_Tap[6];
              Taper = Taper * r + l_Tap[5];
              Taper = Taper * r + l_Tap[4];
              Taper = Taper * r + l_Tap[3];
              Taper = Taper * r + l_Tap[2];
              Taper = Taper * r + l_Tap[1];
              Taper = Taper * r + l_Tap[0];
              denom = r * r * r + gamma;
              denom = p_powd(denom,0.3333333333333);
              
              //H_jlist[fbi[ioff] + foff] = j;
              //H_jlist[fbi[ioff] + foff+tfoff] = j;
              l_Hjlist[tfoff] = j;

              //val[fbi[ioff] + foff] = Taper * EV_TO_KCAL_PER_MOL / denom;
              //val[fbi[ioff] + foff+tfoff] = Taper * EV_TO_KCAL_PER_MOL / denom;
              l_Hval[tfoff] = Taper * EV_TO_KCAL_PER_MOL / denom;
              tfoff++;
              //foff++;
            }//if
          }//for-jj
          if(tfoff>0)
          {
            pe_put(H_jlist+fbi[ioff]+foff, l_Hjlist, sizeof(int)*tfoff);
            pe_put(val+fbi[ioff]+foff,     l_Hval,   sizeof(double)*tfoff);
            dma_syn();
          }
          foff+=tfoff;

        }//for-jst
      }//if
      //if(_MYID == 0)
      //  printf("foff=%d\n",foff);

    }//for-ii
  }//for-ist
}

void compute_H_Full_C_para(fix_qeq_pack_t *param)
{
  dma_init();
  fix_qeq_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(fix_qeq_pack_t));
  dma_syn();

  int i, j, ii, jj, flag,jnum;
  int *jlist;
  double dx, dy, dz, r_sqr;
  double SMALL      = 0.0001;
  double swbsq      = l_pm.swbsq;
  int groupbit      = l_pm.groupbit;
  int m_fill        = l_pm.m_fill;
  double **shld     = l_pm.shld;
  int ntypes        = l_pm.atom.ntypes;
  int *type         = l_pm.atom.type;
  rc_tagint *tag    = l_pm.atom.tag;
  double (*x)[3]    = l_pm.atom.x;
  int *mask         = l_pm.atom.mask;
  int nlocal        = l_pm.atom.nlocal;
  int inum          = l_pm.list.inum;  
  int *numneigh     = l_pm.list.numneigh;
  int **firstneigh  = l_pm.list.firstneigh;
  int *firstnbr     = l_pm.H.firstnbr;
  int *numnbrs      = l_pm.H.numnbrs;
  int *H_jlist      = l_pm.H.jlist;
  double *val       = l_pm.H.val;
  double *Tap       = l_pm.Tap;
  int maxHlist      = l_pm.maxHlist;

  m_fill = 0; r_sqr = 0;
  int ist, isz, ied, ioff;
  int mi[ISTEP];
  int nni[ISTEP],*fni[ISTEP], fbi[ISTEP], nbi[ISTEP];
  int l_Hjlist[SNSTEP];
  double l_Hval[SNSTEP];
  atom_in_t atom_in[ISTEP];
  int jlist_buf[SNSTEP];

  double l_shld[ntypes+1][ntypes+1];
  double l_Tap[8];
  pe_get(l_pm.shld[0], l_shld, sizeof(double)*(ntypes+1)*(ntypes+1));
  pe_get(l_pm.Tap,     l_Tap,  sizeof(double)*8);
  dma_syn();

  int jst, jed, jsz;
  atom_in_t j_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_tag[JCACHE_LINECNT];
  for (i = 0; i < JCACHE_LINECNT; i ++)
    jcache_tag[i] = -1;

  volatile dma_desc get_desc = 0;
  volatile get_reply = 0;
  dma_set_mode(&get_desc, PE_MODE);
  dma_set_size(&get_desc, sizeof(atom_in_t)*JCACHE_LINESIZE);
  dma_set_op(&get_desc, DMA_GET);
  dma_set_reply(&get_desc, &get_reply);

  for(ist = _MYID * ISTEP; ist < nlocal; ist += ISTEP * 64)
  { 
    ied = ist + ISTEP;
    if(ied > inum)
      ied = inum;
    isz = ied - ist;
    pe_get(mask+ist,        mi,     sizeof(int)*isz);
    pe_get(numneigh+ist,    nni,    sizeof(int)*isz);
    pe_get(firstneigh+ist,  fni,    sizeof(int*)*isz);
    pe_get(l_pm.atom_in+ist,atom_in,sizeof(atom_in_t)*isz);
    dma_syn();

    for(ii = ist; ii < ied; ii++) 
    {
      int foff = 0;//
      i = ii;
      ioff = i - ist;
      nbi[ioff] = 0;
      fbi[ioff] = maxHlist * ii;
      if(mi[ioff] & groupbit) 
      {
        jlist = fni[ioff];
        jnum  = nni[ioff];
        for(jst = 0; jst < jnum; jst += SNSTEP)
        {
          jsz = SNSTEP;
          if(jst + SNSTEP > jnum)
            jsz = jnum - jst;
          pe_get(jlist+jst, jlist_buf, sizeof(int)*jsz);
          dma_syn();
          int tfoff = 0;//
          int ti = atom_in[ioff].type;
          for(jj = 0; jj < jsz; jj++) 
          {
            j = jlist_buf[jj];
            int line = j >> JCACHE_SBIT & JCACHE_LMASK;
            int jtag = j >> JCACHE_HBIT;
            if(jcache_tag[line] != jtag)
            {
              int mem = j & ~JCACHE_MMASK;
              get_reply = 0;
              dma_rpl(get_desc, l_pm.atom_in+mem, j_cache[line], get_reply);
              while(get_reply != 1);
              jcache_tag[line] = jtag;
            }
            dx = j_cache[line][j & JCACHE_MMASK].x[0] - atom_in[ioff].x[0];
            dy = j_cache[line][j & JCACHE_MMASK].x[1] - atom_in[ioff].x[1];
            dz = j_cache[line][j & JCACHE_MMASK].x[2] - atom_in[ioff].x[2];
            r_sqr = SQR(dx) + SQR(dy) + SQR(dz);
            int tj = j_cache[line][j & JCACHE_MMASK].type;
            if(r_sqr <= swbsq)
            {                             
              double Taper, denom;
              double r = sqrt(r_sqr);
              double gamma = l_shld[ti][tj];
              Taper = l_Tap[7] * r + l_Tap[6];
              Taper = Taper * r + l_Tap[5];
              Taper = Taper * r + l_Tap[4];
              Taper = Taper * r + l_Tap[3];
              Taper = Taper * r + l_Tap[2];
              Taper = Taper * r + l_Tap[1];
              Taper = Taper * r + l_Tap[0];
              denom = r * r * r + gamma;
              denom = p_powd(denom,0.3333333333333);
              l_Hjlist[tfoff] = j;
              l_Hval[tfoff] = Taper * EV_TO_KCAL_PER_MOL / denom;
              tfoff++;
              nbi[ioff] ++;
            }//if
          }//for-jj
          if(tfoff>0)
          {
            pe_put(H_jlist+fbi[ioff]+foff, l_Hjlist, sizeof(int)*tfoff);
            pe_put(val+fbi[ioff]+foff,     l_Hval,   sizeof(double)*tfoff);
            dma_syn();
          }
          foff+=tfoff;
        }//for-jst
      }//if
    }//for-ii
    pe_put(numnbrs+ist, nbi, sizeof(int)*isz);
    pe_put(firstnbr+ist,fbi, sizeof(int)*isz);
    dma_syn();
  }//for-ist
}
#endif
