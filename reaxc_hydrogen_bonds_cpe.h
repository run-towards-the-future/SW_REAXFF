#if (EVFLAG==0)
void Hydrogen_Bonds_C_ev0_cpe(hb_param_pack_t *param)
#else
void Hydrogen_Bonds_C_ev1_cpe(hb_param_pack_t *param)
#endif
{
  //lwpf_enter(HBOND);
  //lwpf_start(ALL);
  //lwpf_start(BEFORE);

  dma_init();
  hb_param_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(hb_param_pack_t));
  dma_syn();
  
  reax_system_c *system         = l_pm.system;
  control_params *control       = l_pm.control;
  simulation_data *data         = l_pm.data;
  storage *workspace            = l_pm.workspace;
  reax_list **lists             = l_pm.lists;
  output_controls *out_control  = l_pm.out_control;

  int  i, j, k, pi, pk, pi_off, pk_off;
  int  type_i, type_j, type_k;
  int  start_j, end_j, len_j, hb_start_j, hb_end_j, hb_len_j;
  int  hblist[MAX_BONDS];
  int  itr, top;
  ivec rel_jk;
  double r_jk, theta, cos_theta, sin_xhz4, cos_xhz1, sin_theta2;
  double e_hb, exp_hb2, exp_hb3, CEhb1, CEhb2, CEhb3;
  rvec dcos_theta_di, dcos_theta_dj, dcos_theta_dk;
  rvec dvec_jk, force, ext_press;
  hbond_parameters *hbp;
  bond_order_data *bo_ij;
  bond_data *pbond_ij;
  far_neighbor_data_full *nbr_jk;

  // tally variables
  double fi_tmp[3], fk_tmp[3], delij[3], delkj[3];
  
  bond_data *bond_list          = l_pm.bond_list;
  hbond_data *hbond_list        = l_pm.hbond_list;
  double *Cdbo_list             = l_pm.Cdbo_list;
  double *BO_list               = l_pm.BO_list;
  rvec4  *fCdDelta              = l_pm.fCdDelta;
  atom_pack_t *packed_atoms     = l_pm.packed_atoms;
  int n                         = l_pm.n;
  int ntypes                    = l_pm.ntypes;
  int *Hindex                   = l_pm.Hindex;
  int *bindex                   = l_pm.bindex;
  int *end_bindex               = l_pm.end_bindex;
  int *hbindex                  = l_pm.hbindex;
  int *end_hbindex              = l_pm.end_hbindex;


  int jst, jed, jsz, joff;
  rvec4 fcj_tmp, fck_tmp, fci_tmp;
  swcache_lock_t *locks_frc = l_pm.locks_frc;
  SWCACHE_INIT_U(DT_CACHE_F, locks_frc);

  doublev4 eng_virial[2];
  eng_virial[0] = eng_virial[1] = 0;
  double *eng_vdwl = (double*)(void*)eng_virial;
  double *eng_ehb  = eng_vdwl + 1;

  atom_pack_t atom_tmp[BLEN];
  bond_data blist_j[BLEN];
  hbond_data hblist_j[HLEN];
  double Cdbo_list_j[BLEN] ,BO_list_j[BLEN];
  int itype[BLEN];

  atom_pack_t atom_j[JSTEP];
  int bidx[JSTEP], end_bidx[JSTEP], Hidx[JSTEP];

  //read_cache;
  int read_ctag[READ_C_LCNT];
  atom_pack_t atoms_cache[READ_C_LCNT][READ_C_LSZ];
  atom_pack_t atom_k, atom_i;
  for(i = 0; i < READ_C_LCNT; i++)
  {
    read_ctag[i] = -1;
  }
  //read index for posion k;
  int hb_tag[K_C_LCNT];
  int hbidx_cache[K_C_LCNT][K_C_LSZ];
  int end_hbidx_cache[K_C_LCNT][K_C_LSZ];
  for(i = 0; i < K_C_LCNT; i++)
  {
    hb_tag[i] = -1;
  }

  //lwpf_stop(BEFORE);

  //lwpf_start(CMP);
  //if(_MYID == 0)
  for( jst = _MYID * JSTEP; jst < n; jst+=64*JSTEP)
  {
    jed = jst + JSTEP;
    if(jed > n)
      jed = n;
    jsz = jed - jst;
    
    pe_get(l_pm.bindex+jst,       bidx,       sizeof(int) * jsz);
    pe_get(l_pm.end_bindex+jst,   end_bidx,   sizeof(int) * jsz);
    pe_get(l_pm.Hindex+jst,       Hidx,       sizeof(int) * jsz);
    pe_get(packed_atoms+jst,      atom_j,     sizeof(atom_pack_t) * jsz);
    dma_syn();

    for(j = jst; j < jed; j++)
    {
      joff = j - jst;
      type_j     = atom_j[joff].type;
      if(l_pm.sbp[type_j] == 1 ) 
      {
        //lwpf_start(IF);
        if (type_j < 0) continue;
        start_j     = bidx[joff];
        end_j       = end_bidx[joff];
        len_j       = end_j - start_j;
        
        read_cache_index(Hidx[joff], hbidx_cache, end_hbidx_cache, hb_tag, 
                        l_pm.hbindex, l_pm.end_hbindex, &hb_start_j, &hb_end_j);

        hb_len_j    = hb_end_j - hb_start_j;
       
        pe_get(bond_list+start_j,     blist_j,      sizeof(bond_data) * len_j);
        pe_get(BO_list+start_j,       BO_list_j,    sizeof(double) * len_j);
        pe_get(Cdbo_list+start_j,     Cdbo_list_j,  sizeof(double) * len_j);
        pe_get(hbond_list+hb_start_j, hblist_j,     sizeof(hbond_data) * hb_len_j);
        dma_syn();

        fcj_tmp[0] = fcj_tmp[1] = fcj_tmp[2] = fcj_tmp[3] = 0;

        //lwpf_start(BOND);
        top = 0;
        for( pi = start_j; pi < end_j; ++pi )  
        {
          pi_off    = pi - start_j;
          pbond_ij  = &(blist_j[pi_off] );
          i = pbond_ij->nbr;
          read_cache(i, atoms_cache, read_ctag, packed_atoms, &atom_i);
          //type_i = atom_i.type;
          atom_tmp[pi_off] = atom_i;///
	        //if (type_i < 0) continue;
	        if (atom_i.type < 0) continue;
          if(l_pm.sbp[atom_i.type] == 2 && BO_list_j[pi_off] >= HB_THRESHOLD)
            hblist[top++] = pi;
        }
        //lwpf_stop(BOND);

        //lwpf_start(PK);
        for( pk = hb_start_j; pk < hb_end_j; ++pk ) 
        {
          pk_off = pk - hb_start_j;
          /* set k's varibles */
          k = hblist_j[pk_off].nbr;
          read_cache(k, atoms_cache, read_ctag, packed_atoms, &atom_k);
          type_k = atom_k.type;
	        if (type_k < 0) continue;
          nbr_jk  = hblist_j[pk_off].ptr;
          r_jk    = nbr_jk->d;
          rvec_Scale(dvec_jk, hblist_j[pk_off].scl, nbr_jk->dvec );
          
          fck_tmp[0] = fck_tmp[1] = fck_tmp[2] = fck_tmp[3] = 0;

          for( itr = 0; itr < top; ++itr ) 
          {
            pi        = hblist[itr];
            pi_off    = pi - start_j;
            pbond_ij  = &(blist_j[pi_off]);
            i         = pbond_ij->nbr;
            //read_cache(i, atoms_cache, read_ctag, packed_atoms, &atom_i);
            
            //if(atom_i.orig_id != atom_k.orig_id) 
            if(atom_tmp[pi_off].orig_id != atom_k.orig_id) 
            {
              //type_i = atom_i.type;
              type_i = atom_tmp[pi_off].type;
	            if (type_i < 0) continue;

              int pos = (type_i * ntypes + type_j)*ntypes+type_k;
              hbp = &(l_pm.hbp[pos]);
	            if (hbp->r0_hb <= 0.0) continue;

              Calculate_Theta( pbond_ij->dvec, pbond_ij->d, dvec_jk, r_jk,
                               &theta, &cos_theta );
              /* the derivative of cos(theta) */
              Calculate_dCos_Theta(pbond_ij->dvec, pbond_ij->d, dvec_jk, r_jk,
                                   &dcos_theta_di, &dcos_theta_dj, &dcos_theta_dk );

              /* hyrogen bond energy*/
              sin_theta2 = p_sind( theta/2.0 );
              sin_xhz4 = SQR(sin_theta2);
              sin_xhz4 *= sin_xhz4;
              cos_xhz1 = ( 1.0 - cos_theta );
              exp_hb2 = p_expd(-hbp->p_hb2 * BO_list_j[pi_off]);
              exp_hb3 = p_expd(-hbp->p_hb3 * (hbp->r0_hb / r_jk + r_jk / hbp->r0_hb - 2.0));

              *eng_ehb += e_hb = hbp->p_hb1 * (1.0 - exp_hb2) * exp_hb3 * sin_xhz4;

              CEhb1 = hbp->p_hb1 * hbp->p_hb2 * exp_hb2 * exp_hb3 * sin_xhz4;
              CEhb2 = -hbp->p_hb1/2.0 * (1.0 - exp_hb2) * exp_hb3 * cos_xhz1;
              CEhb3 = -hbp->p_hb3 * (-hbp->r0_hb / SQR(r_jk) + 1.0 / hbp->r0_hb) * e_hb;

              /* hydrogen bond forces */
              Cdbo_list_j[pi_off] += CEhb1; // dbo term

              //// dcos terms
              fci_tmp[0] = fci_tmp[1] = fci_tmp[2] = fci_tmp[3] = 0;
              
              rvec4_ScaledAdd(fci_tmp, +CEhb2, dcos_theta_di );
              rvec4_ScaledAdd(fcj_tmp, +CEhb2, dcos_theta_dj );
              rvec4_ScaledAdd(fck_tmp, +CEhb2, dcos_theta_dk );
              rvec4_ScaledAdd(fcj_tmp, -CEhb3/r_jk, dvec_jk );
              rvec4_ScaledAdd(fck_tmp, +CEhb3/r_jk, dvec_jk );

              SWCACHE_UPDATE(DT_CACHE_F, i, (&fci_tmp));
              
              #if(EVFLAG)
              *eng_vdwl += e_hb;
              #endif
            }//if
          }//for-itr

          SWCACHE_UPDATE(DT_CACHE_F, k, (&fck_tmp));
        }//for-pk
        //lwpf_stop(PK);
        
        SWCACHE_UPDATE(DT_CACHE_F, j, (&fcj_tmp));
        pe_put(Cdbo_list+start_j, Cdbo_list_j, sizeof(double) * len_j);
        dma_syn();

        //lwpf_stop(IF);
      }//if

    }//for-j
  }//for-jst
  SWCACHE_FLUSH(DT_CACHE_F);
  
  //lwpf_stop(CMP);
  /******reduce energy*****/
  reg_reduce_inplace_doublev4(eng_virial, 1);
  if(_MYID == 0)
  {
    pe_put(l_pm.packed_eng, eng_virial, sizeof(double) * 4);
    dma_syn();
  }

  //lwpf_stop(ALL);
  //lwpf_exit(HBOND);
}

