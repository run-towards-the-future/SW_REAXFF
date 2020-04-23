#if(EVFLAG)
void vdW_Coulomb_Energy_cpe_evflag1(vdw_coulomb_pack_t * param)
#else
void vdW_Coulomb_Energy_cpe_evflag0(vdw_coulomb_pack_t * param)
#endif
{
  //lwpf_enter(NONBOND);
  //lwpf_start(ALL);
  //lwpf_start(BEFORE);
  dma_init();
  vdw_coulomb_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(vdw_coulomb_pack_t));
  dma_syn();
  reax_system_c *system       = l_pm.system;
  control_params *control     = l_pm.control;
  simulation_data *data       = l_pm.data;
  storage *workspace          = l_pm.workspace;
  reax_list **lists           = l_pm.lists;
  
  int i;
  //int j;
  int pj, natoms, nlocal;
  int start_i, end_i;//, flag;
  rc_tagint orig_i, orig_j;
  double p_vdW1;//, p_vdW1i;
  double powr_vdW1, powgi_vdW1;
  double tmp, r_ij, fn13, exp1, exp2;
  double Tap, dTap, dfn13, CEvd, CEclmb, de_core;
  double dr3gamij_1, dr3gamij_3;
  double e_ele, e_vdW, e_core, SMALL = 0.0001;
  double e_lg, de_lg, r_ij5, r_ij6, re6;
  rvec temp, ext_press;
  far_neighbor_data_full *nbr_pj;
  reax_list *far_nbrs;

  // Tallying variables:
  double pe_vdw, f_tmp, delij[3];
  natoms = l_pm.N;
  nlocal = l_pm.n;
  far_nbrs = (*lists) + FAR_NBRS_FULL;
  p_vdW1 = l_pm.l28;
  //p_vdW1i = 1.0 / p_vdW1;
  e_core = 0;
  e_vdW = 0;
  e_lg = de_lg = 0.0;

  int ntypes = l_pm.num_atom_types;
  tbp_boost_t tbpb[ntypes][ntypes];
  tbp_boost_t *tbpt;
  vdw_coulomb_tbp l_tbpr[ntypes][ntypes];
  vdw_coulomb_tbp *twbp;

  int ist, isz, ied;
  int ii, ioff;
  int idx_st[ISTEP], idx_ed[ISTEP];
  double ei[ISTEP], vi[ISTEP][6];
  double fi[ISTEP][4];
  atom_pack_t packed_atoms[ISTEP];
  //far_neighbor_data_full pj_nbr[SNSTEP];
  doublev4  nbr_v4[SNSTEP*2], *p_nbr_v4;
  far_neighbor_data_full *pj_nbr = nbr_v4;
  
  doublev4 t_powgi_vdW1_v4[NTYP];
  doublev4 t_r_vdWinv_v4[NTYP];
  doublev4 t_alpha_div_rvdW_v4[NTYP];
  doublev4 t_alpha_v4[NTYP];
  doublev4 t_gamma_v4[NTYP];
  doublev4 t_D_v4[NTYP];
  
  double *p_powgi_vdW1      = t_powgi_vdW1_v4;
  double *p_r_vdWinv        = t_r_vdWinv_v4;
  double *p_alpha_div_rvdW  = t_alpha_div_rvdW_v4;
  double *p_alpha           = t_alpha_v4;
  double *p_gamma           = t_gamma_v4;
  double *p_D               = t_D_v4;

  double *params = l_pm.params;

  //pe_get(l_pm.tbpb, tbpb, sizeof(tbp_boost_t) * ntypes * ntypes);
  //pe_get(l_pm.pack_tbpr, l_tbpr, sizeof(vdw_coulomb_tbp) * ntypes * ntypes);
  pe_get(l_pm.powgi_vdW1, p_powgi_vdW1, sizeof(double) * ntypes * ntypes);
  pe_get(l_pm.r_vdWinv      , p_r_vdWinv      , sizeof(double) * ntypes * ntypes);
  pe_get(l_pm.alpha_div_rvdW, p_alpha_div_rvdW, sizeof(double) * ntypes * ntypes);
  pe_get(l_pm.alpha         , p_alpha         , sizeof(double) * ntypes * ntypes);
  pe_get(l_pm.gamma         , p_gamma         , sizeof(double) * ntypes * ntypes);
  pe_get(l_pm.D             , p_D             , sizeof(double) * ntypes * ntypes);
  
  dma_syn();

  doublev4 eng_virial[4];
  for(i = 0; i < 4; i++)
    eng_virial[i] = 0;
  double *eng_vdwl = (double*)(void*)eng_virial;
  double *eng_coul = eng_vdwl + 1; 
  double *virial   = eng_coul + 1;
  double *extp     = virial + 6;
  double *data_vdW = extp + 3;
  double *data_ele = data_vdW + 1;
  
  volatile dma_desc get_desc = 0;
  volatile get_reply = 0;
  dma_set_mode(&get_desc, PE_MODE);
  dma_set_op(&get_desc, DMA_GET);
  dma_set_reply(&get_desc, &get_reply);

  //vectorize
  doublev4 r_ij_v4, rijinv_v4;
  doublev4 fn13_v4, dfn13_v4;
  doublev4 exp1_v4, exp2_v4;
  doublev4 evdW_v4;

  doublev4 dr3gamij_1_v4, dr3gamij_3_v4;
  doublev4 dr3gamij_1inv_v4, dr3gamij_3inv_v4;
  doublev4 CEvd_v4, CEclmb_v4;
  doublev4 tmp_v4, tmpc_v4, fpair_v4;
  doublev4 evdwl_v4, ecoul_v4;


  doublev4 alpha_v4;
  doublev4 D_v4, gamma_v4;
  doublev4 powgi_vdW1_v4;
  doublev4 r_vdWinv_v4;
  doublev4 alpha_div_rvdW_v4;
  doublev4 del0_v4, del1_v4, del2_v4;
  doublev4 qi_v4, qj_v4;
  doublev4 half_v4, one_v4, two_v4, thr_v4, four_v4;
  doublev4 five_v4, six_v4, seven_v4, thrinv_v4;

  //half_v4 = 0.5; one_v4 = 1.0; two_v4 = 2.0; thrinv_v4 = 0.3333333333;
  //thr_v4 = 3; four_v4 = 4; five_v4 = 5; six_v4 = 6; seven_v4 = 7;
  half_v4   =simd_set_doublev4(0.5, 0.5, 0.5, 0.5);
  one_v4    =simd_set_doublev4(1.0, 1.0, 1.0, 1.0);
  two_v4    =simd_set_doublev4(2.0, 2.0, 2.0, 2.0);
  thrinv_v4 =simd_set_doublev4(0.333333, 0.333333, 0.333333, 0.333333);
  thr_v4    =simd_set_doublev4(3.0, 3.0, 3.0, 3.0);
  four_v4   =simd_set_doublev4(4.0, 4.0, 4.0, 4.0);
  five_v4   =simd_set_doublev4(5.0, 5.0, 5.0, 5.0);
  six_v4    =simd_set_doublev4(6.0, 6.0, 6.0, 6.0);
  seven_v4  =simd_set_doublev4(7.0, 7.0, 7.0, 7.0);

  doublev4 remd_v4[4];//remainder;
  remd_v4[0] = simd_set_doublev4(1.0, 1.0, 1.0, 1.0);
  remd_v4[1] = simd_set_doublev4(1.0, 0.0, 0.0, 0.0);
  remd_v4[2] = simd_set_doublev4(1.0, 1.0, 0.0, 0.0);
  remd_v4[3] = simd_set_doublev4(1.0, 1.0, 1.0, 0.0);
  
  //doublev4 test_remd_v4[4];//remainder;
  //test_remd_v4[0] = simd_set_doublev4(0.0, 0.0, 0.0, 0.0);//(1.0, 1.0, 1.0, 1.0);
  //test_remd_v4[1] = simd_set_doublev4(0.0, 1.0, 1.0, 1.0);//(1.0, 0.0, 0.0, 0.0);
  //test_remd_v4[2] = simd_set_doublev4(0.0, 0.0, 1.0, 1.0);//(1.0, 1.0, 0.0, 0.0);
  //test_remd_v4[3] = simd_set_doublev4(0.0, 0.0, 0.0, 1.0);//(1.0, 1.0, 1.0, 0.0);

  doublev4 Cele_v4, p_vdW1_v4, p_vdW1i_v4, powr_vdW1_v4;
  Cele_v4 = simd_set_doublev4(332.06371, 332.06371, 332.06371, 332.06371);//C_ele;  
  p_vdW1_v4 = l_pm.l28;  
  //p_vdW1i_v4 = one_v4 / p_vdW1_v4; //1.0 / l_pm.l28;
  p_vdW1i_v4 = simd_vdivd(one_v4, p_vdW1_v4); //1.0 / l_pm.l28;
  doublev4 power_onev4 = p_vdW1i_v4-one_v4;
  doublev4 power_twov4 = p_vdW1_v4 -two_v4;

  doublev4 Tap_v4, dTap_v4;

  doublev4 Tap0_v4, Tap1_v4, Tap2_v4, Tap3_v4, Tap4_v4, Tap5_v4, Tap6_v4, Tap7_v4; 
  //Tap0_v4 = l_pm.wksp_tap0;
  //Tap1_v4 = l_pm.wksp_tap1;
  //Tap2_v4 = l_pm.wksp_tap2;
  //Tap3_v4 = l_pm.wksp_tap3;
  //Tap4_v4 = l_pm.wksp_tap4;
  //Tap5_v4 = l_pm.wksp_tap5;
  //Tap6_v4 = l_pm.wksp_tap6;
  //Tap7_v4 = l_pm.wksp_tap7;
  
  Tap0_v4 = simd_set_doublev4(l_pm.wksp_tap0,l_pm.wksp_tap0,l_pm.wksp_tap0,l_pm.wksp_tap0);
  Tap1_v4 = simd_set_doublev4(l_pm.wksp_tap1,l_pm.wksp_tap1,l_pm.wksp_tap1,l_pm.wksp_tap1);
  Tap2_v4 = simd_set_doublev4(l_pm.wksp_tap2,l_pm.wksp_tap2,l_pm.wksp_tap2,l_pm.wksp_tap2);
  Tap3_v4 = simd_set_doublev4(l_pm.wksp_tap3,l_pm.wksp_tap3,l_pm.wksp_tap3,l_pm.wksp_tap3);
  Tap4_v4 = simd_set_doublev4(l_pm.wksp_tap4,l_pm.wksp_tap4,l_pm.wksp_tap4,l_pm.wksp_tap4);
  Tap5_v4 = simd_set_doublev4(l_pm.wksp_tap5,l_pm.wksp_tap5,l_pm.wksp_tap5,l_pm.wksp_tap5);
  Tap6_v4 = simd_set_doublev4(l_pm.wksp_tap6,l_pm.wksp_tap6,l_pm.wksp_tap6,l_pm.wksp_tap6);
  Tap7_v4 = simd_set_doublev4(l_pm.wksp_tap7,l_pm.wksp_tap7,l_pm.wksp_tap7,l_pm.wksp_tap7);

  doublev4 eng_vdwl_v4, eng_coul_v4, data_vdW_v4, data_ele_v4;
  eng_vdwl_v4 = eng_coul_v4 = 0;
  data_vdW_v4 = data_ele_v4 = 0;
  doublev4 virial0_v4, virial1_v4, virial2_v4;
  doublev4 virial3_v4, virial4_v4, virial5_v4;
  virial0_v4 = virial1_v4 = virial2_v4 = 0;
  virial3_v4 = virial4_v4 = virial5_v4 = 0;
  //lwpf_stop(BEFORE);

  //lwpf_start(CMP);
  for(ist = _MYID * ISTEP; ist < nlocal; ist += ISTEP * 64)
  {
    ied = ist + ISTEP;
    if(ied > nlocal)
      ied = nlocal;
    isz = ied - ist;

    pe_get(l_pm.index+ist, idx_st, sizeof(int)*isz);
    pe_get(l_pm.end_index+ist, idx_ed, sizeof(int)*isz);
    pe_get(l_pm.packed_atoms+ist, packed_atoms, sizeof(atom_pack_t)*isz);
    pe_get(l_pm.f[ist], fi[0], sizeof(double)*isz*4);
    dma_syn();

    for( ii = ist; ii < ied; ++ii) 
    {
      i = ii;
      int ioff = i - ist;
      if (packed_atoms[ioff].type < 0) continue;
      start_i = idx_st[ioff];
      end_i   = idx_ed[ioff];
      orig_i  = packed_atoms[ioff].orig_id;
      int type_i = packed_atoms[ioff].type;
      int jst, jed, jsz, jsz_to;
      
      doublev4 fi0_v4, fi1_v4, fi2_v4;
      fi0_v4 = fi1_v4 = fi2_v4 = 0;
      
      //qi_v4 = packed_atoms[ioff].q; 
      qi_v4 = simd_set_doublev4(packed_atoms[ioff].q, packed_atoms[ioff].q, 
                                packed_atoms[ioff].q, packed_atoms[ioff].q);

      for(jst = start_i; jst < end_i; jst+=SNSTEP)
      {
        jsz = SNSTEP;
        if(jst + SNSTEP > end_i)
          jsz = end_i - jst;
        
        get_reply = 0;
        dma_set_size(&get_desc, sizeof(far_neighbor_data_full)*jsz);
        dma(get_desc,l_pm.far_nbr_list_full+jst, pj_nbr);
        while(get_reply != 1);
        jsz_to = jsz & ~3;
        int pjst, tt;
            
        //lwpf_start(SEC1);
        for(pjst = 0; pjst < jsz_to; pjst+=4)
        {
          //transpose
          //lwpf_start(TRANS);
          doublev4 trans_nbr[3];//trans_nbr[8];
          p_nbr_v4 = &nbr_v4[pjst * 2];

          transpose8x4( p_nbr_v4[0],p_nbr_v4[2],p_nbr_v4[4],p_nbr_v4[6],
                        p_nbr_v4[1],p_nbr_v4[3],p_nbr_v4[5],p_nbr_v4[7], 
                        trans_nbr[0],trans_nbr[1],trans_nbr[2],r_ij_v4,
                        qj_v4, del0_v4, del1_v4, del2_v4);
          int *pt = &trans_nbr[0];
          int trc = vshuffd_rc(pt[7], pt[5], pt[3], pt[1]);//(typej3, typej2, typej1, typej0);

          //vshuffle
          alpha_v4 = simd_vshff(t_alpha_v4[type_i],t_alpha_v4[type_i],trc); 
          D_v4     = simd_vshff(t_D_v4[type_i],t_D_v4[type_i],trc); 
          gamma_v4 = simd_vshff(t_gamma_v4[type_i],t_gamma_v4[type_i],trc); 
          powgi_vdW1_v4     = simd_vshff(t_powgi_vdW1_v4[type_i],t_powgi_vdW1_v4[type_i],trc); 
          r_vdWinv_v4       = simd_vshff(t_r_vdWinv_v4[type_i],t_r_vdWinv_v4[type_i],trc); 
          alpha_div_rvdW_v4 = simd_vshff(t_alpha_div_rvdW_v4[type_i],t_alpha_div_rvdW_v4[type_i],trc); 
      
          //lwpf_stop(TRANS);
          //lwpf_start(VMAD1);
          //vectorize
          //rijinv_v4 = one_v4 / r_ij_v4;
          rijinv_v4 = simd_vdivd(one_v4, r_ij_v4);
          
          //Tap_v4 = Tap7_v4 * r_ij_v4 + Tap6_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap5_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap4_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap3_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap2_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap1_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap0_v4;//l_pm.wksp_tap0;//Tap0_v4;
          
          Tap_v4 = simd_vmad(Tap7_v4, r_ij_v4, Tap6_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap5_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap4_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap3_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap2_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap1_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap0_v4);//l_pm.wksp_tap0;//Tap0_v4;

          //dTap_v4 = seven_v4 * Tap7_v4 * r_ij_v4 + six_v4  * Tap6_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + five_v4 * Tap5_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + four_v4 * Tap4_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + thr_v4  * Tap3_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + two_v4  * Tap2_v4;
          //dTap_v4 += Tap1_v4 * rijinv_v4;
          dTap_v4 = simd_vmad(seven_v4, Tap7_v4 * r_ij_v4, six_v4  * Tap6_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, five_v4 * Tap5_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, four_v4 * Tap4_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, thr_v4  * Tap3_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, two_v4  * Tap2_v4);
          dTap_v4 = simd_vmad(Tap1_v4, rijinv_v4, dTap_v4);

          //lwpf_stop(VMAD1);

          //lwpf_start(MATH);
         
          //powr_vdW1_v4 = simd_vpowd(r_ij_v4, p_vdW1_v4);
          doublev4 ln_rij_v4 = simd_vlnd(r_ij_v4);
          powr_vdW1_v4 = simd_vexpd(p_vdW1_v4 * ln_rij_v4);

          doublev4 base_vdW = powr_vdW1_v4 + powgi_vdW1_v4;
          
          //fn13_v4 = simd_vpowd(powr_vdW1_v4 + powgi_vdW1_v4, p_vdW1i_v4);
          //fn13_v4 = simd_vpowd(base_vdW, p_vdW1i_v4);
          doublev4 ln_vdW_v4 = simd_vlnd(base_vdW);
          fn13_v4 =  simd_vexpd(p_vdW1i_v4 * ln_vdW_v4);

          exp2_v4 = simd_vexpd(half_v4 * alpha_v4 * (one_v4 - fn13_v4 * r_vdWinv_v4));
          exp1_v4 = exp2_v4 * exp2_v4;
          evdW_v4 = D_v4 * (exp1_v4 - two_v4 * exp2_v4);
          
          //data_vdW_v4 += Tap_v4 * evdW_v4;
          data_vdW_v4 = simd_vmad(Tap_v4, evdW_v4, data_vdW_v4);
          
          //lwpf_start(POW);
          //dfn13_v4 = simd_vpowd(powr_vdW1_v4 + powgi_vdW1_v4, p_vdW1i_v4-one_v4) 
          //            * simd_vpowd(r_ij_v4, p_vdW1_v4 - two_v4);
          
          //dfn13_v4 = simd_vpowd(base_vdW, power_onev4)  * ans_twov4;
          //dfn13_v4 = simd_vpowd(base_vdW, power_onev4)* simd_vpowd(r_ij_v4, power_twov4);
          dfn13_v4 = simd_vexpd(power_onev4 *ln_vdW_v4)* simd_vexpd(power_twov4 * ln_rij_v4);

          //lwpf_stop(POW);
          //lwpf_start(POW2) ;
          CEvd_v4 = dTap_v4 * evdW_v4 - Tap_v4 * D_v4 * alpha_div_rvdW_v4 
                    * (exp1_v4-exp2_v4) *dfn13_v4;
          
          //dr3gamij_1_v4 = (r_ij_v4 * r_ij_v4 * r_ij_v4 + gamma_v4);
          dr3gamij_1_v4 = simd_vmad(r_ij_v4, r_ij_v4 * r_ij_v4, gamma_v4);
          dr3gamij_3_v4 = simd_vpowd(dr3gamij_1_v4, thrinv_v4);
          
          //lwpf_stop(POW2) ;
          //dr3gamij_1inv_v4 = one_v4 / dr3gamij_1_v4;
          //dr3gamij_3inv_v4 = one_v4 / dr3gamij_3_v4;
          //lwpf_start(DIV);
          dr3gamij_1inv_v4 = simd_vdivd(one_v4, dr3gamij_1_v4);
          dr3gamij_3inv_v4 = simd_vdivd(one_v4, dr3gamij_3_v4);
          //lwpf_stop(DIV);

          //lwpf_stop(MATH);
          //
          //lwpf_start(VMAD);
          tmp_v4 = Tap_v4 * dr3gamij_3inv_v4;
          doublev4 Cqij_v4 = Cele_v4 * qi_v4 * qj_v4;
          data_ele_v4 += ecoul_v4 = Cqij_v4 * tmp_v4;
          
          CEclmb_v4 = Cqij_v4 * (dTap_v4 -  Tap_v4 * r_ij_v4 * dr3gamij_1inv_v4) 
                     * dr3gamij_3inv_v4;
            
          evdwl_v4 = Tap_v4 * evdW_v4;
          fpair_v4 = -(CEvd_v4 + CEclmb_v4);
          //tmpc_v4 = CEvd_v4 + CEclmb_v4;
          
          //eng_vdwl_v4 += (half_v4 *evdwl_v4);
          //eng_coul_v4 += (half_v4 *ecoul_v4);
          #if(EVFLAG)
          eng_vdwl_v4 = simd_vmad(half_v4, evdwl_v4, eng_vdwl_v4);
          eng_coul_v4 = simd_vmad(half_v4, ecoul_v4, eng_coul_v4);
          #endif

          //virial0_v4 += (half_v4 * del0_v4 * del0_v4 * fpair_v4);
          //virial1_v4 += (half_v4 * del1_v4 * del1_v4 * fpair_v4);
          //virial2_v4 += (half_v4 * del2_v4 * del2_v4 * fpair_v4);
          //virial3_v4 += (half_v4 * del0_v4 * del1_v4 * fpair_v4);
          //virial4_v4 += (half_v4 * del0_v4 * del2_v4 * fpair_v4);
          //virial5_v4 += (half_v4 * del1_v4 * del2_v4 * fpair_v4);
          doublev4 half_fpair = half_v4 * fpair_v4;
          virial0_v4 = simd_vmad(half_fpair, del0_v4 * del0_v4, virial0_v4);
          virial1_v4 = simd_vmad(half_fpair, del1_v4 * del1_v4, virial1_v4);
          virial2_v4 = simd_vmad(half_fpair, del2_v4 * del2_v4, virial2_v4);
          virial3_v4 = simd_vmad(half_fpair, del0_v4 * del1_v4, virial3_v4);
          virial4_v4 = simd_vmad(half_fpair, del0_v4 * del2_v4, virial4_v4);
          virial5_v4 = simd_vmad(half_fpair, del1_v4 * del2_v4, virial5_v4);

          //fi0_v4 += (-tmpc_v4 * del0_v4);
          //fi1_v4 += (-tmpc_v4 * del1_v4);
          //fi2_v4 += (-tmpc_v4 * del2_v4);
          fi0_v4 = simd_vmad(fpair_v4, del0_v4, fi0_v4);
          fi1_v4 = simd_vmad(fpair_v4, del1_v4, fi1_v4);
          fi2_v4 = simd_vmad(fpair_v4, del2_v4, fi2_v4);
          //lwpf_stop(VMAD);
        }//for-pjst
        //lwpf_stop(SEC1);

        //lwpf_start(REMAINDER);
        if(jsz_to < jsz)
        {
          /**********remainder*********/
          //transpose
          doublev4 trans_nbr[3];//trans_nbr[8];
          p_nbr_v4 = &nbr_v4[jsz_to * 2];

          transpose8x4( p_nbr_v4[0],p_nbr_v4[2],p_nbr_v4[4],p_nbr_v4[6],
                        p_nbr_v4[1],p_nbr_v4[3],p_nbr_v4[5],p_nbr_v4[7], 
                        trans_nbr[0],trans_nbr[1],trans_nbr[2],r_ij_v4,
                        qj_v4, del0_v4, del1_v4, del2_v4);
          int *pt = &trans_nbr[0];
          int trc = vshuffd_rc(pt[7], pt[5], pt[3], pt[1]);//(typej3, typej2, typej1, typej0);

          //vshuffle
          alpha_v4 = simd_vshff(t_alpha_v4[type_i],t_alpha_v4[type_i],trc); 
          D_v4     = simd_vshff(t_D_v4[type_i],t_D_v4[type_i],trc); 
          gamma_v4 = simd_vshff(t_gamma_v4[type_i],t_gamma_v4[type_i],trc); 
          powgi_vdW1_v4     = simd_vshff(t_powgi_vdW1_v4[type_i],t_powgi_vdW1_v4[type_i],trc); 
          r_vdWinv_v4       = simd_vshff(t_r_vdWinv_v4[type_i],t_r_vdWinv_v4[type_i],trc); 
          alpha_div_rvdW_v4 = simd_vshff(t_alpha_div_rvdW_v4[type_i],t_alpha_div_rvdW_v4[type_i],trc); 
                
          //vectorize
          //rijinv_v4 = one_v4 / r_ij_v4;
          rijinv_v4 = simd_vdivd(one_v4, r_ij_v4);
          //rijinv_v4 = simd_vsellt(test_remd_v4[jsz &3], one_v4 / r_ij_v4, test_remd_v4[0]);//* remd_v4[jsz & 3];

          //Tap_v4 = Tap7_v4 * r_ij_v4 + Tap6_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap5_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap4_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap3_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap2_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap1_v4;
          //Tap_v4 = Tap_v4 * r_ij_v4 + Tap0_v4;//l_pm.wksp_tap0;//Tap0_v4;
          //dTap_v4 = seven_v4 * Tap7_v4 * r_ij_v4 + six_v4  * Tap6_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + five_v4 * Tap5_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + four_v4 * Tap4_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + thr_v4  * Tap3_v4;
          //dTap_v4 = dTap_v4 * r_ij_v4 + two_v4  * Tap2_v4;
          //dTap_v4 += Tap1_v4 * rijinv_v4;
          Tap_v4 = simd_vmad(Tap7_v4, r_ij_v4, Tap6_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap5_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap4_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap3_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap2_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap1_v4);
          Tap_v4 = simd_vmad(Tap_v4, r_ij_v4, Tap0_v4);//l_pm.wksp_tap0;//Tap0_v4;
          dTap_v4 = simd_vmad(seven_v4, Tap7_v4 * r_ij_v4, six_v4  * Tap6_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, five_v4 * Tap5_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, four_v4 * Tap4_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, thr_v4  * Tap3_v4);
          dTap_v4 = simd_vmad(dTap_v4,  r_ij_v4, two_v4  * Tap2_v4);
          dTap_v4 = simd_vmad(Tap1_v4, rijinv_v4, dTap_v4);
          
          

          //powr_vdW1_v4 = simd_vpowd(r_ij_v4, p_vdW1_v4);
          doublev4 ln_rij_v4 = simd_vlnd(r_ij_v4);
          powr_vdW1_v4 = simd_vexpd(p_vdW1_v4 * ln_rij_v4);

          doublev4 base_vdW = powr_vdW1_v4 + powgi_vdW1_v4;

          //fn13_v4 = simd_vpowd(powr_vdW1_v4 + powgi_vdW1_v4, p_vdW1i_v4);
          doublev4 ln_vdW_v4 = simd_vlnd(base_vdW);
          fn13_v4 =  simd_vexpd(p_vdW1i_v4 * ln_vdW_v4);


          exp2_v4 = simd_vexpd(half_v4 * alpha_v4 * (one_v4 - fn13_v4 * r_vdWinv_v4));
          exp1_v4 = exp2_v4 * exp2_v4;
          evdW_v4 = D_v4 * (exp1_v4 - two_v4 * exp2_v4);
          //data_vdW_v4 += (Tap_v4 * evdW_v4) * remd_v4[jsz & 3];
          data_vdW_v4 = simd_vmad(Tap_v4, evdW_v4 * remd_v4[jsz & 3], data_vdW_v4);

          //dfn13_v4 = simd_vpowd(powr_vdW1_v4 + powgi_vdW1_v4, p_vdW1i_v4-one_v4) 
          //            * simd_vpowd(r_ij_v4, p_vdW1_v4 - two_v4);
          dfn13_v4 = simd_vexpd(power_onev4 *ln_vdW_v4)* simd_vexpd(power_twov4 * ln_rij_v4);


          CEvd_v4 = dTap_v4 * evdW_v4 - Tap_v4 * D_v4 * alpha_div_rvdW_v4 
                    * (exp1_v4-exp2_v4) *dfn13_v4;
          //dr3gamij_1_v4 = (r_ij_v4 * r_ij_v4 * r_ij_v4 + gamma_v4);
          dr3gamij_1_v4 = simd_vmad(r_ij_v4, r_ij_v4 * r_ij_v4, gamma_v4);
          dr3gamij_3_v4 = simd_vpowd(dr3gamij_1_v4, thrinv_v4);
          //dr3gamij_1inv_v4 = one_v4 / dr3gamij_1_v4;
          //dr3gamij_3inv_v4 = one_v4 / dr3gamij_3_v4;
          dr3gamij_1inv_v4 = simd_vdivd(one_v4, dr3gamij_1_v4);
          dr3gamij_3inv_v4 = simd_vdivd(one_v4, dr3gamij_3_v4);

          
          tmp_v4 = Tap_v4 * dr3gamij_3inv_v4;
          doublev4 Cqij_v4 = Cele_v4 * qi_v4 * qj_v4;
          data_ele_v4 += ecoul_v4 = Cqij_v4 * tmp_v4 * remd_v4[jsz & 3];
          
          CEclmb_v4 = Cqij_v4 * (dTap_v4 -  Tap_v4 * r_ij_v4 * dr3gamij_1inv_v4) 
                     * dr3gamij_3inv_v4;
            
          evdwl_v4 = Tap_v4 * evdW_v4;
          
          fpair_v4 = -(CEvd_v4 + CEclmb_v4)* remd_v4[jsz & 3];
          //tmpc_v4 = (CEvd_v4 + CEclmb_v4)*remd_v4[jsz & 3];
          
          //eng_vdwl_v4 += (half_v4 *evdwl_v4)*remd_v4[jsz & 3];
          //eng_coul_v4 += (half_v4 *ecoul_v4)*remd_v4[jsz & 3];
          #if(EVFLAG)
          eng_vdwl_v4 = simd_vmad(half_v4, evdwl_v4*remd_v4[jsz & 3], eng_vdwl_v4);
          eng_coul_v4 = simd_vmad(half_v4, ecoul_v4*remd_v4[jsz & 3], eng_coul_v4);
          #endif


          //virial0_v4 += (half_v4 * del0_v4 * del0_v4 * fpair_v4);// * remd_v4[jsz & 3];
          //virial1_v4 += (half_v4 * del1_v4 * del1_v4 * fpair_v4);// * remd_v4[jsz & 3];
          //virial2_v4 += (half_v4 * del2_v4 * del2_v4 * fpair_v4);// * remd_v4[jsz & 3];
          //virial3_v4 += (half_v4 * del0_v4 * del1_v4 * fpair_v4);// * remd_v4[jsz & 3];
          //virial4_v4 += (half_v4 * del0_v4 * del2_v4 * fpair_v4);// * remd_v4[jsz & 3];
          //virial5_v4 += (half_v4 * del1_v4 * del2_v4 * fpair_v4);// * remd_v4[jsz & 3];
          doublev4 half_fpair = half_v4 * fpair_v4;
          virial0_v4 = simd_vmad(half_fpair, del0_v4 * del0_v4, virial0_v4);
          virial1_v4 = simd_vmad(half_fpair, del1_v4 * del1_v4, virial1_v4);
          virial2_v4 = simd_vmad(half_fpair, del2_v4 * del2_v4, virial2_v4);
          virial3_v4 = simd_vmad(half_fpair, del0_v4 * del1_v4, virial3_v4);
          virial4_v4 = simd_vmad(half_fpair, del0_v4 * del2_v4, virial4_v4);
          virial5_v4 = simd_vmad(half_fpair, del1_v4 * del2_v4, virial5_v4);


          //fi0_v4 += (-tmpc_v4 * del0_v4);
          //fi1_v4 += (-tmpc_v4 * del1_v4);
          //fi2_v4 += (-tmpc_v4 * del2_v4);
          fi0_v4 = simd_vmad(fpair_v4, del0_v4, fi0_v4);
          fi1_v4 = simd_vmad(fpair_v4, del1_v4, fi1_v4);
          fi2_v4 = simd_vmad(fpair_v4, del2_v4, fi2_v4);

        }
        //lwpf_stop(REMAINDER);
        //lwpf_stop(PJ);
      }//for-jst 
      
      simd_vsumd(fi0_v4);
      simd_vsumd(fi1_v4);
      simd_vsumd(fi2_v4);

      fi[ioff][0] += fi0_v4;
      fi[ioff][1] += fi1_v4;
      fi[ioff][2] += fi2_v4;

      virial[0] += packed_atoms[ioff].x[0] * fi[ioff][0];
      virial[1] += packed_atoms[ioff].x[1] * fi[ioff][1];
      virial[2] += packed_atoms[ioff].x[2] * fi[ioff][2];
      virial[3] += packed_atoms[ioff].x[0] * fi[ioff][1];
      virial[4] += packed_atoms[ioff].x[0] * fi[ioff][2];
      virial[5] += packed_atoms[ioff].x[1] * fi[ioff][2];

    }//for-ii
    pe_put(l_pm.f[ist], fi[0], sizeof(double)*isz*4);
    dma_syn();
  }//for-ist
  //lwpf_stop(CMP);

  //lwpf_start(SUM);
  //sum
  simd_vsumd(data_vdW_v4);
  simd_vsumd(data_ele_v4);
  #if(EVFLAG)
  simd_vsumd(eng_vdwl_v4);
  simd_vsumd(eng_coul_v4);
  #endif
  simd_vsumd(virial0_v4);
  simd_vsumd(virial1_v4);
  simd_vsumd(virial2_v4);
  simd_vsumd(virial3_v4);
  simd_vsumd(virial4_v4);
  simd_vsumd(virial5_v4);
  //reduce
  *data_vdW += data_vdW_v4;
  *data_ele += data_ele_v4;
  #if(EVFLAG)
  *eng_vdwl += eng_vdwl_v4;
  *eng_coul += eng_coul_v4;
  #endif

  virial[0] += virial0_v4;
  virial[1] += virial1_v4;
  virial[2] += virial2_v4;
  virial[3] += virial3_v4;
  virial[4] += virial4_v4;
  virial[5] += virial5_v4;

  //lwpf_stop(SUM);
  //lwpf_start(REDUCE);
  reg_reduce_inplace_doublev4(eng_virial, 4);

  if(_MYID == 0)
  {
    pe_put(l_pm.packed_eng, eng_virial, sizeof(double)*16);
    dma_syn();
  }
  //lwpf_stop(REDUCE);
  //lwpf_stop(ALL);
  //lwpf_exit(NONBOND);
}

