#if (EVFLAG==0)
void Merge_Torsion_Valence_Angles_ev0_cpe(merge_param_pack_t *param)
#else
void Merge_Torsion_Valence_Angles_ev1_cpe(merge_param_pack_t *param)
#endif
{
  //lwpf_enter(TORSION);
  //lwpf_start(ALL);

  dma_init();
  merge_param_pack_t l_pm;
  pe_get(param, &l_pm, sizeof(merge_param_pack_t));
  dma_syn();

  
  reax_system_c *system         = l_pm.system;
  control_params *control       = l_pm.control;
  simulation_data *data         = l_pm.data;
  storage *workspace            = l_pm.workspace;
  reax_list **lists             = l_pm.lists;
  atom_pack_t *packed_atoms     = l_pm.packed_atoms;
 
  int i, j, k, l, pi, pj, pk, pl, pij, plk, natoms, pw, w, h, ph, t_tmp, t;
  int type_i, type_j, type_k, type_l, type_w, type_h;
  int start_j, end_j, len_j, start_k, end_k, len_k;
  int start_pj, end_pj, start_pk, end_pk;
  int toff, pk_off, ph_off, pw_off;

  double Delta_j, Delta_k;
  double r_ij, r_jk, r_kl, r_li;
  double BOA_ij, BOA_jk, BOA_kl, BOA_kw, BOA_hj;

  double exp_tor2_ij, exp_tor2_jk, exp_tor2_kl;
  double exp_tor1, exp_tor3_DjDk, exp_tor4_DjDk, exp_tor34_inv;
  double exp_cot2_jk, exp_cot2_ij, exp_cot2_kl;
  double fn10, f11_DjDk, dfn11, fn12;
  double theta_ijk, theta_jkl, theta_jkw, theta_hjk;
  double cos_theta_jkw, cos_theta_hjk;
  double sin_ijk, sin_jkl;
  double cos_ijk, cos_jkl;
  double tan_ijk_i, tan_jkl_i;
  double omega, cos_omega, cos2omega, cos3omega;
  rvec dcos_omega_di, dcos_omega_dj, dcos_omega_dk, dcos_omega_dl;
  double CV, cmn, CEtors1, CEtors2, CEtors3, CEtors4;
  double CEtors5, CEtors6, CEtors7, CEtors8, CEtors9;
  double Cconj, CEconj1, CEconj2, CEconj3;
  double CEconj4, CEconj5, CEconj6;
  double e_tor, e_con;
  rvec dvec_li;
  rvec force, ext_press;
  ivec rel_box_jl;
  four_body_header *fbh;
  four_body_parameters *fbp;
  bond_data *pbond_ij, *pbond_jk, *pbond_kl, *pbond_kw, *pbond_hj, *pbond_pj, *pbond_jt;
  three_body_interaction_data *p_ijk, *p_jkl, *p_jkw, *p_tmp, *p_hjk;
  three_body_interaction_data thb_intrs_data;
  three_body_interaction_data thb_intrs_data_hjk;
  
  
    
  bond_data *bond_list  = l_pm.bond_list;
  double *Cdbo_list     = l_pm.Cdbo_list;
  double *Cdbopi_list   = l_pm.Cdbopi_list;
  double *Cdbopi2_list  = l_pm.Cdbopi2_list;
  double *BO_list       = l_pm.BO_list;
  rvec2  *BOpi_list     = l_pm.BOpi_list;
  rvec4  *fCdDelta      = l_pm.fCdDelta;
  rvec2  *bo_dboc       = l_pm.bo_dboc;

  // Virial tallying variables
  double eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];

  //valence
  int cnt;
  double temp, temp_bo_jt, pBOjt7;
  double p_val1, p_val2, p_val3, p_val4, p_val5;
  double p_val6, p_val7, p_val8, p_val9, p_val10;
  double p_pen1, p_pen2, p_pen3, p_pen4;
  double p_coa1, p_coa2, p_coa3, p_coa4;
  double trm8, expval6, expval7, expval2theta, expval12theta, exp3hj, exp3jk;
  double exp_pen2hj, exp_pen2jk, exp_pen3, exp_pen4, trm_pen34, exp_coa2;
  double dSBO1, dSBO2, SBO, SBO2, CSBO2, SBOp, prod_SBO, vlpadj;
  double CEval1, CEval2, CEval3, CEval4, CEval5, CEval6, CEval7, CEval8;
  double CEpen1, CEpen2, CEpen3;
  double e_ang, e_coa, e_pen;
  double CEcoa1, CEcoa2, CEcoa3, CEcoa4, CEcoa5;
  double Cf7hj, Cf7jk, Cf8j, Cf9j;
  double f7_hj, f7_jk, f8_Dj, f9_Dj;
  double Ctheta_0, theta_0, theta_00, theta, cos_theta, sin_theta_hjk;
  three_body_header *thbh;
  three_body_parameters *thbp;
  double p_tor2, p_tor3, p_tor4, p_cot2;
  
  p_tor2  = l_pm.p_tor2 ;
  p_tor3  = l_pm.p_tor3 ;
  p_tor4  = l_pm.p_tor4 ;
  p_cot2  = l_pm.p_cot2 ;
  p_val6  = l_pm.p_val6 ;
  p_val8  = l_pm.p_val8 ;
  p_val9  = l_pm.p_val9 ;
  p_val10 = l_pm.p_val10;
  p_pen2  = l_pm.p_pen2 ;
  p_pen3  = l_pm.p_pen3 ;
  p_pen4  = l_pm.p_pen4 ;
  p_coa2  = l_pm.p_coa2 ;
  p_coa3  = l_pm.p_coa3 ;
  p_coa4  = l_pm.p_coa4 ;
  int n             = l_pm.n; 
  int ntypes        = l_pm.ntypes; 
  double thb_cut    = l_pm.thb_cut;
  double thb_cutsq  = l_pm.thb_cutsq;

  int jst, jed, jsz, joff;
 
  doublev4 eng_virial[2];
  eng_virial[0] = eng_virial[1] = 0;
  double *eng_vdwl = (double*)(void*)eng_virial;
  double *eng_ang  = eng_vdwl + 1;
  double *eng_pen  = eng_vdwl + 2;
  double *eng_coa  = eng_vdwl + 3;
  double *eng_tor  = eng_vdwl + 4;
  double *eng_con  = eng_vdwl + 5;

  swcache_lock_t *locks     = l_pm.locks;
  swcache_lock_t *locks_frc = l_pm.locks_frc;
  SWCACHE_INIT_U(DT_CACHE, locks);
  SWCACHE_INIT_U(DT_CACHE_F, locks_frc);

  //sizeof(atom_pack_t) = 40;
  int index[JSTEP], end_index[JSTEP];
  atom_pack_t atom_j[JSTEP];
  rvec2 bo_dboc_j[JSTEP];
  double Delta[JSTEP], Delta_val[JSTEP], dDelta_lp[JSTEP];
  double nlp[JSTEP],  vlpex[JSTEP];
 
  //sizeof(bond_data) = 56Bytes;
  double Cdbopi_j[BLEN], Cdbopi2_j[BLEN];
  double BO_list_j[BLEN], BO_list_k[BLEN];
  rvec2 BOpi_list_j[BLEN];
  bond_data blist_j[BLEN], blist_k[BLEN]; 

  //read_cache;
  int read_ctag[READ_C_LCNT];
  atom_pack_t atoms_cache[READ_C_LCNT][READ_C_LSZ];
  atom_pack_t atom_k, atom_l, atom_h;//, atom_i;
  for(i = 0; i < READ_C_LCNT; i++)
  {
    read_ctag[i] = -1;
  }

  //read total_bond_order and Delta_boc for posion k and h;
  int bo_dboc_tag[DT_C_LCNT];
  rvec2 bo_dboc_cache[DT_C_LCNT][DT_C_LSZ];
  rvec2 bo_dboc_k, bo_dboc_h;
  //for(i = 0; i < DT_C_LCNT; i++)
  //{
  //  bo_dboc_tag[i] = -1;
  //}

  //read index for posion k;
  int k_tag[K_C_LCNT];
  int k_idx_cache[   K_C_LCNT][K_C_LSZ];
  int k_ed_idx_cache[K_C_LCNT][K_C_LSZ];

  //for(i = 0; i < K_C_LCNT; i++)
  //{
  //  k_tag[i] = -1;
  //}

  //read bond_data for pj;
  int bond_tag[BOND_C_LCNT];
  bond_data bond_cache[BOND_C_LCNT][BOND_C_LSZ];
  bond_data bond_pj;
  for(i = 0; i < BOND_C_LCNT; i++)
  {
    bond_tag[i] = -1;
    bo_dboc_tag[i] = -1;
    k_tag[i] = -1;
  }

  rvec4 fdelk, fdelh, fdeli, fdell;

  int thbp_cnt[128], thbp_idx;
  three_body_parameters thbp_prm[100];
  four_body_parameters fbp_prm[256];
  int fbp_cnt[256], fbp_idx;
  pe_get(l_pm.thbp_cnt, thbp_cnt, sizeof(int) * 128);
  pe_get(l_pm.thbp_prm, thbp_prm, sizeof(three_body_parameters)*l_pm.thbp_tot);
  pe_get(l_pm.fbp_cnt,  fbp_cnt,  sizeof(int)*l_pm.fbp_tot);
  pe_get(l_pm.fbp_prm,  fbp_prm,  sizeof(four_body_parameters)*l_pm.fbp_tot);
  dma_syn();

  for( jst = _MYID * JSTEP; jst < n; jst+=JSTEP*64)
  {
    jed = jst + JSTEP;
    if(jed > n)
      jed = n;
    jsz = jed - jst;

    pe_get(l_pm.packed_atoms+jst, atom_j,     sizeof(atom_pack_t) * jsz);
    pe_get(l_pm.index+jst,        index,      sizeof(int) * jsz);
    pe_get(l_pm.end_index+jst,    end_index,  sizeof(int) * jsz);
    pe_get(l_pm.vlpex+jst,        vlpex,      sizeof(double) * jsz);
    pe_get(l_pm.nlp+jst,          nlp,        sizeof(double) * jsz);
    pe_get(l_pm.dDelta_lp+jst,    dDelta_lp,  sizeof(double) * jsz);
    pe_get(l_pm.Delta+jst,        Delta,      sizeof(double) * jsz);
    pe_get(l_pm.Delta_val+jst,    Delta_val,  sizeof(double) * jsz);
    pe_get(l_pm.bo_dboc+jst,      bo_dboc_j,  sizeof(rvec2) * jsz);
    dma_syn();

    for(j = jst; j < jed; j++)
    {
      joff = j - jst;
      type_j  = atom_j[joff].type;
      Delta_j = bo_dboc_j[joff][1];
      start_j = index[joff];
      end_j   = end_index[joff];
      len_j   = end_j - start_j;

      if (type_j < 0 || end_j <= start_j) continue;
      p_val3 = l_pm.sbp[ type_j ].p_val3;
      p_val5 = l_pm.sbp[ type_j ].p_val5;

      pe_get(BO_list + start_j,       BO_list_j,  sizeof(double) * len_j);
      pe_get(BOpi_list + start_j,     BOpi_list_j,sizeof(rvec2) * len_j);
      pe_get(Cdbopi_list + start_j,   Cdbopi_j,   sizeof(double) * len_j);
      pe_get(Cdbopi2_list + start_j,  Cdbopi2_j,  sizeof(double) * len_j);
      pe_get(bond_list + start_j,     blist_j,    sizeof(bond_data) * len_j);
      dma_syn();
      
      rvec4 fdelj;
      fdelj[0] = fdelj[1] = fdelj[2] = fdelj[3] = 0;

      SBOp = 0, prod_SBO = 1;
      for( t = start_j; t < end_j; ++t ) 
      {
        toff = t - start_j;
        SBOp += (BOpi_list_j[toff][0] + BOpi_list_j[toff][1]);
        temp = SQR(BO_list_j[toff]);
        temp *= temp;
        temp *= temp;
        prod_SBO *= p_expd( -temp );
      }

      if (vlpex[joff] >= 0) 
      {
        vlpadj = 0;
        dSBO2 = prod_SBO - 1;
      } 
      else 
      {
        vlpadj = nlp[joff];
        dSBO2 = (prod_SBO - 1) * (1 - p_val8 * dDelta_lp[joff]);
      }

      SBO = SBOp + (1 - prod_SBO) * (-bo_dboc_j[joff][1] - p_val8 * vlpadj);
      dSBO1 = -8 * prod_SBO * (bo_dboc_j[joff][1] + p_val8 * vlpadj );

      if (SBO <= 0)
        SBO2 = 0, CSBO2 = 0;
      else if (SBO > 0 && SBO <= 1) 
      {
          SBO2 = p_powd( SBO, p_val9 );
          CSBO2 = p_val9 * p_powd( SBO, p_val9 - 1 );
      }
      else if (SBO > 1 && SBO < 2) 
      {
        SBO2 = 2 - p_powd( 2-SBO, p_val9 );
        CSBO2 = p_val9 * p_powd( 2 - SBO, p_val9 - 1 );
      }
      else
        SBO2 = 2, CSBO2 = 0;

      expval6 = p_expd( p_val6 * bo_dboc_j[joff][1]);
      
      for( pk = start_j; pk < end_j; ++pk ) 
      {
        pk_off    = pk - start_j;
        pbond_jk  = &(blist_j[pk_off]);
        k         = pbond_jk->nbr;
        BOA_jk    = BO_list_j[pk_off] - thb_cut;

         
        if(BOA_jk > 0.0) 
        {//if(BOA_jk > 0.0 && (j < system->n || pbond_jk->nbr < system->n))
          pj = pbond_jk->sym_index; // pj points to j on k's list
          read_cache_bond(pj, bond_cache, bond_tag, bond_list, &bond_pj);
          read_cache_k_index(k, k_idx_cache, k_ed_idx_cache, k_tag, 
                            l_pm.index, l_pm.end_index, 
                            &start_k, &end_k);
          pbond_pj = &bond_pj;
          len_k = end_k - start_k;


          if(end_k > start_k)//if (Num_Entries(pj, thb_intrs)) 
          {
            fdelk[0] = fdelk[1] = fdelk[2] = fdelk[3] = 0;

            pe_get(BO_list + start_k, BO_list_k, sizeof(double) * len_k);
            pe_get(bond_list + start_k, blist_k, sizeof(bond_data) * len_k);
            dma_syn();

            read_cache(k, atoms_cache, read_ctag, packed_atoms, &atom_k);
            read_cache_bo_dboc(k, bo_dboc_cache, bo_dboc_tag, bo_dboc, bo_dboc_k);
            type_k  = atom_k.type;
            Delta_k = bo_dboc_k[1];
            r_jk    = pbond_jk->d;

            double fck = 0;

            exp_tor2_jk = p_expd( -p_tor2 * BOA_jk );
            exp_cot2_jk = p_expd( -p_cot2 * SQR(BOA_jk - 1.5) );
            exp_tor3_DjDk = p_expd( -p_tor3 * (Delta_j + Delta_k) );
            exp_tor4_DjDk = p_expd( p_tor4  * (Delta_j + Delta_k) );
            exp_tor34_inv = 1.0 / (1.0 + exp_tor3_DjDk + exp_tor4_DjDk);
            f11_DjDk = (2.0 + exp_tor3_DjDk) * exp_tor34_inv;

            for(ph = start_j; ph < end_j; ++ph) 
            {
              ph_off    = ph - start_j;
              pbond_hj  = &(blist_j[ph_off]);
              BOA_hj    = BO_list_j[ph_off] - thb_cut;
              h         = pbond_hj->nbr;
              ph_off    = ph - start_j;

              if(ph != pk)
              {
                double fch = 0;
                fdelh[0] = fdelh[1] = fdelh[2] = fdelh[3] = 0;

                read_cache(h, atoms_cache, read_ctag, packed_atoms, &atom_h);
                read_cache_bo_dboc(h, bo_dboc_cache, bo_dboc_tag, bo_dboc, bo_dboc_h);
                type_h = atom_h.type;

                p_hjk   = &thb_intrs_data_hjk;
                
                Calculate_Theta( pbond_jk->dvec, pbond_jk->d,
                                 pbond_hj->dvec, pbond_hj->d,
                                 &theta_hjk, &cos_theta_hjk);

                Calculate_dCos_Theta( pbond_jk->dvec, pbond_jk->d,
                                      pbond_hj->dvec, pbond_hj->d,
                                      &(p_hjk->dcos_di), &(p_hjk->dcos_dj),
                                      &(p_hjk->dcos_dk));
                    
                p_hjk->thb    = h;
                p_hjk->pthb   = ph;
                p_hjk->theta  = theta_hjk;
                
                sin_theta_hjk = p_sind(theta_hjk);
                if (sin_theta_hjk < 1.0e-5)
                  sin_theta_hjk = 1.0e-5;

                i         = h;
                type_i    = type_h;
                theta_hjk = p_hjk->theta;
                pbond_ij  = &(blist_j[ph_off]);
                p_ijk     = p_hjk;

                //Valence
                if((ph > pk) && (j < n) &&
                  (BO_list_j[pk_off] > thb_cut) &&
                  (BO_list_j[ph_off] > thb_cut) &&
                  (BO_list_j[pk_off] * BO_list_j[ph_off] > thb_cutsq)) 
                {

                  thbp_idx = (type_k * ntypes+type_j)*ntypes+type_h;
                  int pos = thbp_cnt[thbp_idx];

                  for( cnt = 0; cnt < thbp_cnt[NT3 + thbp_idx]; ++cnt ) 
                  {
                    thbp = &(thbp_prm[thbp_cnt[thbp_idx]+cnt]);
                    
                    if (fabs(thbp->p_val1) > 0.001) 
                    {
                      /* ANGLE ENERGY */
                      p_val1 = thbp->p_val1;
                      p_val2 = thbp->p_val2;
                      p_val4 = thbp->p_val4;
                      p_val7 = thbp->p_val7;
                      theta_00 = thbp->theta_00;

                      exp3jk = p_expd( -p_val3 * p_powd( BOA_jk, p_val4 ) );
                      f7_jk = 1.0 - exp3jk;
                      Cf7jk = p_val3 * p_val4 * p_powd( BOA_jk, p_val4 - 1.0 ) * exp3jk;

                      exp3hj = p_expd( -p_val3 * p_powd( BOA_hj, p_val4 ) );
                      f7_hj = 1.0 - exp3hj;
                      Cf7hj = p_val3 * p_val4 * p_powd( BOA_hj, p_val4 - 1.0 ) * exp3hj;

                      expval7 = p_expd( -p_val7 * bo_dboc_j[joff][1]);
                      trm8 = 1.0 + expval6 + expval7;
                      f8_Dj = p_val5 - ( (p_val5 - 1.0) * (2.0 + expval6) / trm8 );
                      Cf8j = ( (1.0 - p_val5) / SQR(trm8) ) *
                        ( p_val6 * expval6 * trm8 -
                          (2.0 + expval6) * ( p_val6*expval6 - p_val7*expval7 ) );

                      theta_0 = 180.0 - theta_00 * (1.0 - p_expd(-p_val10 * (2.0 - SBO2)));
                      theta_0 = DEG2RAD( theta_0 );

                      expval2theta  = p_expd( -p_val2 * SQR(theta_0 - theta_hjk) );
                      if (p_val1 >= 0)
                        expval12theta = p_val1 * (1.0 - expval2theta);
                      else // To avoid linear Me-H-Me angles (6/6/06)
                        expval12theta = p_val1 * -expval2theta;

                      CEval1 = Cf7jk * f7_hj * f8_Dj * expval12theta;
                      CEval2 = Cf7hj * f7_jk * f8_Dj * expval12theta;
                      CEval3 = Cf8j  * f7_hj * f7_jk * expval12theta;
                      CEval4 = -2.0 * p_val1 * p_val2 * f7_jk * f7_hj * f8_Dj *
                        expval2theta * (theta_0 - theta_hjk);

                      Ctheta_0 = p_val10 * DEG2RAD(theta_00) *
                        p_expd( -p_val10 * (2.0 - SBO2) );

                      CEval5 = -CEval4 * Ctheta_0 * CSBO2;
                      CEval6 = CEval5 * dSBO1;
                      CEval7 = CEval5 * dSBO2;
                      CEval8 = -CEval4 / sin_theta_hjk;

                      *eng_ang += e_ang = f7_jk * f7_hj * f8_Dj * expval12theta;
                      /* END ANGLE ENERGY*/

                      /* PENALTY ENERGY */
                      p_pen1 = thbp->p_pen1;
                      
                      exp_pen2jk = p_expd( -p_pen2 * SQR( BOA_jk - 2.0 ) );
                      exp_pen2hj = p_expd( -p_pen2 * SQR( BOA_hj - 2.0 ) );
                      exp_pen3 = p_expd( -p_pen3 * Delta[joff] );
                      exp_pen4 = p_expd(  p_pen4 * Delta[joff] );
                      trm_pen34 = 1.0 + exp_pen3 + exp_pen4;
                      f9_Dj = ( 2.0 + exp_pen3 ) / trm_pen34;
                      Cf9j = ( -p_pen3 * exp_pen3 * trm_pen34 -
                               (2.0 + exp_pen3) * ( -p_pen3 * exp_pen3 +
                                                    p_pen4 * exp_pen4 ) ) / SQR( trm_pen34 );

                      *eng_pen += e_pen = p_pen1 * f9_Dj * exp_pen2jk * exp_pen2hj;

                      CEpen1 = e_pen * Cf9j / f9_Dj;
                      temp   = -2.0 * p_pen2 * e_pen;
                      CEpen2 = temp * (BOA_jk - 2.0);
                      CEpen3 = temp * (BOA_hj - 2.0);
                      /* END PENALTY ENERGY */

                      /* COALITION ENERGY */
                      p_coa1 = thbp->p_coa1;

                      exp_coa2 = p_expd( p_coa2 * Delta_val[joff]);
                      *eng_coa += e_coa =
                        p_coa1 / (1. + exp_coa2) *
                        p_expd( -p_coa3 * SQR(bo_dboc_k[0]-BOA_jk) ) *
                        p_expd( -p_coa3 * SQR(bo_dboc_h[0]-BOA_hj) ) *
                        p_expd( -p_coa4 * SQR(BOA_jk - 1.5) ) *
                        p_expd( -p_coa4 * SQR(BOA_hj - 1.5) );

                      CEcoa1 = -2 * p_coa4 * (BOA_jk - 1.5) * e_coa;
                      CEcoa2 = -2 * p_coa4 * (BOA_hj - 1.5) * e_coa;
                      CEcoa3 = -p_coa2 * exp_coa2 * e_coa / (1 + exp_coa2);
                      CEcoa4 = -2 * p_coa3 * (bo_dboc_k[0]-BOA_jk) * e_coa;
                      CEcoa5 = -2 * p_coa3 * (bo_dboc_h[0]-BOA_hj) * e_coa;
                      /* END COALITION ENERGY */

                      /* FORCES */
                      fck += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));
                      fch += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));

                      fdelj[3] += ((CEval3 + CEval7) + CEpen1 + CEcoa3);//fCdDelta[j][3]
                      fdelk[3] += CEcoa4;//fCdDelta[k][3]
                      fdelh[3] += CEcoa5;//fCdDelta[h][3]

                      for( t = start_j; t < end_j; ++t ) 
                      {
                        toff = t - start_j;
                        temp_bo_jt = BO_list_j[toff];
                        temp = CUBE( temp_bo_jt );
                        pBOjt7 = temp * temp * temp_bo_jt;

                        double fct = (CEval6 * pBOjt7);
                        SWCACHE_UPDATE(DT_CACHE, (t), (&fct));

                        Cdbopi_j[toff] += CEval5;
                        Cdbopi2_j[toff] += CEval5;
                      }
                      rvec4_ScaledAdd( fdelj, CEval8, p_hjk->dcos_dj );
                      rvec4_ScaledAdd( fdelk, CEval8, p_hjk->dcos_di );
                      rvec4_ScaledAdd( fdelh, CEval8, p_hjk->dcos_dk );
                      
                      /* tally into per-atom virials */
                      #if(EVFLAG)
                      eng_tmp = e_ang + e_pen + e_coa;
                      *eng_vdwl += eng_tmp;
                      #endif
                    }
                  }//for-cnt
                  //SWCACHE_UPDATE(DT_CACHE_F, h, (&fdelh));
                }//if
                
                //Torsion
                int flag_c = 0;
                if (atom_j[joff].orig_id >  atom_k.orig_id)
                {
                  flag_c = 1;
                }
                if (atom_j[joff].orig_id == atom_k.orig_id) 
                {
                  if (atom_k.x[2] <  atom_j[joff].x[2]) 
                  {
                    flag_c = 1;
                  }
                  if (atom_k.x[2] == atom_j[joff].x[2] &&
                      atom_k.x[1] <  atom_j[joff].x[1]) 
                  {
                    flag_c = 1;
                  }
                  if (atom_k.x[2] == atom_j[joff].x[2] &&
                      atom_k.x[1] == atom_j[joff].x[1] &&
                      atom_k.x[0] <  atom_j[joff].x[0]) 
                  {
                    flag_c = 1;
                  }
                }
                if(flag_c)
                {
                  SWCACHE_UPDATE(DT_CACHE, (ph), (&fch));
                  SWCACHE_UPDATE(DT_CACHE_F, j, (&fdelj));
                  SWCACHE_UPDATE(DT_CACHE_F, k, (&fdelk));
                  SWCACHE_UPDATE(DT_CACHE_F, h, (&fdelh));
                  continue;
                }


                if (BO_list_j[ph_off] > thb_cut /*0*/ && j < n && end_k > start_k) 
                {
                  i       = p_ijk->thb;
                  type_i  = atom_h.type;
                  r_ij    = pbond_ij->d;
                  BOA_ij  = BO_list_j[ph_off] - thb_cut;

                  theta_ijk = p_ijk->theta;
                  sin_ijk   = p_sind( theta_ijk );
                  cos_ijk   = p_cosd( theta_ijk );
                  if (sin_ijk >= 0 && sin_ijk <= MIN_SINE)
                    tan_ijk_i = cos_ijk / MIN_SINE;
                  else if( sin_ijk <= 0 && sin_ijk >= -MIN_SINE )
                    tan_ijk_i = cos_ijk / -MIN_SINE;
                  else tan_ijk_i = cos_ijk / sin_ijk;

                  exp_tor2_ij = p_expd( -p_tor2 * BOA_ij );
                  exp_cot2_ij = p_expd( -p_cot2 * SQR(BOA_ij -1.5) );
                  
                  fdeli[0] = fdeli[1] = fdeli[2] = fdeli[3] = 0;

                  for(pw = start_k; pw < end_k; pw++)
                  {
                    pw_off    = pw - start_k;
                    pbond_kw  = &(blist_k[pw_off]);
                    BOA_kw    = BO_list_k[pw_off] - thb_cut;
                    w         = pbond_kw->nbr;
                    if(pw != pj)
                    {
                      read_cache(w, atoms_cache, read_ctag, packed_atoms, &atom_l);
                      type_w  = atom_l.type;
                      p_jkw   = &thb_intrs_data;
                      
                      Calculate_Theta( pbond_pj->dvec, pbond_pj->d,
                                       pbond_kw->dvec, pbond_kw->d,
                                       &theta_jkw, &cos_theta_jkw);

                      Calculate_dCos_Theta( pbond_pj->dvec, pbond_pj->d,
                                            pbond_kw->dvec, pbond_kw->d,
                                            &(p_jkw->dcos_di), &(p_jkw->dcos_dj),
                                            &(p_jkw->dcos_dk));
                          
                      p_jkw->thb    = w;
                      p_jkw->pthb   = pw;
                      p_jkw->theta  = theta_jkw;

                      l         = w;
                      type_l    = type_w;
                      theta_jkl = p_jkw->theta;
                      pbond_kl  = &(blist_k[pw_off]);
                      p_jkl     = p_jkw;
                      BOA_kl    = BOA_kw;
                      
                      fbp_idx = ((type_i*ntypes+type_j)*ntypes+type_k) * ntypes;
                      fbp = &(fbp_prm[fbp_idx+type_l]);

                      if(i != l && fbp_cnt[fbp_idx+type_l] && BO_list_k[pw_off] > thb_cut &&
                          BO_list_j[ph_off] * BO_list_j[pk_off] * BO_list_k[pw_off] > thb_cut)

                      {
                        r_kl = pbond_kl->d;
                        BOA_kl = BO_list_k[pw_off] - thb_cut;

                        sin_jkl = p_sind( theta_jkl );
                        cos_jkl = p_cosd( theta_jkl );
                        if (sin_jkl >= 0 && sin_jkl <= MIN_SINE)
                          tan_jkl_i = cos_jkl / MIN_SINE;
                        else if( sin_jkl <= 0 && sin_jkl >= -MIN_SINE )
                          tan_jkl_i = cos_jkl / -MIN_SINE;
                        else tan_jkl_i = cos_jkl /sin_jkl;

                        rvec_ScaledSum( dvec_li, 1., atom_h.x, -1., atom_l.x);

                        double tmp_r, tmp_dev3;
                        tmp_dev3 = SQR(dvec_li[0]) + SQR(dvec_li[1]) + SQR(dvec_li[2]);
                        inv_sqrt(tmp_dev3, tmp_r);
                        r_li = tmp_dev3 * tmp_r;

                        /* omega and its derivative */
                        omega = Calculate_Omega_cpe( pbond_ij->dvec, r_ij,
                                                 pbond_jk->dvec, r_jk,
                                                 pbond_kl->dvec, r_kl,
                                                 dvec_li, r_li,
                                                 p_ijk, p_jkl,
                                                 dcos_omega_di, dcos_omega_dj,
                                                 dcos_omega_dk, dcos_omega_dl);

                        cos_omega = p_cosd( omega );
                        cos2omega = p_cosd( 2. * omega );
                        cos3omega = p_cosd( 3. * omega );
                        /* end omega calculations */

                        /* torsion energy */
                        exp_tor1 = p_expd(fbp->p_tor1 * SQR(2.0 - BOpi_list_j[pk_off][0] - f11_DjDk));
                        exp_tor2_kl = p_expd( -p_tor2 * BOA_kl );
                        exp_cot2_kl = p_expd( -p_cot2 * SQR(BOA_kl - 1.5) );
                        fn10 = (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk) *
                          (1.0 - exp_tor2_kl);

                        CV = 0.5 * ( fbp->V1 * (1.0 + cos_omega) +
                                     fbp->V2 * exp_tor1 * (1.0 - cos2omega) +
                                     fbp->V3 * (1.0 + cos3omega) );

                        *eng_tor += e_tor = fn10 * sin_ijk * sin_jkl * CV;

                        dfn11 = (-p_tor3 * exp_tor3_DjDk +
                                 (p_tor3 * exp_tor3_DjDk - p_tor4 * exp_tor4_DjDk) *
                                 (2.0 + exp_tor3_DjDk) * exp_tor34_inv) *
                          exp_tor34_inv;

                        CEtors1 = sin_ijk * sin_jkl * CV;

                        CEtors2 = -fn10 * 2.0 * fbp->p_tor1 * fbp->V2 * exp_tor1 *
                          (2.0 - BOpi_list_j[pk_off][0] -f11_DjDk) * (1.0 -SQR(cos_omega)) *
                          sin_ijk * sin_jkl;

                        CEtors3 = CEtors2 * dfn11;

                        CEtors4 = CEtors1 * p_tor2 * exp_tor2_ij *
                          (1.0 - exp_tor2_jk) * (1.0 - exp_tor2_kl);
                        CEtors5 = CEtors1 * p_tor2 *
                          (1.0 - exp_tor2_ij) * exp_tor2_jk * (1.0 - exp_tor2_kl);
                        CEtors6 = CEtors1 * p_tor2 *
                          (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk) * exp_tor2_kl;

                        cmn = -fn10 * CV;
                        CEtors7 = cmn * sin_jkl * tan_ijk_i;
                        CEtors8 = cmn * sin_ijk * tan_jkl_i;

                        CEtors9 = fn10 * sin_ijk * sin_jkl *
                          (0.5 * fbp->V1 - 2.0 * fbp->V2 * exp_tor1 * cos_omega +
                           1.5 * fbp->V3 * (cos2omega + 2.0 * SQR(cos_omega)));
                        /* end  of torsion energy */

                        /* 4-body conjugation energy */
                        fn12 = exp_cot2_ij * exp_cot2_jk * exp_cot2_kl;
                        *eng_con += e_con = fbp->p_cot1 * fn12 *
                          (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);

                        Cconj = -2.0 * fn12 * fbp->p_cot1 * p_cot2 *
                          (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);

                        CEconj1 = Cconj * (BOA_ij - 1.5e0);
                        CEconj2 = Cconj * (BOA_jk - 1.5e0);
                        CEconj3 = Cconj * (BOA_kl - 1.5e0);

                        CEconj4 = -fbp->p_cot1 * fn12 *
                          (SQR(cos_omega) - 1.0) * sin_jkl * tan_ijk_i;
                        CEconj5 = -fbp->p_cot1 * fn12 *
                          (SQR(cos_omega) - 1.0) * sin_ijk * tan_jkl_i;
                        CEconj6 = 2.0 * fbp->p_cot1 * fn12 *
                          cos_omega * sin_ijk * sin_jkl;
                        /* end 4-body conjugation energy */

                        /* forces */
                        Cdbopi_j[pk_off] += CEtors2;//bo_jk->Cdbopi;
                        
                        fdell[3] = 0;
                        fdelj[3] += CEtors3;//fCdDelta[j][3]
                        fdelk[3] += CEtors3;//fCdDelta[k][3]
                        
                        
                        double fcw = (CEtors6 + CEconj3);
                        SWCACHE_UPDATE(DT_CACHE, (pw), (&fcw));//Cdbo_list[pw];
                        fck += (CEtors5 + CEconj2);//Cdbo_list[pk];
                        fch += (CEtors4 + CEconj1);//Cdbo_list[ph];
                        
                        rvec fi_tmp, fj_tmp, fk_tmp, fl_tmp;
                        double CE74 = CEtors7+CEconj4;
                        double CE85 = CEtors8+CEconj5;
                        double CE96 = CEtors9+CEconj6;

                        /* dcos_theta_ijk */
                        rvec_Scale(fi_tmp, CE74, p_ijk->dcos_dk );
                        rvec_Scale(fj_tmp, CE74, p_ijk->dcos_dj );
                        rvec_Scale(fk_tmp, CE74, p_ijk->dcos_di );

                        /* dcos_theta_jkl */
                        rvec_ScaledAdd(fj_tmp, CE85, p_jkl->dcos_di);
                        rvec_ScaledAdd(fk_tmp, CE85, p_jkl->dcos_dj);
                        rvec_Scale    (fl_tmp, CE85, p_jkl->dcos_dk);

                        rvec_ScaledAdd(fi_tmp, CE96, dcos_omega_di );
                        rvec_ScaledAdd(fj_tmp, CE96, dcos_omega_dj );
                        rvec_ScaledAdd(fk_tmp, CE96, dcos_omega_dk );
                        rvec_ScaledAdd(fl_tmp, CE96, dcos_omega_dl );
                        
                        rvec4_Add(fdelj, fj_tmp);
                        rvec4_Add(fdelk, fk_tmp);
                        //rvec4_Add(fdeli, fi_tmp);
                        rvec4_Add(fdelh, fi_tmp);
                        rvec4_Copy(fdell, fl_tmp);
                        SWCACHE_UPDATE(DT_CACHE_F, l, (&fdell));

                        /* tally into per-atom virials */
                        #if(EVFLAG)
                        *eng_vdwl += e_tor + e_con;
                        #endif
                      } // pl check ends
                    }//if pw!=pj
                  }//for-pw
                  //SWCACHE_UPDATE(DT_CACHE_F, i, (&fdeli));

                } // pi check ends

                SWCACHE_UPDATE(DT_CACHE, (ph), (&fch));
                SWCACHE_UPDATE(DT_CACHE_F, h, (&fdelh));
              } // if ph!=pi
            }//for-ph
            SWCACHE_UPDATE(DT_CACHE, (pk), (&fck));
            SWCACHE_UPDATE(DT_CACHE_F, k, (&fdelk));
          } // k-j neighbor check ends
        } // j-k neighbor check ends
      } //for-pk 
      
      pe_put(Cdbopi_list + start_j, Cdbopi_j,   sizeof(double) * len_j);
      pe_put(Cdbopi2_list + start_j, Cdbopi2_j,   sizeof(double) * len_j);
      dma_syn();

      SWCACHE_UPDATE(DT_CACHE_F, j, (&fdelj));

    } // j loop
  }//for-jst
  SWCACHE_FLUSH(DT_CACHE);
  SWCACHE_FLUSH(DT_CACHE_F);

  /******reduce energy*****/
  reg_reduce_inplace_doublev4(eng_virial, 2);
  if(_MYID == 0)
  {
    pe_put(l_pm.packed_eng, eng_virial, sizeof(double) * 8);
    dma_syn();
  }
  //lwpf_stop(ALL);
  //lwpf_exit(TORSION);
}

