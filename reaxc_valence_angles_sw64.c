#include "reaxc_defs_sunway.h"
#include "reaxc_ctypes_sunway.h"
#include <math.h>
#include "reaxc_inlines_sw64.h"
//#include"sunway_math.h"
#include "sleef_math.h"

static double Dot( double* v1, double* v2, int k )
{
  double ret = 0.0;
  int i;
  for( i=0; i < k; ++i )
    ret +=  v1[i] * v2[i];

  return ret;
}

void Calculate_Theta( rvec dvec_ji, double d_ji, rvec dvec_jk, double d_jk,
                      double *theta, double *cos_theta )
{
  (*cos_theta) = Dot( dvec_ji, dvec_jk, 3 ) / ( d_ji * d_jk );
  if( *cos_theta > 1. ) *cos_theta  = 1.0;
  if( *cos_theta < -1. ) *cos_theta  = -1.0;

  //(*theta) = acos( *cos_theta );
  (*theta) = xacos( *cos_theta );
  //(*theta) = acos_sw( *cos_theta );
  //if(myrank == 0)
  //  printf("acos(%.10f) = %.10f, %.10f\n", *cos_theta, *theta, acos(*cos_theta));
  
}

void Calculate_dCos_Theta( rvec dvec_ji, double d_ji, rvec dvec_jk, double d_jk,
                           rvec* dcos_theta_di,
                           rvec* dcos_theta_dj,
                           rvec* dcos_theta_dk )
{
  int t;
  double sqr_d_ji = SQR(d_ji);
  double sqr_d_jk = SQR(d_jk);
  double inv_dists = 1.0 / (d_ji * d_jk);
  double inv_dists3 = pow( inv_dists, 3.0 );
  double dot_dvecs = Dot( dvec_ji, dvec_jk, 3 );
  double Cdot_inv3 = dot_dvecs * inv_dists3;

  for( t = 0; t < 3; ++t ) {
    (*dcos_theta_di)[t] = dvec_jk[t] * inv_dists -
      Cdot_inv3 * sqr_d_jk * dvec_ji[t];
    (*dcos_theta_dj)[t] = -(dvec_jk[t] + dvec_ji[t]) * inv_dists +
      Cdot_inv3 * ( sqr_d_jk * dvec_ji[t] + sqr_d_ji * dvec_jk[t] );
    (*dcos_theta_dk)[t] = dvec_ji[t] * inv_dists -
      Cdot_inv3 * sqr_d_ji * dvec_jk[t];
  }
}

#define EXP_PPC0 1.00000000000969824220931059244e+00
#define EXP_PPC1 9.99999996155887638238368708699e-01
#define EXP_PPC2 4.99999998889176233696218787372e-01
#define EXP_PPC3 1.66666783941218449305310400632e-01
#define EXP_PPC4 4.16666866967536769772451066274e-02
#define EXP_PPC5 8.33238238257071571479794869219e-03
#define EXP_PPC6 1.38876391271098567278818869397e-03
#define EXP_PPC7 2.01234277126192251773997843323e-04
#define EXP_PPC8 2.51169087281361337897402780106e-05
#define EXP_PPC9 0.007812500000000000000000000000000



inline double exp_4_sw(double x){
  double x_128th = x * EXP_PPC9;

  double expx = x_128th * EXP_PPC8 + EXP_PPC7;
  expx = expx * x_128th + EXP_PPC6;
  expx = expx * x_128th + EXP_PPC5;
  expx = expx * x_128th + EXP_PPC4;
  expx = expx * x_128th + EXP_PPC3;
  expx = expx * x_128th + EXP_PPC2;
  expx = expx * x_128th + EXP_PPC1;
  expx = expx * x_128th + EXP_PPC0;
  expx *= expx;
  expx *= expx;
  expx *= expx;
  expx *= expx;
  expx *= expx;
  expx *= expx;
  expx *= expx;
  return expx;
}


#ifdef MPE
#include<mpi.h>
#include<athread.h>
#define THIRD 0.333333333333333333
extern SLAVE_FUN(Valence_Angles_C_para)(param_pack_t *);
//void Valence_Angles_C( reax_system_c *system, control_params *control,
//                       simulation_data *data, storage *workspace,
//                       reax_list **lists)

//void Valence_Angles_C(param_pack_t *param)
//{
//  //return;
//
//  //if(athread_idle() == 0)
//  //  athread_init();
//
//  //athread_spawn(Valence_Angle_C_para, pm);
//  //athread_join();
//  reax_system_c *system = param->system;
//  control_params *control = param->control;
//  simulation_data *data = param->data;
//  storage *workspace = param->workspace;
//  reax_list **lists = param->lists;
//
//  int total = 0;
//
//  int i, j, pi, k, pk, t;
//  int type_i, type_j, type_k;
//  int start_j, end_j, start_pk, end_pk;
//  int cnt, num_thb_intrs;
//
//  double temp, temp_bo_jt, pBOjt7;
//  double p_val1, p_val2, p_val3, p_val4, p_val5;
//  double p_val6, p_val7, p_val8, p_val9, p_val10;
//  double p_pen1, p_pen2, p_pen3, p_pen4;
//  double p_coa1, p_coa2, p_coa3, p_coa4;
//  double trm8, expval6, expval7, expval2theta, expval12theta, exp3ij, exp3jk;
//  double exp_pen2ij, exp_pen2jk, exp_pen3, exp_pen4, trm_pen34, exp_coa2;
//  double dSBO1, dSBO2, SBO, SBO2, CSBO2, SBOp, prod_SBO, vlpadj;
//  double CEval1, CEval2, CEval3, CEval4, CEval5, CEval6, CEval7, CEval8;
//  double CEpen1, CEpen2, CEpen3;
//  double e_ang, e_coa, e_pen;
//  double CEcoa1, CEcoa2, CEcoa3, CEcoa4, CEcoa5;
//  double Cf7ij, Cf7jk, Cf8j, Cf9j;
//  double f7_ij, f7_jk, f8_Dj, f9_Dj;
//  double Ctheta_0, theta_0, theta_00, theta, cos_theta, sin_theta;
//  double BOA_ij, BOA_jk;
//  rvec force, ext_press;
//
//  // Tallying variables
//  double eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];
//  double delij[3], delkj[3];
//
//  three_body_header *thbh;
//  three_body_parameters *thbp;
//  three_body_interaction_data *p_ijk, *p_kji;
//  bond_data *pbond_ij, *pbond_jk, *pbond_jt;
//  bond_order_data *bo_ij, *bo_jk, *bo_jt;
//  reax_list *bonds = (*lists) + BONDS;
//  reax_list *thb_intrs =  (*lists) + THREE_BODIES;
//
//  /* global parameters used in these calculations */
//  p_val6 = system->reax_param.gp.l[14];
//  p_val8 = system->reax_param.gp.l[33];
//  p_val9 = system->reax_param.gp.l[16];
//  p_val10 = system->reax_param.gp.l[17];
//
//  p_pen2 = system->reax_param.gp.l[19];
//  p_pen3 = system->reax_param.gp.l[20];
//  p_pen4 = system->reax_param.gp.l[21];
//
//  p_coa2 = system->reax_param.gp.l[2];
//  p_coa3 = system->reax_param.gp.l[38];
//  p_coa4 = system->reax_param.gp.l[30];
//
//
//  num_thb_intrs = 0;////
//  
//  //printf("in valence---\n");
//  
//  int myrank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//  //double s1 = 0, s2 = 0, s3 = 0, s4 = 0;
//  //double s5 = 0, s6 = 0, s7 = 0, s8 = 0, s9 = 0;
//  //double s11 = 0, s12 = 0, s13 = 0;
//  //if(myrank == 0) 
//  //  printf("N = %d\n", system->N);
//
//  for( j = 0; j < system->N; ++j ) 
//  { // Ray: the first one with system->N
//  
//    //s1 = MPI_Wtime();
//
//    type_j = system->my_atoms[j].type;
//    if (type_j < 0) continue;
//    //start_j = Start_Index(j, bonds);
//    //end_j = End_Index(j, bonds);
//    start_j = bonds->index[j];
//    end_j   = bonds->end_index[j];
//
//
//    p_val3 = system->reax_param.sbp[type_j].p_val3;
//    p_val5 = system->reax_param.sbp[type_j].p_val5;
//
//    SBOp = 0, prod_SBO = 1;
//    
//    //s3 = MPI_Wtime();
//    for( t = start_j; t < end_j; ++t ) 
//    {
//      bo_jt = &(bonds->select.bond_list[t].bo_data);
//      SBOp += (bo_jt->BO_pi + bo_jt->BO_pi2);
//      temp = SQR( bo_jt->BO );
//      temp *= temp;
//      temp *= temp;
//#define SW_MATH
//#ifdef SW_MATH
//      prod_SBO *= exp_4_sw( -temp );
//#else
//      prod_SBO *= exp( -temp );
//#endif
//    }//for
//    //s4 = MPI_Wtime();
//
//    if( workspace->vlpex[j] >= 0 )
//    {
//      vlpadj = 0;
//      dSBO2 = prod_SBO - 1;
//    }
//    else
//    {
//      vlpadj = workspace->nlp[j];
//      dSBO2 = (prod_SBO - 1) * (1 - p_val8 * workspace->dDelta_lp[j]);
//    }
//
//    SBO = SBOp + (1 - prod_SBO) * (-workspace->Delta_boc[j] - p_val8 * vlpadj);
//    dSBO1 = -8 * prod_SBO * ( workspace->Delta_boc[j] + p_val8 * vlpadj );
//
//    if( SBO <= 0 )
//      SBO2 = 0, CSBO2 = 0;
//    else if( SBO > 0 && SBO <= 1 ) 
//    {
//      SBO2 = pow( SBO, p_val9 );
//      CSBO2 = p_val9 * pow( SBO, p_val9 - 1 );
//    }
//    else if( SBO > 1 && SBO < 2 ) 
//    {
//      SBO2 = 2 - pow( 2-SBO, p_val9 );
//      CSBO2 = p_val9 * pow( 2 - SBO, p_val9 - 1 );
//    }
//    else
//      SBO2 = 2, CSBO2 = 0;
//
//#ifdef SW_MATH
//    expval6 = exp_4_sw( p_val6 * workspace->Delta_boc[j] );
//#else
//    expval6 = exp( p_val6 * workspace->Delta_boc[j] );
//#endif
//    
//    int pi_a, pi_b;
//    //s5 = MPI_Wtime();
//    for( pi = start_j; pi < end_j; ++pi ) 
//    {
//      //pi_a = num_thb_intrs;
//
//      //Set_Start_Index( pi, num_thb_intrs, thb_intrs );
//      thb_intrs->index[pi] = num_thb_intrs;//////
//      pbond_ij = &(bonds->select.bond_list[pi]);
//      bo_ij = &(pbond_ij->bo_data);
//      BOA_ij = bo_ij->BO - control->thb_cut;
//      //pi_a = thb_intrs->index[pi];
//
//      if(BOA_ij > 0.0 && (j < system->n || pbond_ij->nbr < system->n)) 
//      {
//        i = pbond_ij->nbr;
//        type_i = system->my_atoms[i].type;///
//
//        //s7 = MPI_Wtime();
//        for( pk = start_j; pk < pi; ++pk ) 
//        {
//          //start_pk = Start_Index( pk, thb_intrs );
//          //end_pk = End_Index( pk, thb_intrs );
//          start_pk = thb_intrs->index[pk];
//          end_pk   = thb_intrs->end_index[pk];
//          //if(myrank == 0)
//          //  printf("val-1: start_j=%d,end_j=%d,start_pk=%d,end_pk=%d\n", start_j,end_j,start_pk,end_pk);
//
//          for(t = start_pk; t < end_pk; ++t)
//          {
//            if(thb_intrs->select.three_body_list[t].thb == i) 
//            {
//              p_ijk = &(thb_intrs->select.three_body_list[num_thb_intrs] );
//              p_kji = &(thb_intrs->select.three_body_list[t]);
//              
//
//              p_ijk->thb = bonds->select.bond_list[pk].nbr;
//              p_ijk->pthb  = pk;
//              p_ijk->theta = p_kji->theta;
//              rvec_Copy( p_ijk->dcos_di, p_kji->dcos_dk );
//              rvec_Copy( p_ijk->dcos_dj, p_kji->dcos_dj );
//              rvec_Copy( p_ijk->dcos_dk, p_kji->dcos_di );
//
//              ++num_thb_intrs;
//              break;
//            }//if
//          }//for-t
//        }//for-pk
//        //s8 = MPI_Wtime();
//        //if(myrank == 0)
//        //    printf("val-1: start_j=%d,end_j=%d,pi+1=%d,end_j=%d\n", start_j,end_j,pi+1,end_j);
//
//        for( pk = pi+1; pk < end_j; ++pk ) 
//        {
//          
//          //s11 = MPI_Wtime();
//
//          pbond_jk = &(bonds->select.bond_list[pk]);
//          bo_jk    = &(pbond_jk->bo_data);
//          BOA_jk   = bo_jk->BO - control->thb_cut;
//          k        = pbond_jk->nbr;
//          //total++;
//          //if(myrank == 0)
//          //{
//          //  printf("(j, i, k) = (%d, %d, %d)\n", j, i, k);
//          //}
//          //if(myrank == 0)
//          //printf("i = %d, j = %d, k = %d\n", i, j, k);
//          type_k   = system->my_atoms[k].type;///
//          p_ijk    = &( thb_intrs->select.three_body_list[num_thb_intrs] );
//          
//
//          Calculate_Theta( pbond_ij->dvec, pbond_ij->d,
//                           pbond_jk->dvec, pbond_jk->d,
//                           &theta, &cos_theta );
//
//          Calculate_dCos_Theta( pbond_ij->dvec, pbond_ij->d,
//                                pbond_jk->dvec, pbond_jk->d,
//                                &(p_ijk->dcos_di), &(p_ijk->dcos_dj),
//                                &(p_ijk->dcos_dk) );
//          p_ijk->thb = k;
//          p_ijk->pthb = pk;
//          p_ijk->theta = theta;
//
//          sin_theta = sin( theta );
//          if( sin_theta < 1.0e-5 )
//            sin_theta = 1.0e-5;
//
//          ++num_thb_intrs;
//          
//          //s13 = MPI_Wtime();
//
//          if((j < system->n) && (BOA_jk > 0.0) && (bo_ij->BO > control->thb_cut) &&
//            (bo_jk->BO>control->thb_cut)&&(bo_ij->BO*bo_jk->BO > control->thb_cutsq)) 
//          {
//            thbh = &( system->reax_param.thbp[type_i][type_j][type_k]);///
//          //if(myrank == 1)
//          //{
//          //  printf("(j, i, k) = (%d, %d, %d)\n", j, i, k);
//          //}
//
//            for( cnt = 0; cnt < thbh->cnt; ++cnt ) 
//            {
//              if( fabs(thbh->prm[cnt].p_val1) > 0.001 ) 
//              {
//                thbp = &(thbh->prm[cnt]);
//
//                /* ANGLE ENERGY */
//                p_val1    = thbp->p_val1;
//                p_val2    = thbp->p_val2;
//                p_val4    = thbp->p_val4;
//                p_val7    = thbp->p_val7;
//                theta_00  = thbp->theta_00;
//
//#ifdef SW_MATH
//                exp3ij = exp_4_sw( -p_val3 * pow( BOA_ij, p_val4 ) );
//#else
//                exp3ij = exp( -p_val3 * pow( BOA_ij, p_val4 ) );
//#endif
//
//                f7_ij = 1.0 - exp3ij;
//                Cf7ij = p_val3 * p_val4 * pow( BOA_ij, p_val4 - 1.0 ) * exp3ij;
//
//#ifdef SW_MATH
//                exp3jk = exp_4_sw( -p_val3 * pow( BOA_jk, p_val4 ) );
//#else
//                exp3jk = exp( -p_val3 * pow( BOA_jk, p_val4 ) );
//#endif
//
//                f7_jk = 1.0 - exp3jk;
//                Cf7jk = p_val3 * p_val4 * pow( BOA_jk, p_val4 - 1.0 ) * exp3jk;
//
//#ifdef SW_MATH
//                expval7 = exp_4_sw( -p_val7 * workspace->Delta_boc[j] );
//#else
//                expval7 = exp( -p_val7 * workspace->Delta_boc[j] );
//#endif
//
//                trm8 = 1.0 + expval6 + expval7;
//                f8_Dj = p_val5 - ( (p_val5 - 1.0) * (2.0 + expval6) / trm8 );
//                Cf8j = ( (1.0 - p_val5) / SQR(trm8) ) *( p_val6 * expval6 * trm8 -
//                    (2.0 + expval6) * ( p_val6*expval6 - p_val7*expval7 ) );
//
//#ifdef SW_MATH
//                theta_0 = 180.0 - theta_00 * (1.0 - exp_4_sw(-p_val10 * (2.0-SBO2)));
//#else
//                theta_0 = 180.0 - theta_00 * (1.0 - exp(-p_val10 * (2.0 - SBO2)));
//#endif
//
//                theta_0 = DEG2RAD( theta_0 );
//
//#ifdef SW_MATH
//                expval2theta  = exp_4_sw( -p_val2 * SQR(theta_0 - theta) );
//#else
//                expval2theta  = exp( -p_val2 * SQR(theta_0 - theta) );
//#endif
//
//                if( p_val1 >= 0 )
//                  expval12theta = p_val1 * (1.0 - expval2theta);
//                else // To avoid linear Me-H-Me angles (6/6/06)
//                  expval12theta = p_val1 * -expval2theta;
//
//                CEval1 = Cf7ij * f7_jk * f8_Dj * expval12theta;
//                CEval2 = Cf7jk * f7_ij * f8_Dj * expval12theta;
//                CEval3 = Cf8j  * f7_ij * f7_jk * expval12theta;
//                CEval4 = -2.0 * p_val1 * p_val2 * f7_ij * f7_jk * f8_Dj *
//                  expval2theta * (theta_0 - theta);
//
//#ifdef SW_MATH
//                Ctheta_0 = p_val10 * DEG2RAD(theta_00) * exp_4_sw(-p_val10 * (2.0-SBO2));
//#else
//                Ctheta_0 = p_val10 * DEG2RAD(theta_00) * exp(-p_val10 * (2.0-SBO2));
//#endif
//
//                CEval5 = -CEval4 * Ctheta_0 * CSBO2;
//                CEval6 = CEval5 * dSBO1;
//                CEval7 = CEval5 * dSBO2;
//                CEval8 = -CEval4 / sin_theta;
//
//                data->my_en.e_ang += e_ang = f7_ij * f7_jk * f8_Dj * expval12theta;
//                /* END ANGLE ENERGY*/
//
//                /* PENALTY ENERGY */
//                p_pen1 = thbp->p_pen1;
//                //p_pen2 = system->reax_param.gp.l[19];
//                //p_pen3 = system->reax_param.gp.l[20];
//                //p_pen4 = system->reax_param.gp.l[21];
//
//#ifdef SW_MATH
//                exp_pen2ij = exp_4_sw( -p_pen2 * SQR( BOA_ij - 2.0 ) );
//                exp_pen2jk = exp_4_sw( -p_pen2 * SQR( BOA_jk - 2.0 ) );
//                exp_pen3   = exp_4_sw( -p_pen3 * workspace->Delta[j] );
//                exp_pen4   = exp_4_sw(  p_pen4 * workspace->Delta[j] );
//                
//#else
//                exp_pen2ij = exp( -p_pen2 * SQR( BOA_ij - 2.0 ) );
//                exp_pen2jk = exp( -p_pen2 * SQR( BOA_jk - 2.0 ) );
//                exp_pen3 = exp( -p_pen3 * workspace->Delta[j] );
//                exp_pen4 = exp(  p_pen4 * workspace->Delta[j] );
//#endif
//
//
//
//                trm_pen34 = 1.0 + exp_pen3 + exp_pen4;
//                f9_Dj = ( 2.0 + exp_pen3 ) / trm_pen34;
//                Cf9j = ( -p_pen3 * exp_pen3 * trm_pen34 -
//                         (2.0 + exp_pen3) * ( -p_pen3 * exp_pen3 +
//                        p_pen4 * exp_pen4 ) ) / SQR( trm_pen34 );
//
//                data->my_en.e_pen += e_pen = p_pen1 *f9_Dj *exp_pen2ij *exp_pen2jk;
//
//                CEpen1 = e_pen * Cf9j / f9_Dj;
//                temp   = -2.0 * p_pen2 * e_pen;
//                CEpen2 = temp * (BOA_ij - 2.0);
//                CEpen3 = temp * (BOA_jk - 2.0);
//                /* END PENALTY ENERGY */
//
//                /* COALITION ENERGY */
//                p_coa1 = thbp->p_coa1;
//                //p_coa2 = system->reax_param.gp.l[2];
//                //p_coa3 = system->reax_param.gp.l[38];
//                //p_coa4 = system->reax_param.gp.l[30];
//
//#ifdef SW_MATH
//                exp_coa2 = exp_4_sw( p_coa2 * workspace->Delta_val[j] );
//                data->my_en.e_coa += e_coa = p_coa1 / (1. + exp_coa2) *
//                  exp_4_sw( -p_coa3 * SQR(workspace->total_bond_order[i]-BOA_ij) ) *
//                  exp_4_sw( -p_coa3 * SQR(workspace->total_bond_order[k]-BOA_jk) ) *
//                  exp_4_sw( -p_coa4 * SQR(BOA_ij - 1.5) ) *
//                  exp_4_sw( -p_coa4 * SQR(BOA_jk - 1.5) );
//
//#else
//                exp_coa2 = exp( p_coa2 * workspace->Delta_val[j] );
//                data->my_en.e_coa += e_coa = p_coa1 / (1. + exp_coa2) *
//                  exp( -p_coa3 * SQR(workspace->total_bond_order[i]-BOA_ij) ) *
//                  exp( -p_coa3 * SQR(workspace->total_bond_order[k]-BOA_jk) ) *
//                  exp( -p_coa4 * SQR(BOA_ij - 1.5) ) *
//                  exp( -p_coa4 * SQR(BOA_jk - 1.5) );
//
//#endif
//
//                CEcoa1 = -2 * p_coa4 * (BOA_ij - 1.5) * e_coa;
//                CEcoa2 = -2 * p_coa4 * (BOA_jk - 1.5) * e_coa;
//                CEcoa3 = -p_coa2 * exp_coa2 * e_coa / (1 + exp_coa2);
//                CEcoa4 = -2 * p_coa3 * (workspace->total_bond_order[i]-BOA_ij) *e_coa;
//                CEcoa5 = -2 * p_coa3 * (workspace->total_bond_order[k]-BOA_jk) * e_coa;
//
//                /* END COALITION ENERGY */
//
//                /* FORCES */
//                bo_ij->Cdbo += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));
//                bo_jk->Cdbo += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));
//                workspace->CdDelta[j] += ((CEval3 + CEval7) + CEpen1 + CEcoa3);
//                workspace->CdDelta[i] += CEcoa4;
//                workspace->CdDelta[k] += CEcoa5;
//
//                for( t = start_j; t < end_j; ++t ) 
//                {
//                  pbond_jt = &( bonds->select.bond_list[t] );
//                  bo_jt = &(pbond_jt->bo_data);
//                  temp_bo_jt = bo_jt->BO;
//                  temp = CUBE( temp_bo_jt );
//                  pBOjt7 = temp * temp * temp_bo_jt;
//
//                  bo_jt->Cdbo += (CEval6 * pBOjt7);
//                  bo_jt->Cdbopi += CEval5;
//                  bo_jt->Cdbopi2 += CEval5;
//                }//for-t
//
//                if( control->virial == 0 )//virial == 0 
//                {
//                  rvec_ScaledAdd( workspace->f[i], CEval8, p_ijk->dcos_di );
//                  rvec_ScaledAdd( workspace->f[j], CEval8, p_ijk->dcos_dj );
//                  rvec_ScaledAdd( workspace->f[k], CEval8, p_ijk->dcos_dk );
//                }
//                else 
//                {
//                  rvec_Scale( force, CEval8, p_ijk->dcos_di );
//                  rvec_Add( workspace->f[i], force );
//                  rvec_iMultiply( ext_press, pbond_ij->rel_box, force );
//                  rvec_Add( data->my_ext_press, ext_press );
//
//                  rvec_ScaledAdd( workspace->f[j], CEval8, p_ijk->dcos_dj );
//
//                  rvec_Scale( force, CEval8, p_ijk->dcos_dk );
//                  rvec_Add( workspace->f[k], force );
//                  rvec_iMultiply( ext_press, pbond_jk->rel_box, force );
//                  rvec_Add( data->my_ext_press, ext_press );
//                }
//
//                /* tally into per-atom virials */
//                if( system->vflag_atom || system->evflag) 
//                {
//
//                  /* Acquire vectors */
//                  rvec_ScaledSum( delij, 1., system->my_atoms[i].x,
//                                  -1., system->my_atoms[j].x );
//                  rvec_ScaledSum( delkj, 1., system->my_atoms[k].x,
//                                  -1., system->my_atoms[j].x );
//
//                  rvec_Scale( fi_tmp, -CEval8, p_ijk->dcos_di );
//                  rvec_Scale( fj_tmp, -CEval8, p_ijk->dcos_dj );
//                  rvec_Scale( fk_tmp, -CEval8, p_ijk->dcos_dk );
//
//                  eng_tmp = e_ang + e_pen + e_coa;
//                  if (system->evflag)
//                  {
//                    //e_tally1_sys(system, i, eng_tmp);
//                    if (system->eflag_global){
//                      system->eng_vdwl += eng_tmp;
//                    }
//                    if (system->eflag_atom){
//                      system->eatom[i] += eng_tmp;
//                    }
//
//                    if (system->vflag_atom)
//                    {
//                      //v_tally3_sys(system, i, j, k, fi_tmp, fk_tmp, delij, delkj);
//
//                      double v[6];
//
//                      v[0] = THIRD * (delij[0]*fi_tmp[0] + delkj[0]*fj_tmp[0]);
//                      v[1] = THIRD * (delij[1]*fi_tmp[1] + delkj[1]*fj_tmp[1]);
//                      v[2] = THIRD * (delij[2]*fi_tmp[2] + delkj[2]*fj_tmp[2]);
//                      v[3] = THIRD * (delij[0]*fi_tmp[1] + delkj[0]*fj_tmp[1]);
//                      v[4] = THIRD * (delij[0]*fi_tmp[2] + delkj[0]*fj_tmp[2]);
//                      v[5] = THIRD * (delij[1]*fi_tmp[2] + delkj[1]*fj_tmp[2]);
//
//                      system->vatom[i][0] += v[0]; system->vatom[i][1] += v[1]; system->vatom[i][2] += v[2];
//                      system->vatom[i][3] += v[3]; system->vatom[i][4] += v[4]; system->vatom[i][5] += v[5];
//                      system->vatom[j][0] += v[0]; system->vatom[j][1] += v[1]; system->vatom[j][2] += v[2];
//                      system->vatom[j][3] += v[3]; system->vatom[j][4] += v[4]; system->vatom[j][5] += v[5];
//                      system->vatom[k][0] += v[0]; system->vatom[k][1] += v[1]; system->vatom[k][2] += v[2];
//                      system->vatom[k][3] += v[3]; system->vatom[k][4] += v[4]; system->vatom[k][5] += v[5];
//
//
//                    }
//                  }
//                }//if
//              }//if
//            }//for-cnt
//
//          }//if
//
//          //s12 = MPI_Wtime();
//        //if(myrank == 0)
//        //{
//        //  printf("%.10f, %.10f, %.10f\n", s12-s11, s12-s13, (s12-s13) / (s12 - s11));
//        //}
//
//        }//for-pk
//        //s9 = MPI_Wtime();
//        //if(myrank == 0)
//        //{
//        //  printf("%.10f, %.10f\n", s8-s7, s9-s8);
//        //}
//
//      }//if
//      thb_intrs->end_index[pi] = num_thb_intrs;
//      //pi_b = thb_intrs->end_index[pi];
//
//      //pi_b = num_thb_intrs;
//      //if(myrank == 0)
//      //  //if(pi_b - pi_a > 10)
//      //printf("j = %d, pi_b - pi_a = %d\n", j, pi_b - pi_a);
//    }//for-pi
//    //s6 = MPI_Wtime();
//
//
//    //s2 = MPI_Wtime();
//    //if(myrank == 0)
//    //{
//    //  //printf("N = %d\n", system->N);
//    //  printf("total:%.10f, for-t: %.10f , for-pi: %.10f, %.10f, %.10f\n", s2-s1, s4-s3, s6-s5, s8-s7, s9-s8);
//    //}
//
//  }//for-j
//  
//  //if(myrank == 0) 
//  //{
//  //  printf("num_thb_intrs: %d, %d\n", num_thb_intrs, thb_intrs->num_intrs);
//  //  printf("interaction_data:%d, total = %d\n", sizeof(three_body_interaction_data),  sizeof(three_body_interaction_data) * thb_intrs->num_intrs);
//
//  //}
//  //if(myrank == 0)
//  //  printf("total = %d\n\n\n", total);
//
//  if( num_thb_intrs >= thb_intrs->num_intrs * DANGER_ZONE ) 
//  {
//    workspace->realloc.num_3body = num_thb_intrs;/////
//    if( num_thb_intrs > thb_intrs->num_intrs ) 
//    {
//      fprintf( stderr, "step%d-ran out of space on angle_list: top=%d, max=%d",
//               data->step, num_thb_intrs, thb_intrs->num_intrs );
//      MPI_Abort( MPI_COMM_WORLD, INSUFFICIENT_MEMORY );
//    }
//  }
//
//}
//
//void Valence_Angles_C_part_compute_nopack(param_pack_t *param)
//{
//  //return ;
//  //if(athread_idle() == 0)
//  //  athread_init();
//
//  //athread_spawn(Valence_Angle_C_para, pm);
//  //athread_join();
//  int myrank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//  if(myrank == 0)
//    printf("part_compute\n");
//  reax_system_c *system = param->system;
//  control_params *control = param->control;
//  simulation_data *data = param->data;
//  storage *workspace = param->workspace;
//  reax_list **lists = param->lists;
//
//  int total_angles = 0;
//
//  int i, j, pi, k, pk, t;
//  int type_i, type_j, type_k;
//  int start_j, end_j;
//  int cnt;
//
//  double temp, temp_bo_jt, pBOjt7;
//  double p_val1, p_val2, p_val3, p_val4, p_val5;
//  double p_val6, p_val7, p_val8, p_val9, p_val10;
//  double p_pen1, p_pen2, p_pen3, p_pen4;
//  double p_coa1, p_coa2, p_coa3, p_coa4;
//  double trm8, expval6, expval7, expval2theta, expval12theta, exp3ij, exp3jk;
//  double exp_pen2ij, exp_pen2jk, exp_pen3, exp_pen4, trm_pen34, exp_coa2;
//  double dSBO1, dSBO2, SBO, SBO2, CSBO2, SBOp, prod_SBO, vlpadj;
//  double CEval1, CEval2, CEval3, CEval4, CEval5, CEval6, CEval7, CEval8;
//  double CEpen1, CEpen2, CEpen3;
//  double e_ang, e_coa, e_pen;
//  double CEcoa1, CEcoa2, CEcoa3, CEcoa4, CEcoa5;
//  double Cf7ij, Cf7jk, Cf8j, Cf9j;
//  double f7_ij, f7_jk, f8_Dj, f9_Dj;
//  double Ctheta_0, theta_0, theta_00,  sin_theta;
//  double BOA_ij, BOA_jk;
//  rvec force, ext_press;
//
//  // Tallying variables
//  double eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];
//  double delij[3], delkj[3];
//
//  three_body_header *thbh;
//  three_body_parameters *thbp;
//  three_body_interaction_data *p_ijk, *p_kji;
//  bond_data *pbond_ij, *pbond_jk, *pbond_jt;
//  bond_order_data *bo_ij, *bo_jk, *bo_jt;
//  reax_list *bonds = (*lists) + BONDS;
//  reax_list *thb_intrs =  (*lists) + THREE_BODIES;
//
//  /* global parameters used in these calculations */
//  p_val6 = system->reax_param.gp.l[14];
//  p_val8 = system->reax_param.gp.l[33];
//  p_val9 = system->reax_param.gp.l[16];
//  p_val10 = system->reax_param.gp.l[17];
//  
//  p_pen2 = system->reax_param.gp.l[19];
//  p_pen3 = system->reax_param.gp.l[20];
//  p_pen4 = system->reax_param.gp.l[21];
//
//  p_coa2 = system->reax_param.gp.l[2];
//  p_coa3 = system->reax_param.gp.l[38];
//  p_coa4 = system->reax_param.gp.l[30];
//  
//  for( j = 0; j < system->n; ++j ) 
//  { // Ray: the first one with system->N
//    type_j = system->my_atoms[j].type;
//    if (type_j < 0) continue;
//    start_j = bonds->index[j];
//    end_j   = bonds->end_index[j];
//
//    p_val3 = system->reax_param.sbp[type_j].p_val3;
//    p_val5 = system->reax_param.sbp[type_j].p_val5;
//
//    SBOp = 0, prod_SBO = 1;
//    for( t = start_j; t < end_j; ++t ) //start_j~end_j < 15
//    {
//      bo_jt = &(bonds->select.bond_list[t].bo_data);
//      SBOp += (bo_jt->BO_pi + bo_jt->BO_pi2);
//      temp = SQR( bo_jt->BO );
//      temp *= temp;
//      temp *= temp;
//#define SW_MATH
//#ifdef SW_MATH
//      prod_SBO *= exp_4_sw( -temp );
//#else
//      prod_SBO *= exp( -temp );
//#endif
//    }//for
//
//    if( workspace->vlpex[j] >= 0 )
//    {
//      vlpadj = 0;
//      dSBO2 = prod_SBO - 1;
//    }
//    else
//    {
//      vlpadj = workspace->nlp[j];
//      dSBO2 = (prod_SBO - 1) * (1 - p_val8 * workspace->dDelta_lp[j]);
//    }
//
//    SBO = SBOp + (1 - prod_SBO) * (-workspace->Delta_boc[j] - p_val8 * vlpadj);
//    dSBO1 = -8 * prod_SBO * ( workspace->Delta_boc[j] + p_val8 * vlpadj );
//
//    if( SBO <= 0 )
//      SBO2 = 0, CSBO2 = 0;
//    else if( SBO > 0 && SBO <= 1 ) 
//    {
//      SBO2 = pow( SBO, p_val9 );
//      CSBO2 = p_val9 * pow( SBO, p_val9 - 1 );
//    }
//    else if( SBO > 1 && SBO < 2 ) 
//    {
//      SBO2 = 2 - pow( 2-SBO, p_val9 );
//      CSBO2 = p_val9 * pow( 2 - SBO, p_val9 - 1 );
//    }
//    else
//      SBO2 = 2, CSBO2 = 0;
//
//#ifdef SW_MATH
//    expval6 = exp_4_sw( p_val6 * workspace->Delta_boc[j] );
//#else
//    expval6 = exp( p_val6 * workspace->Delta_boc[j] );
//#endif
//    //if(myrank == 0 && end_j - start_j >= 15)
//    //  printf("end_j - start_j = %d\n", end_j - start_j);
//    for( pi = start_j; pi < end_j; ++pi ) 
//    {
//      pbond_ij = &(bonds->select.bond_list[pi]);
//      bo_ij = &(pbond_ij->bo_data);
//      BOA_ij = bo_ij->BO - control->thb_cut;
//      
//      if(BOA_ij <= 0.0) continue;
//
//      
//      i = pbond_ij->nbr;
//      type_i = system->my_atoms[i].type;///
//
//      int num_thb_intrs = thb_intrs->middle_index[pi];
//      //if(myrank == 0)
//      //  printf("%d, %d\n", end_j - start_j, end_j-pi);
//
//      for( pk = pi+1; pk < end_j; ++pk ) 
//      {
//        pbond_jk = &(bonds->select.bond_list[pk]);
//        bo_jk    = &(pbond_jk->bo_data);
//        BOA_jk   = bo_jk->BO - control->thb_cut;
//        k        = pbond_jk->nbr;
//        type_k   = system->my_atoms[k].type;///
//        p_ijk    = &( thb_intrs->select.three_body_list[num_thb_intrs] );
//
//        //if(myrank == 0) 
//        //  printf("j=%d, i = %d, k = %d\n", j, i, k);
//
//        sin_theta = sin(p_ijk->theta);
//        if(sin_theta < 1.0e-5) sin_theta = 1.0e-5;
//
//        ++num_thb_intrs;
//
//        if((BOA_jk <= 0.0) || (bo_ij->BO*bo_jk->BO <= control->thb_cutsq)) continue;
//        
//        thbh = &( system->reax_param.thbp[type_i][type_j][type_k]);///
//        for( cnt = 0; cnt < thbh->cnt; ++cnt ) //thbh->cnt < 3
//        {
//          if( fabs(thbh->prm[cnt].p_val1) <= 0.001 ) continue;
//          total_angles++;
//          
//          thbp = &(thbh->prm[cnt]);
//
//          /* ANGLE ENERGY */
//          p_val1    = thbp->p_val1;
//          p_val2    = thbp->p_val2;
//          p_val4    = thbp->p_val4;
//          p_val7    = thbp->p_val7;
//          theta_00  = thbp->theta_00;
//
//#ifdef SW_MATH
//          exp3ij = exp_4_sw( -p_val3 * pow( BOA_ij, p_val4 ) );
//#else
//          exp3ij = exp( -p_val3 * pow( BOA_ij, p_val4 ) );
//#endif
//
//          f7_ij = 1.0 - exp3ij;
//          Cf7ij = p_val3 * p_val4 * pow( BOA_ij, p_val4 - 1.0 ) * exp3ij;
//
//#ifdef SW_MATH
//          exp3jk = exp_4_sw( -p_val3 * pow( BOA_jk, p_val4 ) );
//#else
//          exp3jk = exp( -p_val3 * pow( BOA_jk, p_val4 ) );
//#endif
//
//          f7_jk = 1.0 - exp3jk;
//          Cf7jk = p_val3 * p_val4 * pow( BOA_jk, p_val4 - 1.0 ) * exp3jk;
//
//#ifdef SW_MATH
//          expval7 = exp_4_sw( -p_val7 * workspace->Delta_boc[j] );
//#else
//          expval7 = exp( -p_val7 * workspace->Delta_boc[j] );
//#endif
//
//          trm8 = 1.0 + expval6 + expval7;
//          f8_Dj = p_val5 - ( (p_val5 - 1.0) * (2.0 + expval6) / trm8 );
//          Cf8j = ( (1.0 - p_val5) / SQR(trm8) ) *( p_val6 * expval6 * trm8 -
//              (2.0 + expval6) * ( p_val6*expval6 - p_val7*expval7 ) );
//
//#ifdef SW_MATH
//          theta_0 = 180.0 - theta_00 * (1.0 - exp_4_sw(-p_val10 * (2.0-SBO2)));
//#else
//          theta_0 = 180.0 - theta_00 * (1.0 - exp(-p_val10 * (2.0 - SBO2)));
//#endif
//
//          theta_0 = DEG2RAD( theta_0 );
//
//#ifdef SW_MATH
//          expval2theta  = exp_4_sw( -p_val2 * SQR(theta_0 - p_ijk->theta) );
//#else
//          expval2theta  = exp( -p_val2 * SQR(theta_0 - p_ijk->theta) );
//#endif
//
//          if( p_val1 >= 0 )
//            expval12theta = p_val1 * (1.0 - expval2theta);
//          else // To avoid linear Me-H-Me angles (6/6/06)
//            expval12theta = p_val1 * -expval2theta;
//
//          CEval1 = Cf7ij * f7_jk * f8_Dj * expval12theta;
//          CEval2 = Cf7jk * f7_ij * f8_Dj * expval12theta;
//          CEval3 = Cf8j  * f7_ij * f7_jk * expval12theta;
//          CEval4 = -2.0 * p_val1 * p_val2 * f7_ij * f7_jk * f8_Dj *
//            expval2theta * (theta_0 - p_ijk->theta);
//
//#ifdef SW_MATH
//          Ctheta_0 = p_val10 * DEG2RAD(theta_00) * exp_4_sw(-p_val10 * (2.0-SBO2));
//#else
//          Ctheta_0 = p_val10 * DEG2RAD(theta_00) * exp(-p_val10 * (2.0-SBO2));
//#endif
//
//          CEval5 = -CEval4 * Ctheta_0 * CSBO2;
//          CEval6 = CEval5 * dSBO1;
//          CEval7 = CEval5 * dSBO2;
//          CEval8 = -CEval4 / sin_theta;
//
//          data->my_en.e_ang += e_ang = f7_ij * f7_jk * f8_Dj * expval12theta;
//          /* END ANGLE ENERGY*/
//
//          /* PENALTY ENERGY */
//          p_pen1 = thbp->p_pen1;
//
//#ifdef SW_MATH
//          exp_pen2ij = exp_4_sw( -p_pen2 * SQR( BOA_ij - 2.0 ) );
//          exp_pen2jk = exp_4_sw( -p_pen2 * SQR( BOA_jk - 2.0 ) );
//          exp_pen3   = exp_4_sw( -p_pen3 * workspace->Delta[j] );
//          exp_pen4   = exp_4_sw(  p_pen4 * workspace->Delta[j] );
//          
//#else
//          exp_pen2ij = exp( -p_pen2 * SQR( BOA_ij - 2.0 ) );
//          exp_pen2jk = exp( -p_pen2 * SQR( BOA_jk - 2.0 ) );
//          exp_pen3 = exp( -p_pen3 * workspace->Delta[j] );
//          exp_pen4 = exp(  p_pen4 * workspace->Delta[j] );
//#endif
//
//          trm_pen34 = 1.0 + exp_pen3 + exp_pen4;
//          f9_Dj = ( 2.0 + exp_pen3 ) / trm_pen34;
//          Cf9j = ( -p_pen3 * exp_pen3 * trm_pen34 -
//                   (2.0 + exp_pen3) * ( -p_pen3 * exp_pen3 +
//                  p_pen4 * exp_pen4 ) ) / SQR( trm_pen34 );
//
//          data->my_en.e_pen += e_pen = p_pen1 *f9_Dj *exp_pen2ij *exp_pen2jk;
//
//          CEpen1 = e_pen * Cf9j / f9_Dj;
//          temp   = -2.0 * p_pen2 * e_pen;
//          CEpen2 = temp * (BOA_ij - 2.0);
//          CEpen3 = temp * (BOA_jk - 2.0);
//          /* END PENALTY ENERGY */
//
//          /* COALITION ENERGY */
//          p_coa1 = thbp->p_coa1;
//
//#ifdef SW_MATH
//          exp_coa2 = exp_4_sw( p_coa2 * workspace->Delta_val[j] );
//          data->my_en.e_coa += e_coa = p_coa1 / (1. + exp_coa2) *
//            exp_4_sw( -p_coa3 * SQR(workspace->total_bond_order[i]-BOA_ij) ) *
//            exp_4_sw( -p_coa3 * SQR(workspace->total_bond_order[k]-BOA_jk) ) *
//            exp_4_sw( -p_coa4 * SQR(BOA_ij - 1.5) ) *
//            exp_4_sw( -p_coa4 * SQR(BOA_jk - 1.5) );
//#else
//          exp_coa2 = exp( p_coa2 * workspace->Delta_val[j] );
//          data->my_en.e_coa += e_coa = p_coa1 / (1. + exp_coa2) *
//            exp( -p_coa3 * SQR(workspace->total_bond_order[i]-BOA_ij) ) *
//            exp( -p_coa3 * SQR(workspace->total_bond_order[k]-BOA_jk) ) *
//            exp( -p_coa4 * SQR(BOA_ij - 1.5) ) *
//            exp( -p_coa4 * SQR(BOA_jk - 1.5) );
//
//#endif
//
//          CEcoa1 = -2 * p_coa4 * (BOA_ij - 1.5) * e_coa;
//          CEcoa2 = -2 * p_coa4 * (BOA_jk - 1.5) * e_coa;
//          CEcoa3 = -p_coa2 * exp_coa2 * e_coa / (1 + exp_coa2);
//          CEcoa4 = -2 * p_coa3 * (workspace->total_bond_order[i]-BOA_ij) *e_coa;
//          CEcoa5 = -2 * p_coa3 * (workspace->total_bond_order[k]-BOA_jk) * e_coa;
//          /* END COALITION ENERGY */
//
//          /* FORCES */
//          bo_ij->Cdbo += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));
//          bo_jk->Cdbo += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));
//          workspace->CdDelta[j] += ((CEval3 + CEval7) + CEpen1 + CEcoa3);
//          workspace->CdDelta[i] += CEcoa4;
//          workspace->CdDelta[k] += CEcoa5;
//          //if(myrank == 0)
//          //  printf("%lf, %lf, %lf\n", workspace->CdDelta[j], workspace->CdDelta[k], workspace->CdDelta[k]);
//
//          for( t = start_j; t < end_j; ++t ) 
//          {
//            pbond_jt = &( bonds->select.bond_list[t] );
//            bo_jt = &(pbond_jt->bo_data);
//            temp_bo_jt = bo_jt->BO;
//            temp = CUBE( temp_bo_jt );
//            pBOjt7 = temp * temp * temp_bo_jt;
//
//            bo_jt->Cdbo += (CEval6 * pBOjt7);
//            bo_jt->Cdbopi += CEval5;
//            bo_jt->Cdbopi2 += CEval5;
//          }//for-t
//
//          //////////// control->virial == 0/////////////
//          workspace->f[i][0] += CEval8 * (p_ijk->dcos_di)[0];
//          workspace->f[i][1] += CEval8 * (p_ijk->dcos_di)[1];
//          workspace->f[i][2] += CEval8 * (p_ijk->dcos_di)[2];
//
//          workspace->f[j][0] += CEval8 * (p_ijk->dcos_dj)[0];
//          workspace->f[j][1] += CEval8 * (p_ijk->dcos_dj)[1];
//          workspace->f[j][2] += CEval8 * (p_ijk->dcos_dj)[2];
//
//          workspace->f[k][0] += CEval8 * (p_ijk->dcos_dk)[0];
//          workspace->f[k][1] += CEval8 * (p_ijk->dcos_dk)[1];
//          workspace->f[k][2] += CEval8 * (p_ijk->dcos_dk)[2];
//          
//          /* tally into per-atom virials */
//          if( system->vflag_atom || system->evflag) 
//          {
//            delij[0] = system->my_atoms[i].x[0] - system->my_atoms[j].x[0];
//            delij[1] = system->my_atoms[i].x[1] - system->my_atoms[j].x[1];
//            delij[2] = system->my_atoms[i].x[2] - system->my_atoms[j].x[2];
//            
//            delkj[0] = system->my_atoms[k].x[0] - system->my_atoms[j].x[0];
//            delkj[1] = system->my_atoms[k].x[1] - system->my_atoms[j].x[1];
//            delkj[2] = system->my_atoms[k].x[2] - system->my_atoms[j].x[2];
//
//            fi_tmp[0] = -CEval8, p_ijk->dcos_di[0];
//            fi_tmp[1] = -CEval8, p_ijk->dcos_di[1];
//            fi_tmp[2] = -CEval8, p_ijk->dcos_di[2];
//
//            fj_tmp[0] = -CEval8, p_ijk->dcos_dj[0];
//            fj_tmp[1] = -CEval8, p_ijk->dcos_dj[1];
//            fj_tmp[2] = -CEval8, p_ijk->dcos_dj[2];
//
//            fk_tmp[0] = -CEval8, p_ijk->dcos_dk[0];
//            fk_tmp[1] = -CEval8, p_ijk->dcos_dk[1];
//            fk_tmp[2] = -CEval8, p_ijk->dcos_dk[2];
//
//            eng_tmp = e_ang + e_pen + e_coa;
//            if (system->evflag)
//            {
//              if (system->eflag_global){
//                system->eng_vdwl += eng_tmp;
//              }
//              if (system->eflag_atom){
//                system->eatom[i] += eng_tmp;
//              }
//
//              if (system->vflag_atom)
//              {
//                double v[6];
//                v[0] = THIRD * (delij[0]*fi_tmp[0] + delkj[0]*fj_tmp[0]);
//                v[1] = THIRD * (delij[1]*fi_tmp[1] + delkj[1]*fj_tmp[1]);
//                v[2] = THIRD * (delij[2]*fi_tmp[2] + delkj[2]*fj_tmp[2]);
//                v[3] = THIRD * (delij[0]*fi_tmp[1] + delkj[0]*fj_tmp[1]);
//                v[4] = THIRD * (delij[0]*fi_tmp[2] + delkj[0]*fj_tmp[2]);
//                v[5] = THIRD * (delij[1]*fi_tmp[2] + delkj[1]*fj_tmp[2]);
//
//                system->vatom[i][0] += v[0]; 
//                system->vatom[i][1] += v[1]; 
//                system->vatom[i][2] += v[2];
//                system->vatom[i][3] += v[3]; 
//                system->vatom[i][4] += v[4]; 
//                system->vatom[i][5] += v[5];
//
//                system->vatom[j][0] += v[0]; 
//                system->vatom[j][1] += v[1]; 
//                system->vatom[j][2] += v[2];
//                system->vatom[j][3] += v[3]; 
//                system->vatom[j][4] += v[4]; 
//                system->vatom[j][5] += v[5];
//
//                system->vatom[k][0] += v[0]; 
//                system->vatom[k][1] += v[1]; 
//                system->vatom[k][2] += v[2];
//                system->vatom[k][3] += v[3]; 
//                system->vatom[k][4] += v[4]; 
//                system->vatom[k][5] += v[5];
//
//              }
//            }
//          }//if
//        }//for-cnt
//      }//for-pk
//      //if(myrank == 0)
//      //  printf("\n\n");
//    }//for-pi
//    
//    //if(myrank == 0)
//    //    printf("\n\n");
//
//  }//for-j
//  //if(myrank == 0)
//  //  printf("total_angles = %d\n", total_angles);
//}
//
////////////////////////////////////////////////////////////////////////////////
//
//
//void Valence_Angles_C_part_init(param_pack_t *param)
//{
//  //return;
//  
//  reax_system_c *system = param->system;
//  control_params *control = param->control;
//  simulation_data *data = param->data;
//  storage *workspace = param->workspace;
//  reax_list **lists = param->lists;
//  
//  int myrank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//
//  if(myrank == 0)  
//  {
//        
//    
//  //  PP[0]   = 3.1415926535897932384626433832795/2 ;
//  //  PP[1]   = 3.0/6 ;
//  //  PP[2]   = 3.0/40 ;
//  //  PP[3]   = 5.0/112 ;
//  //  PP[4]   = 35.0/1152 ;
//  //  PP[5]   = 63.0/2816 ;
//  //  PP[6]   = 231.0/13312 ;
//  //  PP[7]   = 143.0/10240 ;
//  //  PP[8]   = 6435.0/557056 ;
//  //  PP[9]   = 12155.0/1245184 ;
//  //  PP[10]  = 46189.0/5505024 ;
//  //  PP[11]  = 88179.0/12058624 ;
//  //  PP[12]  = 676039.0/104857600 ;
//  //  PP[13]  = 1300075.0/226492416;
//  //  PP[14]  = 5014575.0/973078528;
//  //  PP[15]  = 9694845.0/2080374784;
//  //  PP[16]  = 100180065.0/23622320128;
//  //  PP[17]  = 116680311.0/30064771072;
//  //  PP[18]  = 2268783825.0/635655159808;
//  //  PP[19]  = 1472719325.0/446676598784;
//  //  PP[20]  = 34461632205.0/11269994184704;
//  //  PP[21]  = 67282234305.0/23639499997184;
//  //  PP[22]  = 17534158031.0/6597069766656;
//  //  PP[23]  = 514589420475.0/206708186021888;
//  //  PP[24]  = 8061900920775.0/3448068464705536;
//  //  PP[25]  = 5267108601573.0/2392537302040576; 
//  //  PP[26] = 61989816618513.0/29836347531329536;
//  //  PP[27] = 121683714103007.0/61924494876344320;
//  //  PP[28] = 956086325095055.0/513410357520236544;
//  //  PP[29] = 1879204156221315.0/1062849512059437056;
//  //  PP[30] = 7391536347803839.0/4395513236313604096;
//  //  PP[31] =2077805148460987.0/1297036692682702848 ;
//  //  //PP[32] =916312070471295267.0/599519182395560427520 ;
//  //  //PP[33] =1804857108504066435.0/1235931852938539958272 ;
//  //  //PP[34] =2371086789603381395.0/1697100454781278748672 ;
//  //  
//
//  //  printf("P0:   %.40f\n",PP[0]);
//  //  printf("P1:   %.40f\n",PP[1]);
//  //  printf("P2:   %.40f\n",PP[2]);
//  //  printf("P3:   %.40f\n",PP[3]);
//  //  printf("P4:   %.40f\n",PP[4]);
//  //  printf("P5:   %.40f\n",PP[5]);
//  //  printf("P6:   %.40f\n",PP[6]);
//  //  printf("P7:   %.40f\n",PP[7]);
//  //  printf("P8:   %.40f\n",PP[8]);
//  //  printf("P9:   %.40f\n",PP[9]);
//  //  printf("P10:  %.40f\n",PP[10]);
//  //  printf("P11:  %.40f\n",PP[11]);
//  //  printf("P12:  %.40f\n",PP[12]);
//  //  printf("P13:  %.40f\n",PP[13]);
//  //  printf("P14:  %.40f\n",PP[14]);
//  //  printf("P15:  %.40f\n",PP[15]);
//  //  printf("P16:  %.40f\n",PP[16]);
//  //  printf("P17:  %.40f\n",PP[17]);
//  //  printf("P18:  %.40f\n",PP[18]);
//  //  printf("P19:  %.40f\n",PP[19]);
//  //  printf("P20:  %.40f\n",PP[20]);
//  //  printf("P21:  %.40f\n",PP[21]);
//  //  printf("P22:  %.40f\n",PP[22]);
//  //  printf("P23:  %.40f\n",PP[23]);
//  //  printf("P24:  %.40f\n",PP[24]);
//  //  printf("P25:  %.40f\n",PP[25]);
//  //  printf("P26:  %.40f\n",PP[26]);
//  //  printf("P27:  %.40f\n",PP[27]);
//  //  printf("P28:  %.40f\n",PP[28]);
//  //  printf("P29:  %.40f\n",PP[29]);
//  //  printf("P30:  %.40f\n",PP[30]);
//  //  printf("P31:  %.40f\n",PP[31]);
//  //  //printf("P32:  %.40f\n",PP[32]);
//  //  //printf("P33:  %.40f\n",PP[33]);
//  //  //printf("P34:  %.40f\n",PP[34]);
//
//  //  //TT[1]   = PP[3] / PP[1];//p3/p1;
//  //  //TT[2]   = PP[5] / PP[3];//p3/p1;
//  //  //TT[3]   = PP[7] / PP[5];//p3/p1;
//  //  //TT[4]   = PP[9] / PP[7];//p3/p1;
//  //  //TT[5]   = PP[11] / PP[9];//p3/p1;
//  //  //TT[6]   = PP[13] / PP[11];//p3/p1;
//  //  //TT[7]   = PP[15] / PP[13];//p3/p1;
//  //  //TT[8]   = PP[17] / PP[15];//p3/p1;
//  //  //TT[9]   = PP[19] / PP[17];//p3/p1;
//  //  //TT[10]  = PP[21] / PP[19];//p3/p1;
//  //  //TT[11]  = PP[23] / PP[21];//p3/p1;
//  //  //TT[12]  = PP[25] / PP[23];//p3/p1;
//  //  //
//  //  //printf("TT0:   %.40f\n",TT[0]);
//  //  //printf("TT1:   %.40f\n",TT[1]);
//  //  //printf("TT2:   %.40f\n",TT[2]);
//  //  //printf("TT3:   %.40f\n",TT[3]);
//  //  //printf("TT4:   %.40f\n",TT[4]);
//  //  //printf("TT5:   %.40f\n",TT[5]);
//  //  //printf("TT6:   %.40f\n",TT[6]);
//  //  //printf("TT7:   %.40f\n",TT[7]);
//  //  //printf("TT8:   %.40f\n",TT[8]);
//  //  //printf("TT9:   %.40f\n",TT[9]);
//  //  //printf("TT10:  %.40f\n",TT[10]);
//  //  //printf("TT11:  %.40f\n",TT[11]);
//  //  //printf("TT12:  %.40f\n",TT[12]);
//
//  }
//
//
//  int total = 0;
//
//  int i, j, pi, k, pk, t, z;//, pz;
//  int type_i, type_j, type_k;
//  int start_j, end_j, start_pk, end_pk, start_pz, end_pz, start_pi, end_pi;
//  int cnt, num_thb_intrs;
//
//  double theta, cos_theta;
//  double BOA_ij, BOA_jk;//, BOA_jz;
//
//  three_body_header *thbh;
//  three_body_parameters *thbp;
//  three_body_interaction_data *p_ijk, *p_kji;
//  bond_data *pbond_ij, *pbond_jk, *pbond_jt, *pbond_jz;
//  bond_order_data *bo_ij, *bo_jk, *bo_jt, *bo_jz;
//  reax_list *bonds      = (*lists) + BONDS;
//  reax_list *thb_intrs =  (*lists) + THREE_BODIES;
//
//  num_thb_intrs = 0;
//  
//  //printf("in valence  init---\n");
//  //if(myrank == 0) 
//  //{
//  //  printf("n = %d, N = %d\n", system->n, system->N);
//  //}
//
//  for( j = 0; j < system->N; ++j ) 
//  { // Ray: the first one with system->N
//    type_j = system->my_atoms[j].type;
//    if (type_j < 0) continue;
//    start_j = bonds->index[j];
//    end_j   = bonds->end_index[j];
//
//    for( pk = start_j; pk < end_j; ++pk) //pk->pz
//    {
//      thb_intrs->index[pk] = num_thb_intrs;
//      pbond_jk  = &(bonds->select.bond_list[pk]);
//      bo_jk     = &(pbond_jk->bo_data);
//      BOA_jk    = bo_jk->BO - control->thb_cut;
//
//      if(BOA_jk > 0.0 && (j < system->n || pbond_jk->nbr < system->n)) 
//      {
//        k = pbond_jk->nbr;
//
//        //if(myrank == 0)
//        //  printf("before: num_thb_intrs=%d, %d\n", num_thb_intrs, pk-start_j+1);
//        int test_num = num_thb_intrs;
//        for( pi = start_j; pi < pk; ++pi) 
//        {
//          start_pi = thb_intrs->index[pi];
//          end_pi   = thb_intrs->end_index[pi];
//          
//          for(t = start_pi; t < end_pi; ++t)
//          {
//            if(thb_intrs->select.three_body_list[t].thb == k) 
//            {
//              p_ijk = &(thb_intrs->select.three_body_list[num_thb_intrs] );
//              p_kji = &(thb_intrs->select.three_body_list[t]);
//
//              p_ijk->thb    = bonds->select.bond_list[pi].nbr;
//              p_ijk->pthb   = pi;
//              p_ijk->theta  = p_kji->theta;
//              rvec_Copy(p_ijk->dcos_di, p_kji->dcos_dk);
//              rvec_Copy(p_ijk->dcos_dj, p_kji->dcos_dj);
//              rvec_Copy(p_ijk->dcos_dk, p_kji->dcos_di);
//
//              ++num_thb_intrs;
//              break;
//            }//if
//          }//for-t
//        }//for-pi
//        //if(myrank == 0)
//        //  if(test_num + pk - start_j + 1 != num_thb_intrs)
//        //    printf("after: num_thb_intrs=%d, %d, %d\n", num_thb_intrs, test_num, pk-start_j+1);
//
//
//        thb_intrs->middle_index[pk] = num_thb_intrs;//////
//
//        for(pi = pk+1; pi < end_j; ++pi) 
//        {
//          pbond_ij = &(bonds->select.bond_list[pi]);
//          i        = pbond_ij->nbr;
//          p_ijk    = &(thb_intrs->select.three_body_list[num_thb_intrs] );
//
//          Calculate_Theta( pbond_jk->dvec, pbond_jk->d,
//                           pbond_ij->dvec, pbond_ij->d,
//                           &theta, &cos_theta );
//
//          Calculate_dCos_Theta( pbond_jk->dvec, pbond_jk->d,
//                                pbond_ij->dvec, pbond_ij->d,
//                                &(p_ijk->dcos_di), &(p_ijk->dcos_dj),
//                                &(p_ijk->dcos_dk) );
//          p_ijk->thb = i;
//          p_ijk->pthb = pi;
//          p_ijk->theta = theta;
//          
//          ++num_thb_intrs;
//        }//for-pi
//
//      }//if
//      thb_intrs->end_index[pk] = num_thb_intrs;
//
//    }//for-pk
//
//  }//for-j
//    
//  //if(myrank == 0)
//  //  printf("total: num_thb_intrs = %d\n", num_thb_intrs);
//
//  if( num_thb_intrs >= thb_intrs->num_intrs * DANGER_ZONE ) 
//  {
//    workspace->realloc.num_3body = num_thb_intrs;/////
//    if( num_thb_intrs > thb_intrs->num_intrs ) 
//    {
//      fprintf( stderr, "step%d-ran out of space on angle_list: top=%d, max=%d",
//               data->step, num_thb_intrs, thb_intrs->num_intrs );
//      MPI_Abort( MPI_COMM_WORLD, INSUFFICIENT_MEMORY );
//    }
//  }
//
//}
//
//
//
//////////////////////////////////////////////////////
//void Valence_Angles_C_part_compute(param_pack_t *param)
//{
//  //return ;
//  //if(athread_idle() == 0)
//  //  athread_init();
//
//  //athread_spawn(Valence_Angle_C_para, pm);
//  //athread_join();
//  int myrank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//  //if(myrank == 0)
//  //  printf("part_compute\n");
//  reax_system_c *system   = param->system;
//  control_params *control = param->control;
//  simulation_data *data   = param->data;
//  storage *workspace      = param->workspace;
//  reax_list **lists       = param->lists;
//
//  int total_angles = 0;
//
//  int i, j, pi, k, pk, t;
//  int type_i, type_j, type_k;
//  int start_j, end_j;
//  int cnt;
//
//  double temp, temp_bo_jt, pBOjt7;
//  double p_val1, p_val2, p_val3, p_val4, p_val5;
//  double p_val6, p_val7, p_val8, p_val9, p_val10;
//  double p_pen1, p_pen2, p_pen3, p_pen4;
//  double p_coa1, p_coa2, p_coa3, p_coa4;
//  double trm8, expval6, expval7, expval2theta, expval12theta, exp3ij, exp3jk;
//  double exp_pen2ij, exp_pen2jk, exp_pen3, exp_pen4, trm_pen34, exp_coa2;
//  double dSBO1, dSBO2, SBO, SBO2, CSBO2, SBOp, prod_SBO, vlpadj;
//  double CEval1, CEval2, CEval3, CEval4, CEval5, CEval6, CEval7, CEval8;
//  double CEpen1, CEpen2, CEpen3;
//  double e_ang, e_coa, e_pen;
//  double CEcoa1, CEcoa2, CEcoa3, CEcoa4, CEcoa5;
//  double Cf7ij, Cf7jk, Cf8j, Cf9j;
//  double f7_ij, f7_jk, f8_Dj, f9_Dj;
//  double Ctheta_0, theta_0, theta_00,  sin_theta;
//  double BOA_ij, BOA_jk;
//  rvec force, ext_press;
//
//  // Tallying variables
//  double eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];
//  double delij[3], delkj[3];
//
//  three_body_header *thbh;
//  three_body_parameters *thbp;
//  three_body_interaction_data *p_ijk, *p_kji;
//  bond_data *pbond_ij, *pbond_jk, *pbond_jt;
//  bond_order_data *bo_ij, *bo_jk, *bo_jt;
//  reax_list *bonds = (*lists) + BONDS;
//  reax_list *thb_intrs =  (*lists) + THREE_BODIES;
//
//  /* global parameters used in these calculations */
//  p_val6 = system->reax_param.gp.l[14];
//  p_val8 = system->reax_param.gp.l[33];
//  p_val9 = system->reax_param.gp.l[16];
//  p_val10 = system->reax_param.gp.l[17];
//  
//  p_pen2 = system->reax_param.gp.l[19];
//  p_pen3 = system->reax_param.gp.l[20];
//  p_pen4 = system->reax_param.gp.l[21];
//
//  p_coa2 = system->reax_param.gp.l[2];
//  p_coa3 = system->reax_param.gp.l[38];
//  p_coa4 = system->reax_param.gp.l[30];
//  
//  for( j = 0; j < system->n; ++j ) 
//  { // Ray: the first one with system->N
//    type_j = system->my_atoms[j].type;
//    if (type_j < 0) continue;
//    start_j = bonds->index[j];
//    end_j   = bonds->end_index[j];
//
//    p_val3 = system->reax_param.sbp[type_j].p_val3;
//    p_val5 = system->reax_param.sbp[type_j].p_val5;
//
//    SBOp = 0, prod_SBO = 1;
//    for( t = start_j; t < end_j; ++t ) //start_j~end_j < 15
//    {
//      bo_jt = &(bonds->select.bond_list[t].bo_data);
//      SBOp += (bo_jt->BO_pi + bo_jt->BO_pi2);
//      temp = SQR( bo_jt->BO );
//      temp *= temp;
//      temp *= temp;
//#define SW_MATH
//#ifdef SW_MATH
//      prod_SBO *= exp_4_sw( -temp );
//#else
//      prod_SBO *= exp( -temp );
//#endif
//    }//for
//
//    if( workspace->vlpex[j] >= 0 )
//    {
//      vlpadj = 0;
//      dSBO2 = prod_SBO - 1;
//    }
//    else
//    {
//      vlpadj = workspace->nlp[j];
//      dSBO2 = (prod_SBO - 1) * (1 - p_val8 * workspace->dDelta_lp[j]);
//    }
//
//    SBO = SBOp + (1 - prod_SBO) * (-workspace->Delta_boc[j] - p_val8 * vlpadj);
//    dSBO1 = -8 * prod_SBO * ( workspace->Delta_boc[j] + p_val8 * vlpadj );
//
//    if( SBO <= 0 )
//      SBO2 = 0, CSBO2 = 0;
//    else if( SBO > 0 && SBO <= 1 ) 
//    {
//      SBO2 = pow( SBO, p_val9 );
//      CSBO2 = p_val9 * pow( SBO, p_val9 - 1 );
//    }
//    else if( SBO > 1 && SBO < 2 ) 
//    {
//      SBO2 = 2 - pow( 2-SBO, p_val9 );
//      CSBO2 = p_val9 * pow( 2 - SBO, p_val9 - 1 );
//    }
//    else
//      SBO2 = 2, CSBO2 = 0;
//
//#ifdef SW_MATH
//    expval6 = exp_4_sw( p_val6 * workspace->Delta_boc[j] );
//#else
//    expval6 = exp( p_val6 * workspace->Delta_boc[j] );
//#endif
//    
//    for( pi = start_j; pi < end_j; ++pi ) 
//    {
//      pbond_ij = &(bonds->select.bond_list[pi]);
//      bo_ij = &(pbond_ij->bo_data);
//      BOA_ij = bo_ij->BO - control->thb_cut;
//      
//      if(BOA_ij <= 0.0) continue;
//
//      
//      i = pbond_ij->nbr;
//      type_i = system->my_atoms[i].type;///
//
//      int num_thb_intrs = thb_intrs->middle_index[pi];
//
//      for( pk = pi+1; pk < end_j; ++pk ) 
//      {
//        pbond_jk = &(bonds->select.bond_list[pk]);
//        bo_jk    = &(pbond_jk->bo_data);
//        BOA_jk   = bo_jk->BO - control->thb_cut;
//        k        = pbond_jk->nbr;
//        type_k   = system->my_atoms[k].type;///
//        p_ijk    = &( thb_intrs->select.three_body_list[num_thb_intrs] );
//
//
//        sin_theta = sin(p_ijk->theta);
//        if(sin_theta < 1.0e-5) sin_theta = 1.0e-5;
//
//        ++num_thb_intrs;
//
//        if((BOA_jk <= 0.0) || (bo_ij->BO*bo_jk->BO <= control->thb_cutsq)) continue;
//        
//        thbh = &( system->reax_param.thbp[type_i][type_j][type_k]);///
//        for( cnt = 0; cnt < thbh->cnt; ++cnt ) //thbh->cnt < 3
//        {
//          if( fabs(thbh->prm[cnt].p_val1) <= 0.001 ) continue;
//          total_angles++;
//          
//          thbp = &(thbh->prm[cnt]);
//
//          /* ANGLE ENERGY */
//          p_val1    = thbp->p_val1;
//          p_val2    = thbp->p_val2;
//          p_val4    = thbp->p_val4;
//          p_val7    = thbp->p_val7;
//          theta_00  = thbp->theta_00;
//
//#ifdef SW_MATH
//          exp3ij = exp_4_sw( -p_val3 * pow( BOA_ij, p_val4 ) );
//#else
//          exp3ij = exp( -p_val3 * pow( BOA_ij, p_val4 ) );
//#endif
//
//          f7_ij = 1.0 - exp3ij;
//          Cf7ij = p_val3 * p_val4 * pow( BOA_ij, p_val4 - 1.0 ) * exp3ij;
//
//#ifdef SW_MATH
//          exp3jk = exp_4_sw( -p_val3 * pow( BOA_jk, p_val4 ) );
//#else
//          exp3jk = exp( -p_val3 * pow( BOA_jk, p_val4 ) );
//#endif
//
//          f7_jk = 1.0 - exp3jk;
//          Cf7jk = p_val3 * p_val4 * pow( BOA_jk, p_val4 - 1.0 ) * exp3jk;
//
//#ifdef SW_MATH
//          expval7 = exp_4_sw( -p_val7 * workspace->Delta_boc[j] );
//#else
//          expval7 = exp( -p_val7 * workspace->Delta_boc[j] );
//#endif
//
//          trm8 = 1.0 + expval6 + expval7;
//          f8_Dj = p_val5 - ( (p_val5 - 1.0) * (2.0 + expval6) / trm8 );
//          Cf8j = ( (1.0 - p_val5) / SQR(trm8) ) *( p_val6 * expval6 * trm8 -
//              (2.0 + expval6) * ( p_val6*expval6 - p_val7*expval7 ) );
//
//#ifdef SW_MATH
//          theta_0 = 180.0 - theta_00 * (1.0 - exp_4_sw(-p_val10 * (2.0-SBO2)));
//#else
//          theta_0 = 180.0 - theta_00 * (1.0 - exp(-p_val10 * (2.0 - SBO2)));
//#endif
//
//          theta_0 = DEG2RAD( theta_0 );
//
//#ifdef SW_MATH
//          expval2theta  = exp_4_sw( -p_val2 * SQR(theta_0 - p_ijk->theta) );
//#else
//          expval2theta  = exp( -p_val2 * SQR(theta_0 - p_ijk->theta) );
//#endif
//
//          if( p_val1 >= 0 )
//            expval12theta = p_val1 * (1.0 - expval2theta);
//          else // To avoid linear Me-H-Me angles (6/6/06)
//            expval12theta = p_val1 * -expval2theta;
//
//          CEval1 = Cf7ij * f7_jk * f8_Dj * expval12theta;
//          CEval2 = Cf7jk * f7_ij * f8_Dj * expval12theta;
//          CEval3 = Cf8j  * f7_ij * f7_jk * expval12theta;
//          CEval4 = -2.0 * p_val1 * p_val2 * f7_ij * f7_jk * f8_Dj *
//            expval2theta * (theta_0 - p_ijk->theta);
//
//#ifdef SW_MATH
//          Ctheta_0 = p_val10 * DEG2RAD(theta_00) * exp_4_sw(-p_val10 * (2.0-SBO2));
//#else
//          Ctheta_0 = p_val10 * DEG2RAD(theta_00) * exp(-p_val10 * (2.0-SBO2));
//#endif
//
//          CEval5 = -CEval4 * Ctheta_0 * CSBO2;
//          CEval6 = CEval5 * dSBO1;
//          CEval7 = CEval5 * dSBO2;
//          CEval8 = -CEval4 / sin_theta;
//
//          data->my_en.e_ang += e_ang = f7_ij * f7_jk * f8_Dj * expval12theta;
//          /* END ANGLE ENERGY*/
//
//          /* PENALTY ENERGY */
//          p_pen1 = thbp->p_pen1;
//
//#ifdef SW_MATH
//          exp_pen2ij = exp_4_sw( -p_pen2 * SQR( BOA_ij - 2.0 ) );
//          exp_pen2jk = exp_4_sw( -p_pen2 * SQR( BOA_jk - 2.0 ) );
//          exp_pen3   = exp_4_sw( -p_pen3 * workspace->Delta[j] );
//          exp_pen4   = exp_4_sw(  p_pen4 * workspace->Delta[j] );
//          
//#else
//          exp_pen2ij = exp( -p_pen2 * SQR( BOA_ij - 2.0 ) );
//          exp_pen2jk = exp( -p_pen2 * SQR( BOA_jk - 2.0 ) );
//          exp_pen3 = exp( -p_pen3 * workspace->Delta[j] );
//          exp_pen4 = exp(  p_pen4 * workspace->Delta[j] );
//#endif
//
//          trm_pen34 = 1.0 + exp_pen3 + exp_pen4;
//          f9_Dj = ( 2.0 + exp_pen3 ) / trm_pen34;
//          Cf9j = ( -p_pen3 * exp_pen3 * trm_pen34 -
//                   (2.0 + exp_pen3) * ( -p_pen3 * exp_pen3 +
//                  p_pen4 * exp_pen4 ) ) / SQR( trm_pen34 );
//
//          data->my_en.e_pen += e_pen = p_pen1 *f9_Dj *exp_pen2ij *exp_pen2jk;
//
//          CEpen1 = e_pen * Cf9j / f9_Dj;
//          temp   = -2.0 * p_pen2 * e_pen;
//          CEpen2 = temp * (BOA_ij - 2.0);
//          CEpen3 = temp * (BOA_jk - 2.0);
//          /* END PENALTY ENERGY */
//
//          /* COALITION ENERGY */
//          p_coa1 = thbp->p_coa1;
//
//#ifdef SW_MATH
//          exp_coa2 = exp_4_sw( p_coa2 * workspace->Delta_val[j] );
//          data->my_en.e_coa += e_coa = p_coa1 / (1. + exp_coa2) *
//            exp_4_sw( -p_coa3 * SQR(workspace->total_bond_order[i]-BOA_ij) ) *
//            exp_4_sw( -p_coa3 * SQR(workspace->total_bond_order[k]-BOA_jk) ) *
//            exp_4_sw( -p_coa4 * SQR(BOA_ij - 1.5) ) *
//            exp_4_sw( -p_coa4 * SQR(BOA_jk - 1.5) );
//
//#else
//          exp_coa2 = exp( p_coa2 * workspace->Delta_val[j] );
//          data->my_en.e_coa += e_coa = p_coa1 / (1. + exp_coa2) *
//            exp( -p_coa3 * SQR(workspace->total_bond_order[i]-BOA_ij) ) *
//            exp( -p_coa3 * SQR(workspace->total_bond_order[k]-BOA_jk) ) *
//            exp( -p_coa4 * SQR(BOA_ij - 1.5) ) *
//            exp( -p_coa4 * SQR(BOA_jk - 1.5) );
//
//#endif
//
//          CEcoa1 = -2 * p_coa4 * (BOA_ij - 1.5) * e_coa;
//          CEcoa2 = -2 * p_coa4 * (BOA_jk - 1.5) * e_coa;
//          CEcoa3 = -p_coa2 * exp_coa2 * e_coa / (1 + exp_coa2);
//          CEcoa4 = -2 * p_coa3 * (workspace->total_bond_order[i]-BOA_ij) *e_coa;
//          CEcoa5 = -2 * p_coa3 * (workspace->total_bond_order[k]-BOA_jk) * e_coa;
//
//          /* END COALITION ENERGY */
//
//          /* FORCES */
//          bo_ij->Cdbo += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));
//          bo_jk->Cdbo += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));
//          workspace->CdDelta[j] += ((CEval3 + CEval7) + CEpen1 + CEcoa3);
//          workspace->CdDelta[i] += CEcoa4;
//          workspace->CdDelta[k] += CEcoa5;
//
//          for( t = start_j; t < end_j; ++t ) 
//          {
//            pbond_jt = &( bonds->select.bond_list[t] );
//            bo_jt = &(pbond_jt->bo_data);
//            temp_bo_jt = bo_jt->BO;
//            temp = CUBE( temp_bo_jt );
//            pBOjt7 = temp * temp * temp_bo_jt;
//
//            bo_jt->Cdbo += (CEval6 * pBOjt7);
//            bo_jt->Cdbopi += CEval5;
//            bo_jt->Cdbopi2 += CEval5;
//          }//for-t
//
//          //////////// control->virial == 0/////////////
//          workspace->f[i][0] += CEval8 * (p_ijk->dcos_di)[0];
//          workspace->f[i][1] += CEval8 * (p_ijk->dcos_di)[1];
//          workspace->f[i][2] += CEval8 * (p_ijk->dcos_di)[2];
//
//          workspace->f[j][0] += CEval8 * (p_ijk->dcos_dj)[0];
//          workspace->f[j][1] += CEval8 * (p_ijk->dcos_dj)[1];
//          workspace->f[j][2] += CEval8 * (p_ijk->dcos_dj)[2];
//
//          workspace->f[k][0] += CEval8 * (p_ijk->dcos_dk)[0];
//          workspace->f[k][1] += CEval8 * (p_ijk->dcos_dk)[1];
//          workspace->f[k][2] += CEval8 * (p_ijk->dcos_dk)[2];
//          
//          /* tally into per-atom virials */
//          if( system->vflag_atom || system->evflag) 
//          {
//            delij[0] = system->my_atoms[i].x[0] - system->my_atoms[j].x[0];
//            delij[1] = system->my_atoms[i].x[1] - system->my_atoms[j].x[1];
//            delij[2] = system->my_atoms[i].x[2] - system->my_atoms[j].x[2];
//            
//            delkj[0] = system->my_atoms[k].x[0] - system->my_atoms[j].x[0];
//            delkj[1] = system->my_atoms[k].x[1] - system->my_atoms[j].x[1];
//            delkj[2] = system->my_atoms[k].x[2] - system->my_atoms[j].x[2];
//            
//            
//            fi_tmp[0] = -CEval8, p_ijk->dcos_di[0];
//            fi_tmp[1] = -CEval8, p_ijk->dcos_di[1];
//            fi_tmp[2] = -CEval8, p_ijk->dcos_di[2];
//
//            fj_tmp[0] = -CEval8, p_ijk->dcos_dj[0];
//            fj_tmp[1] = -CEval8, p_ijk->dcos_dj[1];
//            fj_tmp[2] = -CEval8, p_ijk->dcos_dj[2];
//
//            fk_tmp[0] = -CEval8, p_ijk->dcos_dk[0];
//            fk_tmp[1] = -CEval8, p_ijk->dcos_dk[1];
//            fk_tmp[2] = -CEval8, p_ijk->dcos_dk[2];
//
//            eng_tmp = e_ang + e_pen + e_coa;
//            if (system->evflag)
//            {
//              if (system->eflag_global){
//                system->eng_vdwl += eng_tmp;
//              }
//              if (system->eflag_atom){
//                system->eatom[i] += eng_tmp;
//              }
//
//              if (system->vflag_atom)
//              {
//                double v[6];
//                v[0] = THIRD * (delij[0]*fi_tmp[0] + delkj[0]*fj_tmp[0]);
//                v[1] = THIRD * (delij[1]*fi_tmp[1] + delkj[1]*fj_tmp[1]);
//                v[2] = THIRD * (delij[2]*fi_tmp[2] + delkj[2]*fj_tmp[2]);
//                v[3] = THIRD * (delij[0]*fi_tmp[1] + delkj[0]*fj_tmp[1]);
//                v[4] = THIRD * (delij[0]*fi_tmp[2] + delkj[0]*fj_tmp[2]);
//                v[5] = THIRD * (delij[1]*fi_tmp[2] + delkj[1]*fj_tmp[2]);
//            
//
//                system->vatom[i][0] += v[0]; 
//                system->vatom[i][1] += v[1]; 
//                system->vatom[i][2] += v[2];
//                system->vatom[i][3] += v[3]; 
//                system->vatom[i][4] += v[4]; 
//                system->vatom[i][5] += v[5];
//
//                system->vatom[j][0] += v[0]; 
//                system->vatom[j][1] += v[1]; 
//                system->vatom[j][2] += v[2];
//                system->vatom[j][3] += v[3]; 
//                system->vatom[j][4] += v[4]; 
//                system->vatom[j][5] += v[5];
//
//                system->vatom[k][0] += v[0]; 
//                system->vatom[k][1] += v[1]; 
//                system->vatom[k][2] += v[2];
//                system->vatom[k][3] += v[3]; 
//                system->vatom[k][4] += v[4]; 
//                system->vatom[k][5] += v[5];
//
//              }
//            }
//          }//if
//        }//for-cnt
//      }//for-pk
//    }//for-pi
//  }//for-j
//
//}


#endif


#ifdef CPE
#include "STUBS/mpi.h"
#include "slave.h"
#define THIRD 0.333333333333333333

#endif



