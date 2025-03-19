//
// File: multiModeQuad_ROS.h
//
// Code generated for Simulink model 'multiModeQuad_ROS'.
//
// Model version                  : 1.67
// Simulink Coder version         : 9.6 (R2021b) 14-May-2021
// C/C++ source code generated on : Wed Mar 19 11:44:06 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: Generic->Unspecified (assume 32-bit Generic)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#ifndef RTW_HEADER_multiModeQuad_ROS_h_
#define RTW_HEADER_multiModeQuad_ROS_h_
#include <string.h>
#include <math.h>
#include <stddef.h>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "slros_initialize.h"
#include "multiModeQuad_ROS_types.h"
#include "rtGetNaN.h"
#include "rt_nonfinite.h"
#include "rt_assert.h"
#include "rtGetInf.h"

// Macros for accessing real-time model data structure
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

// Block signals (default storage)
struct B_multiModeQuad_ROS_T {
  SL_Bus_multiModeQuad_ROS_sensor_msgs_Imu BusAssignment1;// '<Root>/Bus Assignment1' 
  creal_T V[16];
  creal_T At[16];
  SL_Bus_multiModeQuad_ROS_geometry_msgs_PoseStamped In1;// '<S14>/In1'
  SL_Bus_multiModeQuad_ROS_geometry_msgs_PoseStamped BusAssignment;// '<Root>/Bus Assignment' 
  SL_Bus_multiModeQuad_ROS_geometry_msgs_TwistStamped In1_o;// '<S15>/In1'
  SL_Bus_multiModeQuad_ROS_geometry_msgs_TwistStamped b_varargout_2;
  real_T VectorConcatenate[18];        // '<S36>/Vector Concatenate'
  real_T A[16];
  real_T T[16];
  real_T U[16];
  real_T Selector1[9];                 // '<S35>/Selector1'
  real_T Selector2[9];                 // '<S35>/Selector2'
  real_T VectorConcatenate_m[9];
  real_T Product_tmp[9];
  real_T Product_tmp_c[9];
  real_T A_k[9];
  creal_T D[4];
  creal_T y[4];
  creal_T work1[4];
  SL_Bus_multiModeQuad_ROS_geometry_msgs_Twist In1_h;// '<S17>/In1'
  real_T Selector[9];                  // '<S35>/Selector'
  SL_Bus_multiModeQuad_ROS_geometry_msgs_Twist b_varargout_2_c;
  char_T b_zeroDelimTopic[36];
  real_T feedback_quat[4];
  real_T quat_diff[4];
  real_T work[4];
  real_T rworka[4];
  char_T b_zeroDelimTopic_b[27];
  real_T Product[3];                   // '<S41>/Product'
  real_T Sum[3];                       // '<S18>/Sum'
  real_T TmpSignalConversionAtphithetaps[3];// '<S34>/phidot thetadot psidot'
  real_T Merge[4];                     // '<S20>/Merge'
  real_T att_sp_Body_temp[3];
  real_T sincos_o1_l[3];               // '<S43>/sincos'
  real_T ubvbwb[3];                    // '<S18>/ub,vb,wb'
  real_T att_sp_Body[3];               // '<S10>/R_EW1'
  real_T sincos_o1[3];                 // '<S31>/sincos'
  real_T rtb_Saturation2_p[3];
  real_T v[3];
  real_T v_c[3];
  creal_T tmp;
  creal_T ctemp;
  creal_T shift;
  creal_T ascale;
  int32_T rscale[4];
  boolean_T Compare[9];                // '<S87>/Compare'
  real_T FilterCoefficient;            // '<S124>/Filter Coefficient'
  real_T FilterCoefficient_m;          // '<S172>/Filter Coefficient'
  real_T Product2[3];                  // '<S35>/Product2'
  real_T Reshape[3];                   // '<S21>/Reshape'
  real_T ProportionalGain;             // '<S359>/Proportional Gain'
  real_T ProportionalGain_o;           // '<S311>/Proportional Gain'
  real_T ProportionalGain_l;           // '<S263>/Proportional Gain'
  real_T FilterCoefficient_j;          // '<S220>/Filter Coefficient'
  real_T Komega;                       // '<S10>/Komega'
  real_T Komega1;                      // '<S10>/Komega1'
  real_T Komega2;                      // '<S10>/Komega2'
  real_T Komega3;                      // '<S10>/Komega3'
  real_T IntegralGain;                 // '<S118>/Integral Gain'
  real_T IntegralGain_i;               // '<S166>/Integral Gain'
  real_T IntegralGain_g;               // '<S214>/Integral Gain'
  real_T K12;
  real_T K14;
  real_T K23;
  real_T K24;
  real_T K34;
  real_T rtb_Switch_idx_0;
  real_T rtb_Switch_idx_1;
  real_T rtb_Switch_idx_3;
  real_T rtb_Sum2_idx_0;
  real_T rtb_Sum2_idx_1;
  real_T rtb_Sum2_idx_2;
  real_T rtb_sincos_o2_idx_0;
  real_T Selector_tmp;
  real_T K12_f;
  real_T FilterCoefficient_i;          // '<S268>/Filter Coefficient'
  real_T FilterCoefficient_a;          // '<S316>/Filter Coefficient'
  real_T FilterCoefficient_g;          // '<S364>/Filter Coefficient'
  real_T IntegralGain_g1;              // '<S262>/Integral Gain'
  real_T IntegralGain_i_m;             // '<S310>/Integral Gain'
  real_T IntegralGain_j;               // '<S358>/Integral Gain'
  real_T anrm;
  real_T vtemp;
  real_T b_absxk;
  real_T cfromc;
  real_T ctoc;
  real_T cto1;
  real_T stemp_im_tmp;
  real_T anorm;
  real_T scale;
  real_T ssq;
  real_T colscale;
  real_T absxk;
  real_T t;
  real_T ar;
  real_T ai;
  real_T t1_re;
  real_T t1_im;
  real_T shift_im;
  real_T eshift_re;
  real_T eshift_im;
  real_T shift_tmp;
  real_T g2;
  real_T f2s;
  real_T di;
  real_T scale_n;
  real_T x;
  real_T fs_re;
  real_T fs_im;
  real_T gs_re;
  real_T gs_im;
  real_T a;
  real_T alpha1;
  real_T beta1;
  real_T tau_idx_2;
  real_T tau_idx_1;
  real_T tau_idx_0;
  real_T xnorm_tmp_tmp;
  real_T tst;
  real_T htmp1;
  real_T ab;
  real_T ba;
  real_T aa;
  real_T h12;
  real_T h21s;
  real_T a__4;
  real_T s_tmp;
  real_T p;
  real_T bcmax;
  real_T bcmis;
  real_T scale_p;
  real_T z;
  real_T tau;
  real_T anorm_l;
  real_T ascale_j;
  real_T temp;
  real_T acoeff;
  real_T scale_d;
  real_T dmin;
  real_T e_y;
  real_T salpha_re;
  real_T salpha_im;
  real_T work2_idx_2_im;
  real_T work2_idx_3_re;
  real_T g2_g;
  real_T d;
  real_T f2s_l;
  real_T scale_dh;
  real_T x_d;
  real_T fs_re_l;
  real_T fs_im_o;
  real_T gs_re_b;
  real_T xnorm;
  real_T a_n;
  real_T scale_b;
  real_T absxk_l;
  real_T t_h;
  real_T temp_b;
  real_T absxr;
  real_T temp_d;
  SL_Bus_multiModeQuad_ROS_std_msgs_Int16 In1_e;// '<S13>/In1'
  SL_Bus_multiModeQuad_ROS_std_msgs_Float32 In1_ol;// '<S16>/In1'
  real32_T Gain;                       // '<S10>/Gain'
  int32_T c_i;
  int32_T c_j;
  int32_T sgn;
  int32_T rtb_VectorConcatenate_tmp;
  int32_T b_k;
  int32_T jrow;
  int32_T jcol;
  int32_T b_j;
  int32_T c_i_e;
  int32_T stemp_re_tmp;
  int32_T stemp_re_tmp_b;
  int32_T j;
  int32_T ifirst;
  int32_T knt;
  int32_T lastc;
  int32_T ix;
  int32_T iac;
  int32_T g;
  int32_T b_ia;
  int32_T jy;
  int32_T b_ix;
  int32_T i;
  int32_T L;
  SL_Bus_multiModeQuad_ROS_std_msgs_Float32 b_varargout_2_j;
};

// Block states (default storage) for system '<Root>'
struct DW_multiModeQuad_ROS_T {
  ros_slros_internal_block_GetP_T obj; // '<Root>/Get Parameter2'
  ros_slros_internal_block_GetP_T obj_a;// '<Root>/Get Parameter1'
  ros_slros_internal_block_GetP_T obj_p;// '<Root>/Get Parameter'
  ros_slroscpp_internal_block_S_T obj_c;// '<S9>/SourceBlock'
  ros_slroscpp_internal_block_S_T obj_ps;// '<S8>/SourceBlock'
  ros_slroscpp_internal_block_S_T obj_b;// '<S7>/SourceBlock'
  ros_slroscpp_internal_block_S_T obj_ap;// '<S6>/SourceBlock'
  ros_slroscpp_internal_block_S_T obj_l;// '<S3>/SourceBlock'
  ros_slroscpp_internal_block_P_T obj_k;// '<S5>/SinkBlock'
  ros_slroscpp_internal_block_P_T obj_m;// '<S4>/SinkBlock'
  real_T Filter_DSTATE;                // '<S356>/Filter'
  real_T Integrator_DSTATE;            // '<S361>/Integrator'
  real_T Integrator_DSTATE_f;          // '<S313>/Integrator'
  real_T Filter_DSTATE_a;              // '<S308>/Filter'
  real_T Integrator_DSTATE_l;          // '<S265>/Integrator'
  real_T Filter_DSTATE_c;              // '<S260>/Filter'
  real_T Product2_DWORK4[9];           // '<S35>/Product2'
  int8_T If_ActiveSubsystem;           // '<S20>/If'
  int8_T If1_ActiveSubsystem;          // '<S54>/If1'
  int8_T FindMaximumDiagonalValue_Active;// '<S52>/Find Maximum Diagonal Value'
};

// Continuous states (default storage)
struct X_multiModeQuad_ROS_T {
  real_T phithetapsi_CSTATE[3];        // '<S34>/phi theta psi'
  real_T ubvbwb_CSTATE[3];             // '<S18>/ub,vb,wb'
  real_T Integrator_CSTATE;            // '<S121>/Integrator'
  real_T Filter_CSTATE;                // '<S116>/Filter'
  real_T Integrator_CSTATE_o;          // '<S169>/Integrator'
  real_T Filter_CSTATE_k;              // '<S164>/Filter'
  real_T pqr_CSTATE[3];                // '<S18>/p,q,r '
  real_T TransferFcn2_CSTATE;          // '<S10>/Transfer Fcn2'
  real_T TransferFcn1_CSTATE;          // '<S10>/Transfer Fcn1'
  real_T TransferFcn3_CSTATE;          // '<S10>/Transfer Fcn3'
  real_T TransferFcn4_CSTATE;          // '<S10>/Transfer Fcn4'
  real_T xeyeze_CSTATE[3];             // '<S18>/xe,ye,ze'
  real_T Integrator_CSTATE_h;          // '<S217>/Integrator'
  real_T Filter_CSTATE_h;              // '<S212>/Filter'
};

// Periodic continuous state vector (global)
typedef int_T PeriodicIndX_multiModeQuad_RO_T[3];
typedef real_T PeriodicRngX_multiModeQuad_RO_T[6];

// State derivatives (default storage)
struct XDot_multiModeQuad_ROS_T {
  real_T phithetapsi_CSTATE[3];        // '<S34>/phi theta psi'
  real_T ubvbwb_CSTATE[3];             // '<S18>/ub,vb,wb'
  real_T Integrator_CSTATE;            // '<S121>/Integrator'
  real_T Filter_CSTATE;                // '<S116>/Filter'
  real_T Integrator_CSTATE_o;          // '<S169>/Integrator'
  real_T Filter_CSTATE_k;              // '<S164>/Filter'
  real_T pqr_CSTATE[3];                // '<S18>/p,q,r '
  real_T TransferFcn2_CSTATE;          // '<S10>/Transfer Fcn2'
  real_T TransferFcn1_CSTATE;          // '<S10>/Transfer Fcn1'
  real_T TransferFcn3_CSTATE;          // '<S10>/Transfer Fcn3'
  real_T TransferFcn4_CSTATE;          // '<S10>/Transfer Fcn4'
  real_T xeyeze_CSTATE[3];             // '<S18>/xe,ye,ze'
  real_T Integrator_CSTATE_h;          // '<S217>/Integrator'
  real_T Filter_CSTATE_h;              // '<S212>/Filter'
};

// State disabled
struct XDis_multiModeQuad_ROS_T {
  boolean_T phithetapsi_CSTATE[3];     // '<S34>/phi theta psi'
  boolean_T ubvbwb_CSTATE[3];          // '<S18>/ub,vb,wb'
  boolean_T Integrator_CSTATE;         // '<S121>/Integrator'
  boolean_T Filter_CSTATE;             // '<S116>/Filter'
  boolean_T Integrator_CSTATE_o;       // '<S169>/Integrator'
  boolean_T Filter_CSTATE_k;           // '<S164>/Filter'
  boolean_T pqr_CSTATE[3];             // '<S18>/p,q,r '
  boolean_T TransferFcn2_CSTATE;       // '<S10>/Transfer Fcn2'
  boolean_T TransferFcn1_CSTATE;       // '<S10>/Transfer Fcn1'
  boolean_T TransferFcn3_CSTATE;       // '<S10>/Transfer Fcn3'
  boolean_T TransferFcn4_CSTATE;       // '<S10>/Transfer Fcn4'
  boolean_T xeyeze_CSTATE[3];          // '<S18>/xe,ye,ze'
  boolean_T Integrator_CSTATE_h;       // '<S217>/Integrator'
  boolean_T Filter_CSTATE_h;           // '<S212>/Filter'
};

#ifndef ODE3_INTG
#define ODE3_INTG

// ODE3 Integration Data
struct ODE3_IntgData {
  real_T *y;                           // output
  real_T *f[3];                        // derivatives
};

#endif

// Parameters (default storage)
struct P_multiModeQuad_ROS_T_ {
  real_T PIDVelocityx_D;               // Mask Parameter: PIDVelocityx_D
                                          //  Referenced by: '<S115>/Derivative Gain'

  real_T PIDVelocityy_D;               // Mask Parameter: PIDVelocityy_D
                                          //  Referenced by: '<S163>/Derivative Gain'

  real_T PIDangularroll_D;             // Mask Parameter: PIDangularroll_D
                                          //  Referenced by: '<S355>/Derivative Gain'

  real_T PIDangularpitch_D;            // Mask Parameter: PIDangularpitch_D
                                          //  Referenced by: '<S307>/Derivative Gain'

  real_T PIDangulayaw_D;               // Mask Parameter: PIDangulayaw_D
                                          //  Referenced by: '<S259>/Derivative Gain'

  real_T PIDVelocityz_D;               // Mask Parameter: PIDVelocityz_D
                                          //  Referenced by: '<S211>/Derivative Gain'

  real_T PIDangulayaw_I;               // Mask Parameter: PIDangulayaw_I
                                          //  Referenced by: '<S262>/Integral Gain'

  real_T PIDangularpitch_I;            // Mask Parameter: PIDangularpitch_I
                                          //  Referenced by: '<S310>/Integral Gain'

  real_T PIDangularroll_I;             // Mask Parameter: PIDangularroll_I
                                          //  Referenced by: '<S358>/Integral Gain'

  real_T PIDVelocityx_I;               // Mask Parameter: PIDVelocityx_I
                                          //  Referenced by: '<S118>/Integral Gain'

  real_T PIDVelocityy_I;               // Mask Parameter: PIDVelocityy_I
                                          //  Referenced by: '<S166>/Integral Gain'

  real_T PIDVelocityz_I;               // Mask Parameter: PIDVelocityz_I
                                          //  Referenced by: '<S214>/Integral Gain'

  real_T PIDVelocityx_InitialConditionFo;
                              // Mask Parameter: PIDVelocityx_InitialConditionFo
                                 //  Referenced by: '<S116>/Filter'

  real_T PIDVelocityy_InitialConditionFo;
                              // Mask Parameter: PIDVelocityy_InitialConditionFo
                                 //  Referenced by: '<S164>/Filter'

  real_T PIDangularroll_InitialCondition;
                              // Mask Parameter: PIDangularroll_InitialCondition
                                 //  Referenced by: '<S356>/Filter'

  real_T PIDangularpitch_InitialConditio;
                              // Mask Parameter: PIDangularpitch_InitialConditio
                                 //  Referenced by: '<S308>/Filter'

  real_T PIDangulayaw_InitialConditionFo;
                              // Mask Parameter: PIDangulayaw_InitialConditionFo
                                 //  Referenced by: '<S260>/Filter'

  real_T PIDVelocityz_InitialConditionFo;
                              // Mask Parameter: PIDVelocityz_InitialConditionFo
                                 //  Referenced by: '<S212>/Filter'

  real_T PIDVelocityx_InitialCondition_o;
                              // Mask Parameter: PIDVelocityx_InitialCondition_o
                                 //  Referenced by: '<S121>/Integrator'

  real_T PIDVelocityy_InitialCondition_i;
                              // Mask Parameter: PIDVelocityy_InitialCondition_i
                                 //  Referenced by: '<S169>/Integrator'

  real_T PIDangularroll_InitialConditi_o;
                              // Mask Parameter: PIDangularroll_InitialConditi_o
                                 //  Referenced by: '<S361>/Integrator'

  real_T PIDangularpitch_InitialCondit_b;
                              // Mask Parameter: PIDangularpitch_InitialCondit_b
                                 //  Referenced by: '<S313>/Integrator'

  real_T PIDangulayaw_InitialCondition_m;
                              // Mask Parameter: PIDangulayaw_InitialCondition_m
                                 //  Referenced by: '<S265>/Integrator'

  real_T PIDVelocityz_InitialCondition_g;
                              // Mask Parameter: PIDVelocityz_InitialCondition_g
                                 //  Referenced by: '<S217>/Integrator'

  real_T PIDVelocityx_N;               // Mask Parameter: PIDVelocityx_N
                                          //  Referenced by: '<S124>/Filter Coefficient'

  real_T PIDVelocityy_N;               // Mask Parameter: PIDVelocityy_N
                                          //  Referenced by: '<S172>/Filter Coefficient'

  real_T PIDangularroll_N;             // Mask Parameter: PIDangularroll_N
                                          //  Referenced by: '<S364>/Filter Coefficient'

  real_T PIDangularpitch_N;            // Mask Parameter: PIDangularpitch_N
                                          //  Referenced by: '<S316>/Filter Coefficient'

  real_T PIDangulayaw_N;               // Mask Parameter: PIDangulayaw_N
                                          //  Referenced by: '<S268>/Filter Coefficient'

  real_T PIDVelocityz_N;               // Mask Parameter: PIDVelocityz_N
                                          //  Referenced by: '<S220>/Filter Coefficient'

  real_T PIDVelocityx_P;               // Mask Parameter: PIDVelocityx_P
                                          //  Referenced by: '<S126>/Proportional Gain'

  real_T PIDVelocityy_P;               // Mask Parameter: PIDVelocityy_P
                                          //  Referenced by: '<S174>/Proportional Gain'

  real_T PIDangularroll_P;             // Mask Parameter: PIDangularroll_P
                                          //  Referenced by: '<S359>/Proportional Gain'

  real_T PIDangularpitch_P;            // Mask Parameter: PIDangularpitch_P
                                          //  Referenced by: '<S311>/Proportional Gain'

  real_T PIDangulayaw_P;               // Mask Parameter: PIDangulayaw_P
                                          //  Referenced by: '<S263>/Proportional Gain'

  real_T PIDVelocityz_P;               // Mask Parameter: PIDVelocityz_P
                                          //  Referenced by: '<S222>/Proportional Gain'

  real_T uDOFEulerAngles2_Vm_0[3];     // Mask Parameter: uDOFEulerAngles2_Vm_0
                                          //  Referenced by: '<S18>/ub,vb,wb'

  real_T DirectionCosineMatrixtoQuaterni;
                              // Mask Parameter: DirectionCosineMatrixtoQuaterni
                                 //  Referenced by:
                                 //    '<S54>/Constant'
                                 //    '<S79>/Constant'
                                 //    '<S80>/Constant'

  real_T uDOFEulerAngles2_eul_0[3];    // Mask Parameter: uDOFEulerAngles2_eul_0
                                          //  Referenced by: '<S34>/phi theta psi'

  real_T uDOFEulerAngles2_inertia[9];// Mask Parameter: uDOFEulerAngles2_inertia
                                        //  Referenced by: '<S36>/Constant1'

  real_T uDOFEulerAngles2_mass_0;     // Mask Parameter: uDOFEulerAngles2_mass_0
                                         //  Referenced by: '<S36>/Constant'

  real_T uDOFEulerAngles2_pm_0[3];     // Mask Parameter: uDOFEulerAngles2_pm_0
                                          //  Referenced by: '<S18>/p,q,r '

  real_T DirectionCosineMatrixtoQuater_h;
                              // Mask Parameter: DirectionCosineMatrixtoQuater_h
                                 //  Referenced by:
                                 //    '<S87>/Constant'
                                 //    '<S89>/Constant'

  real_T uDOFEulerAngles2_xme_0[3];    // Mask Parameter: uDOFEulerAngles2_xme_0
                                          //  Referenced by: '<S18>/xe,ye,ze'

  SL_Bus_multiModeQuad_ROS_sensor_msgs_Imu Constant_Value;// Computed Parameter: Constant_Value
                                                             //  Referenced by: '<S2>/Constant'

  SL_Bus_multiModeQuad_ROS_geometry_msgs_PoseStamped Constant_Value_f;// Computed Parameter: Constant_Value_f
                                                                      //  Referenced by: '<S1>/Constant'

  SL_Bus_multiModeQuad_ROS_geometry_msgs_PoseStamped Out1_Y0;// Computed Parameter: Out1_Y0
                                                                //  Referenced by: '<S14>/Out1'

  SL_Bus_multiModeQuad_ROS_geometry_msgs_PoseStamped Constant_Value_a;// Computed Parameter: Constant_Value_a
                                                                      //  Referenced by: '<S6>/Constant'

  SL_Bus_multiModeQuad_ROS_geometry_msgs_TwistStamped Out1_Y0_k;// Computed Parameter: Out1_Y0_k
                                                                   //  Referenced by: '<S15>/Out1'

  SL_Bus_multiModeQuad_ROS_geometry_msgs_TwistStamped Constant_Value_c;// Computed Parameter: Constant_Value_c
                                                                      //  Referenced by: '<S7>/Constant'

  SL_Bus_multiModeQuad_ROS_geometry_msgs_Twist Out1_Y0_g;// Computed Parameter: Out1_Y0_g
                                                            //  Referenced by: '<S17>/Out1'

  SL_Bus_multiModeQuad_ROS_geometry_msgs_Twist Constant_Value_fz;// Computed Parameter: Constant_Value_fz
                                                                    //  Referenced by: '<S9>/Constant'

  SL_Bus_multiModeQuad_ROS_std_msgs_Float32 Out1_Y0_kv;// Computed Parameter: Out1_Y0_kv
                                                          //  Referenced by: '<S16>/Out1'

  SL_Bus_multiModeQuad_ROS_std_msgs_Float32 Constant_Value_p;// Computed Parameter: Constant_Value_p
                                                                //  Referenced by: '<S8>/Constant'

  SL_Bus_multiModeQuad_ROS_std_msgs_Int16 Out1_Y0_p;// Computed Parameter: Out1_Y0_p
                                                       //  Referenced by: '<S13>/Out1'

  SL_Bus_multiModeQuad_ROS_std_msgs_Int16 Constant_Value_e;// Computed Parameter: Constant_Value_e
                                                              //  Referenced by: '<S3>/Constant'

  real_T Constant_Value_o;             // Expression: 1
                                          //  Referenced by: '<S53>/Constant'

  real_T Gain_Gain;                    // Expression: 0.5
                                          //  Referenced by: '<S53>/Gain'

  real_T Gain1_Gain;                   // Expression: 2
                                          //  Referenced by: '<S53>/Gain1'

  real_T Constant_Value_n;             // Expression: 1
                                          //  Referenced by: '<S69>/Constant'

  real_T Constant1_Value;              // Expression: 0.5
                                          //  Referenced by: '<S68>/Constant1'

  real_T Constant2_Value[2];           // Expression: [0 1]
                                          //  Referenced by: '<S68>/Constant2'

  real_T Gain1_Gain_k;                 // Expression: 1
                                          //  Referenced by: '<S57>/Gain1'

  real_T Gain3_Gain;                   // Expression: 1
                                          //  Referenced by: '<S57>/Gain3'

  real_T Gain4_Gain;                   // Expression: 1
                                          //  Referenced by: '<S57>/Gain4'

  real_T Gain_Gain_g;                  // Expression: 0.5
                                          //  Referenced by: '<S57>/Gain'

  real_T Constant_Value_pu;            // Expression: 1
                                          //  Referenced by: '<S74>/Constant'

  real_T Constant1_Value_h;            // Expression: 0.5
                                          //  Referenced by: '<S73>/Constant1'

  real_T Constant2_Value_k[2];         // Expression: [0 1]
                                          //  Referenced by: '<S73>/Constant2'

  real_T Gain1_Gain_l;                 // Expression: 1
                                          //  Referenced by: '<S58>/Gain1'

  real_T Gain2_Gain;                   // Expression: 1
                                          //  Referenced by: '<S58>/Gain2'

  real_T Gain3_Gain_a;                 // Expression: 1
                                          //  Referenced by: '<S58>/Gain3'

  real_T Gain_Gain_d;                  // Expression: 0.5
                                          //  Referenced by: '<S58>/Gain'

  real_T Constant_Value_fv;            // Expression: 1
                                          //  Referenced by: '<S64>/Constant'

  real_T Constant1_Value_b;            // Expression: 0.5
                                          //  Referenced by: '<S63>/Constant1'

  real_T Constant2_Value_o[2];         // Expression: [0 1]
                                          //  Referenced by: '<S63>/Constant2'

  real_T Gain1_Gain_b;                 // Expression: 1
                                          //  Referenced by: '<S56>/Gain1'

  real_T Gain2_Gain_d;                 // Expression: 1
                                          //  Referenced by: '<S56>/Gain2'

  real_T Gain3_Gain_az;                // Expression: 1
                                          //  Referenced by: '<S56>/Gain3'

  real_T Gain_Gain_i;                  // Expression: 0.5
                                          //  Referenced by: '<S56>/Gain'

  real_T Constant1_Value_d;            // Expression: 0
                                          //  Referenced by: '<S80>/Constant1'

  real_T Constant1_Value_j;            // Expression: 0
                                          //  Referenced by: '<S79>/Constant1'

  real_T Bias1_Bias[9];                // Expression: -eye(3)
                                          //  Referenced by: '<S81>/Bias1'

  real_T Bias_Bias;                    // Expression: -1
                                          //  Referenced by: '<S82>/Bias'

  real_T u2_Gain;                      // Expression: 0.5
                                          //  Referenced by: '<S31>/1//2'

  real_T Constant2_Value_n[9];         // Expression: zeros(3)
                                          //  Referenced by: '<S36>/Constant2'

  real_T phithetapsi_WrappedStateUpperVa;// Expression: pi
                                            //  Referenced by: '<S34>/phi theta psi'

  real_T phithetapsi_WrappedStateLowerVa;// Expression: -pi
                                            //  Referenced by: '<S34>/phi theta psi'

  real_T Gain_Gain_ip;                 // Expression: -1
                                          //  Referenced by: '<S32>/Gain'

  real_T TransferFcn2_A;               // Computed Parameter: TransferFcn2_A
                                          //  Referenced by: '<S10>/Transfer Fcn2'

  real_T TransferFcn2_C;               // Computed Parameter: TransferFcn2_C
                                          //  Referenced by: '<S10>/Transfer Fcn2'

  real_T Saturation_UpperSat;          // Expression: 10
                                          //  Referenced by: '<S10>/Saturation'

  real_T Saturation_LowerSat;          // Expression: 0
                                          //  Referenced by: '<S10>/Saturation'

  real_T TransferFcn1_A;               // Computed Parameter: TransferFcn1_A
                                          //  Referenced by: '<S10>/Transfer Fcn1'

  real_T TransferFcn1_C;               // Computed Parameter: TransferFcn1_C
                                          //  Referenced by: '<S10>/Transfer Fcn1'

  real_T Saturation1_UpperSat;         // Expression: 10
                                          //  Referenced by: '<S10>/Saturation1'

  real_T Saturation1_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S10>/Saturation1'

  real_T TransferFcn3_A;               // Computed Parameter: TransferFcn3_A
                                          //  Referenced by: '<S10>/Transfer Fcn3'

  real_T TransferFcn3_C;               // Computed Parameter: TransferFcn3_C
                                          //  Referenced by: '<S10>/Transfer Fcn3'

  real_T Saturation2_UpperSat;         // Expression: 10
                                          //  Referenced by: '<S10>/Saturation2'

  real_T Saturation2_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S10>/Saturation2'

  real_T TransferFcn4_A;               // Computed Parameter: TransferFcn4_A
                                          //  Referenced by: '<S10>/Transfer Fcn4'

  real_T TransferFcn4_C;               // Computed Parameter: TransferFcn4_C
                                          //  Referenced by: '<S10>/Transfer Fcn4'

  real_T Saturation3_UpperSat;         // Expression: 10
                                          //  Referenced by: '<S10>/Saturation3'

  real_T Saturation3_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S10>/Saturation3'

  real_T Constant4_Value;              // Expression: 0
                                          //  Referenced by: '<S21>/Constant4'

  real_T Constant5_Value;              // Expression: 0
                                          //  Referenced by: '<S21>/Constant5'

  real_T Mass_Value;                   // Expression: 2
                                          //  Referenced by: '<S21>/Mass'

  real_T g_Gain;                       // Expression: 9.8
                                          //  Referenced by: '<S21>/g'

  real_T Gain_Gain_m;                  // Expression: -1
                                          //  Referenced by: '<S33>/Gain'

  real_T Merge_InitialOutput[4];       // Expression: [1 0 0 0]
                                          //  Referenced by: '<S20>/Merge'

  real_T Filter_gainval;               // Computed Parameter: Filter_gainval
                                          //  Referenced by: '<S356>/Filter'

  real_T Integrator_gainval;           // Computed Parameter: Integrator_gainval
                                          //  Referenced by: '<S361>/Integrator'

  real_T Integrator_gainval_e;       // Computed Parameter: Integrator_gainval_e
                                        //  Referenced by: '<S313>/Integrator'

  real_T Filter_gainval_j;             // Computed Parameter: Filter_gainval_j
                                          //  Referenced by: '<S308>/Filter'

  real_T Integrator_gainval_a;       // Computed Parameter: Integrator_gainval_a
                                        //  Referenced by: '<S265>/Integrator'

  real_T Filter_gainval_c;             // Computed Parameter: Filter_gainval_c
                                          //  Referenced by: '<S260>/Filter'

  real_T Hoverthrottle_Value;          // Expression: 3.13
                                          //  Referenced by: '<S10>/Hover throttle'

  real_T Saturation4_UpperSat;         // Expression: inf
                                          //  Referenced by: '<S10>/Saturation4'

  real_T Saturation4_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S10>/Saturation4'

  real_T Komega_Gain;                  // Expression: 0.5
                                          //  Referenced by: '<S10>/Komega'

  real_T Saturation5_UpperSat;         // Expression: inf
                                          //  Referenced by: '<S10>/Saturation5'

  real_T Saturation5_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S10>/Saturation5'

  real_T Komega1_Gain;                 // Expression: 0.5
                                          //  Referenced by: '<S10>/Komega1'

  real_T Saturation6_UpperSat;         // Expression: inf
                                          //  Referenced by: '<S10>/Saturation6'

  real_T Saturation6_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S10>/Saturation6'

  real_T Komega2_Gain;                 // Expression: 0.5
                                          //  Referenced by: '<S10>/Komega2'

  real_T Saturation7_UpperSat;         // Expression: inf
                                          //  Referenced by: '<S10>/Saturation7'

  real_T Saturation7_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S10>/Saturation7'

  real_T Komega3_Gain;                 // Expression: 0.5
                                          //  Referenced by: '<S10>/Komega3'

  real32_T Gain_Gain_dx;               // Computed Parameter: Gain_Gain_dx
                                          //  Referenced by: '<S10>/Gain'

  int16_T Switch3_Threshold;           // Computed Parameter: Switch3_Threshold
                                          //  Referenced by: '<S10>/Switch3'

  int16_T Switch_Threshold;            // Computed Parameter: Switch_Threshold
                                          //  Referenced by: '<S10>/Switch'

  int16_T Switch1_Threshold;           // Computed Parameter: Switch1_Threshold
                                          //  Referenced by: '<S10>/Switch1'

  int16_T Switch2_Threshold;           // Computed Parameter: Switch2_Threshold
                                          //  Referenced by: '<S10>/Switch2'

};

// Real-time Model Data Structure
struct tag_RTM_multiModeQuad_ROS_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_multiModeQuad_ROS_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[22];
  real_T odeF[3][22];
  ODE3_IntgData intgData;

  //
  //  Sizes:
  //  The following substructure contains sizes information
  //  for many of the model attributes such as inputs, outputs,
  //  dwork, sample times, etc.

  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  //
  //  Timing:
  //  The following substructure contains information regarding
  //  the timing information for the model.

  struct {
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

// Block parameters (default storage)
#ifdef __cplusplus

extern "C" {

#endif

  extern P_multiModeQuad_ROS_T multiModeQuad_ROS_P;

#ifdef __cplusplus

}
#endif

// Block signals (default storage)
#ifdef __cplusplus

extern "C" {

#endif

  extern struct B_multiModeQuad_ROS_T multiModeQuad_ROS_B;

#ifdef __cplusplus

}
#endif

// Continuous states (default storage)
extern X_multiModeQuad_ROS_T multiModeQuad_ROS_X;

// Block states (default storage)
extern struct DW_multiModeQuad_ROS_T multiModeQuad_ROS_DW;

#ifdef __cplusplus

extern "C" {

#endif

  // Model entry point functions
  extern void multiModeQuad_ROS_initialize(void);
  extern void multiModeQuad_ROS_step(void);
  extern void multiModeQuad_ROS_terminate(void);

#ifdef __cplusplus

}
#endif

// Real-time Model object
#ifdef __cplusplus

extern "C" {

#endif

  extern RT_MODEL_multiModeQuad_ROS_T *const multiModeQuad_ROS_M;

#ifdef __cplusplus

}
#endif

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S38>/Unit Conversion' : Unused code path elimination
//  Block '<S40>/Unit Conversion' : Unused code path elimination
//  Block '<S10>/Scope2' : Unused code path elimination
//  Block '<S10>/attitude' : Unused code path elimination
//  Block '<S10>/omega' : Unused code path elimination
//  Block '<S10>/pose' : Unused code path elimination
//  Block '<S10>/velocity' : Unused code path elimination
//  Block '<S10>/velocity1' : Unused code path elimination
//  Block '<S44>/Reshape (9) to [3x3] column-major' : Reshape block reduction
//  Block '<S46>/Reshape1' : Reshape block reduction
//  Block '<S46>/Reshape2' : Reshape block reduction
//  Block '<S47>/Reshape1' : Reshape block reduction
//  Block '<S47>/Reshape2' : Reshape block reduction
//  Block '<S35>/Reshape' : Reshape block reduction
//  Block '<S35>/Reshape1' : Reshape block reduction
//  Block '<S39>/Unit Conversion' : Eliminated nontunable gain of 1
//  Block '<S41>/Reshape1' : Reshape block reduction
//  Block '<S41>/Reshape2' : Reshape block reduction
//  Block '<S20>/Reshape 3x3 -> 9' : Reshape block reduction
//  Block '<S81>/Reshape' : Reshape block reduction


//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Use the MATLAB hilite_system command to trace the generated code back
//  to the model.  For example,
//
//  hilite_system('<S3>')    - opens system 3
//  hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'multiModeQuad_ROS'
//  '<S1>'   : 'multiModeQuad_ROS/Blank Message'
//  '<S2>'   : 'multiModeQuad_ROS/Blank Message1'
//  '<S3>'   : 'multiModeQuad_ROS/Flight mode'
//  '<S4>'   : 'multiModeQuad_ROS/Publish'
//  '<S5>'   : 'multiModeQuad_ROS/Publish1'
//  '<S6>'   : 'multiModeQuad_ROS/Sub setpoint attitude'
//  '<S7>'   : 'multiModeQuad_ROS/Sub setpoint rate'
//  '<S8>'   : 'multiModeQuad_ROS/Sub setpoint thrust'
//  '<S9>'   : 'multiModeQuad_ROS/Sub setpoint velocity'
//  '<S10>'  : 'multiModeQuad_ROS/Subsystem1'
//  '<S11>'  : 'multiModeQuad_ROS/time to sec & nsec'
//  '<S12>'  : 'multiModeQuad_ROS/time to sec & nsec1'
//  '<S13>'  : 'multiModeQuad_ROS/Flight mode/Enabled Subsystem'
//  '<S14>'  : 'multiModeQuad_ROS/Sub setpoint attitude/Enabled Subsystem'
//  '<S15>'  : 'multiModeQuad_ROS/Sub setpoint rate/Enabled Subsystem'
//  '<S16>'  : 'multiModeQuad_ROS/Sub setpoint thrust/Enabled Subsystem'
//  '<S17>'  : 'multiModeQuad_ROS/Sub setpoint velocity/Enabled Subsystem'
//  '<S18>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2'
//  '<S19>'  : 'multiModeQuad_ROS/Subsystem1/Attitude control'
//  '<S20>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions'
//  '<S21>'  : 'multiModeQuad_ROS/Subsystem1/Gravity1'
//  '<S22>'  : 'multiModeQuad_ROS/Subsystem1/Mapping'
//  '<S23>'  : 'multiModeQuad_ROS/Subsystem1/Mixer'
//  '<S24>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x'
//  '<S25>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity y'
//  '<S26>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity z'
//  '<S27>'  : 'multiModeQuad_ROS/Subsystem1/PID angula yaw'
//  '<S28>'  : 'multiModeQuad_ROS/Subsystem1/PID angular pitch'
//  '<S29>'  : 'multiModeQuad_ROS/Subsystem1/PID angular roll'
//  '<S30>'  : 'multiModeQuad_ROS/Subsystem1/R_EW1'
//  '<S31>'  : 'multiModeQuad_ROS/Subsystem1/Rotation Angles to Quaternions'
//  '<S32>'  : 'multiModeQuad_ROS/Subsystem1/Subsystem'
//  '<S33>'  : 'multiModeQuad_ROS/Subsystem1/Subsystem1'
//  '<S34>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Calculate DCM & Euler Angles'
//  '<S35>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot'
//  '<S36>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Determine Force,  Mass & Inertia'
//  '<S37>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Vbxw'
//  '<S38>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Velocity Conversion'
//  '<S39>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Velocity Conversion1'
//  '<S40>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Velocity Conversion2'
//  '<S41>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/transform to Inertial axes '
//  '<S42>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix'
//  '<S43>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Calculate DCM & Euler Angles/phidot thetadot psidot'
//  '<S44>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
//  '<S45>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot/3x3 Cross Product'
//  '<S46>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot/I x w'
//  '<S47>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot/I x w1'
//  '<S48>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot/3x3 Cross Product/Subsystem'
//  '<S49>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot/3x3 Cross Product/Subsystem1'
//  '<S50>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Vbxw/Subsystem'
//  '<S51>'  : 'multiModeQuad_ROS/Subsystem1/6DOF (Euler Angles)2/Vbxw/Subsystem1'
//  '<S52>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace'
//  '<S53>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Positive Trace'
//  '<S54>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM'
//  '<S55>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/trace(DCM)'
//  '<S56>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)'
//  '<S57>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)'
//  '<S58>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)'
//  '<S59>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/diag(DCM)'
//  '<S60>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) -sin(theta)'
//  '<S61>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(theta)sin(phi) - (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
//  '<S62>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(theta)sin(psi) + (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
//  '<S63>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/if s~=0; s=0.5//s'
//  '<S64>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/u(1) -(u(5)+u(9)) +1'
//  '<S65>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) +sin(theta)'
//  '<S66>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(theta)sin(phi) + (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
//  '<S67>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(theta)sin(psi) + (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
//  '<S68>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/if s~=0; s=0.5//s'
//  '<S69>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/u(5) -(u(1)+u(9)) +1'
//  '<S70>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) -sin(theta)'
//  '<S71>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(theta)sin(phi) + (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
//  '<S72>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(theta)sin(psi) - (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
//  '<S73>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/if s~=0; s=0.5//s'
//  '<S74>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/u(9) -(u(1)+u(5)) +1'
//  '<S75>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) +sin(theta)'
//  '<S76>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(theta)sin(phi) - (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
//  '<S77>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(theta)sin(psi) - (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
//  '<S78>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error'
//  '<S79>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal'
//  '<S80>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper'
//  '<S81>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotOrthogonal'
//  '<S82>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper'
//  '<S83>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal/Error'
//  '<S84>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal/Warning'
//  '<S85>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper/Error'
//  '<S86>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper/Warning'
//  '<S87>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotOrthogonal/transpose*dcm ~= eye(3)'
//  '<S88>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper/Determinant of 3x3 Matrix'
//  '<S89>'  : 'multiModeQuad_ROS/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper/determinant does not equal 1'
//  '<S90>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Anti-windup'
//  '<S91>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/D Gain'
//  '<S92>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Filter'
//  '<S93>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Filter ICs'
//  '<S94>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/I Gain'
//  '<S95>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Ideal P Gain'
//  '<S96>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Ideal P Gain Fdbk'
//  '<S97>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Integrator'
//  '<S98>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Integrator ICs'
//  '<S99>'  : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/N Copy'
//  '<S100>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/N Gain'
//  '<S101>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/P Copy'
//  '<S102>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Parallel P Gain'
//  '<S103>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Reset Signal'
//  '<S104>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Saturation'
//  '<S105>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Saturation Fdbk'
//  '<S106>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Sum'
//  '<S107>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Sum Fdbk'
//  '<S108>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Tracking Mode'
//  '<S109>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Tracking Mode Sum'
//  '<S110>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Tsamp - Integral'
//  '<S111>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Tsamp - Ngain'
//  '<S112>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/postSat Signal'
//  '<S113>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/preSat Signal'
//  '<S114>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Anti-windup/Passthrough'
//  '<S115>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/D Gain/Internal Parameters'
//  '<S116>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Filter/Cont. Filter'
//  '<S117>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Filter ICs/Internal IC - Filter'
//  '<S118>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/I Gain/Internal Parameters'
//  '<S119>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Ideal P Gain/Passthrough'
//  '<S120>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Ideal P Gain Fdbk/Disabled'
//  '<S121>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Integrator/Continuous'
//  '<S122>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Integrator ICs/Internal IC'
//  '<S123>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/N Copy/Disabled'
//  '<S124>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/N Gain/Internal Parameters'
//  '<S125>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/P Copy/Disabled'
//  '<S126>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Parallel P Gain/Internal Parameters'
//  '<S127>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Reset Signal/Disabled'
//  '<S128>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Saturation/Passthrough'
//  '<S129>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Saturation Fdbk/Disabled'
//  '<S130>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Sum/Sum_PID'
//  '<S131>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Sum Fdbk/Disabled'
//  '<S132>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Tracking Mode/Disabled'
//  '<S133>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Tracking Mode Sum/Passthrough'
//  '<S134>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Tsamp - Integral/Passthrough'
//  '<S135>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/Tsamp - Ngain/Passthrough'
//  '<S136>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/postSat Signal/Forward_Path'
//  '<S137>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity x/preSat Signal/Forward_Path'
//  '<S138>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Anti-windup'
//  '<S139>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/D Gain'
//  '<S140>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Filter'
//  '<S141>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Filter ICs'
//  '<S142>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/I Gain'
//  '<S143>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Ideal P Gain'
//  '<S144>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Ideal P Gain Fdbk'
//  '<S145>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Integrator'
//  '<S146>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Integrator ICs'
//  '<S147>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/N Copy'
//  '<S148>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/N Gain'
//  '<S149>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/P Copy'
//  '<S150>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Parallel P Gain'
//  '<S151>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Reset Signal'
//  '<S152>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Saturation'
//  '<S153>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Saturation Fdbk'
//  '<S154>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Sum'
//  '<S155>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Sum Fdbk'
//  '<S156>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Tracking Mode'
//  '<S157>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Tracking Mode Sum'
//  '<S158>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Tsamp - Integral'
//  '<S159>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Tsamp - Ngain'
//  '<S160>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/postSat Signal'
//  '<S161>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/preSat Signal'
//  '<S162>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Anti-windup/Passthrough'
//  '<S163>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/D Gain/Internal Parameters'
//  '<S164>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Filter/Cont. Filter'
//  '<S165>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Filter ICs/Internal IC - Filter'
//  '<S166>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/I Gain/Internal Parameters'
//  '<S167>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Ideal P Gain/Passthrough'
//  '<S168>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Ideal P Gain Fdbk/Disabled'
//  '<S169>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Integrator/Continuous'
//  '<S170>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Integrator ICs/Internal IC'
//  '<S171>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/N Copy/Disabled'
//  '<S172>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/N Gain/Internal Parameters'
//  '<S173>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/P Copy/Disabled'
//  '<S174>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Parallel P Gain/Internal Parameters'
//  '<S175>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Reset Signal/Disabled'
//  '<S176>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Saturation/Passthrough'
//  '<S177>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Saturation Fdbk/Disabled'
//  '<S178>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Sum/Sum_PID'
//  '<S179>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Sum Fdbk/Disabled'
//  '<S180>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Tracking Mode/Disabled'
//  '<S181>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Tracking Mode Sum/Passthrough'
//  '<S182>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Tsamp - Integral/Passthrough'
//  '<S183>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/Tsamp - Ngain/Passthrough'
//  '<S184>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/postSat Signal/Forward_Path'
//  '<S185>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity y/preSat Signal/Forward_Path'
//  '<S186>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Anti-windup'
//  '<S187>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/D Gain'
//  '<S188>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Filter'
//  '<S189>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Filter ICs'
//  '<S190>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/I Gain'
//  '<S191>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Ideal P Gain'
//  '<S192>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Ideal P Gain Fdbk'
//  '<S193>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Integrator'
//  '<S194>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Integrator ICs'
//  '<S195>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/N Copy'
//  '<S196>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/N Gain'
//  '<S197>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/P Copy'
//  '<S198>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Parallel P Gain'
//  '<S199>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Reset Signal'
//  '<S200>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Saturation'
//  '<S201>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Saturation Fdbk'
//  '<S202>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Sum'
//  '<S203>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Sum Fdbk'
//  '<S204>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Tracking Mode'
//  '<S205>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Tracking Mode Sum'
//  '<S206>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Tsamp - Integral'
//  '<S207>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Tsamp - Ngain'
//  '<S208>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/postSat Signal'
//  '<S209>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/preSat Signal'
//  '<S210>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Anti-windup/Passthrough'
//  '<S211>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/D Gain/Internal Parameters'
//  '<S212>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Filter/Cont. Filter'
//  '<S213>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Filter ICs/Internal IC - Filter'
//  '<S214>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/I Gain/Internal Parameters'
//  '<S215>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Ideal P Gain/Passthrough'
//  '<S216>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Ideal P Gain Fdbk/Disabled'
//  '<S217>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Integrator/Continuous'
//  '<S218>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Integrator ICs/Internal IC'
//  '<S219>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/N Copy/Disabled'
//  '<S220>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/N Gain/Internal Parameters'
//  '<S221>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/P Copy/Disabled'
//  '<S222>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Parallel P Gain/Internal Parameters'
//  '<S223>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Reset Signal/Disabled'
//  '<S224>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Saturation/Passthrough'
//  '<S225>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Saturation Fdbk/Disabled'
//  '<S226>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Sum/Sum_PID'
//  '<S227>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Sum Fdbk/Disabled'
//  '<S228>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Tracking Mode/Disabled'
//  '<S229>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Tracking Mode Sum/Passthrough'
//  '<S230>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Tsamp - Integral/Passthrough'
//  '<S231>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/Tsamp - Ngain/Passthrough'
//  '<S232>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/postSat Signal/Forward_Path'
//  '<S233>' : 'multiModeQuad_ROS/Subsystem1/PID Velocity z/preSat Signal/Forward_Path'
//  '<S234>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Anti-windup'
//  '<S235>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/D Gain'
//  '<S236>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Filter'
//  '<S237>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Filter ICs'
//  '<S238>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/I Gain'
//  '<S239>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Ideal P Gain'
//  '<S240>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Ideal P Gain Fdbk'
//  '<S241>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Integrator'
//  '<S242>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Integrator ICs'
//  '<S243>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/N Copy'
//  '<S244>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/N Gain'
//  '<S245>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/P Copy'
//  '<S246>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Parallel P Gain'
//  '<S247>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Reset Signal'
//  '<S248>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Saturation'
//  '<S249>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Saturation Fdbk'
//  '<S250>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Sum'
//  '<S251>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Sum Fdbk'
//  '<S252>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Tracking Mode'
//  '<S253>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Tracking Mode Sum'
//  '<S254>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Tsamp - Integral'
//  '<S255>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Tsamp - Ngain'
//  '<S256>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/postSat Signal'
//  '<S257>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/preSat Signal'
//  '<S258>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Anti-windup/Passthrough'
//  '<S259>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/D Gain/Internal Parameters'
//  '<S260>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Filter/Disc. Forward Euler Filter'
//  '<S261>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Filter ICs/Internal IC - Filter'
//  '<S262>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/I Gain/Internal Parameters'
//  '<S263>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Ideal P Gain/Internal Parameters'
//  '<S264>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Ideal P Gain Fdbk/Disabled'
//  '<S265>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Integrator/Discrete'
//  '<S266>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Integrator ICs/Internal IC'
//  '<S267>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/N Copy/Disabled'
//  '<S268>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/N Gain/Internal Parameters'
//  '<S269>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/P Copy/Disabled'
//  '<S270>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Parallel P Gain/Passthrough'
//  '<S271>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Reset Signal/Disabled'
//  '<S272>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Saturation/Passthrough'
//  '<S273>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Saturation Fdbk/Disabled'
//  '<S274>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Sum/Sum_PID'
//  '<S275>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Sum Fdbk/Disabled'
//  '<S276>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Tracking Mode/Disabled'
//  '<S277>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Tracking Mode Sum/Passthrough'
//  '<S278>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Tsamp - Integral/Passthrough'
//  '<S279>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/Tsamp - Ngain/Passthrough'
//  '<S280>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/postSat Signal/Forward_Path'
//  '<S281>' : 'multiModeQuad_ROS/Subsystem1/PID angula yaw/preSat Signal/Forward_Path'
//  '<S282>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Anti-windup'
//  '<S283>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/D Gain'
//  '<S284>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Filter'
//  '<S285>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Filter ICs'
//  '<S286>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/I Gain'
//  '<S287>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Ideal P Gain'
//  '<S288>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Ideal P Gain Fdbk'
//  '<S289>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Integrator'
//  '<S290>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Integrator ICs'
//  '<S291>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/N Copy'
//  '<S292>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/N Gain'
//  '<S293>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/P Copy'
//  '<S294>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Parallel P Gain'
//  '<S295>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Reset Signal'
//  '<S296>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Saturation'
//  '<S297>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Saturation Fdbk'
//  '<S298>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Sum'
//  '<S299>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Sum Fdbk'
//  '<S300>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Tracking Mode'
//  '<S301>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Tracking Mode Sum'
//  '<S302>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Tsamp - Integral'
//  '<S303>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Tsamp - Ngain'
//  '<S304>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/postSat Signal'
//  '<S305>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/preSat Signal'
//  '<S306>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Anti-windup/Passthrough'
//  '<S307>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/D Gain/Internal Parameters'
//  '<S308>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Filter/Disc. Forward Euler Filter'
//  '<S309>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Filter ICs/Internal IC - Filter'
//  '<S310>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/I Gain/Internal Parameters'
//  '<S311>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Ideal P Gain/Internal Parameters'
//  '<S312>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Ideal P Gain Fdbk/Disabled'
//  '<S313>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Integrator/Discrete'
//  '<S314>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Integrator ICs/Internal IC'
//  '<S315>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/N Copy/Disabled'
//  '<S316>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/N Gain/Internal Parameters'
//  '<S317>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/P Copy/Disabled'
//  '<S318>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Parallel P Gain/Passthrough'
//  '<S319>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Reset Signal/Disabled'
//  '<S320>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Saturation/Passthrough'
//  '<S321>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Saturation Fdbk/Disabled'
//  '<S322>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Sum/Sum_PID'
//  '<S323>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Sum Fdbk/Disabled'
//  '<S324>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Tracking Mode/Disabled'
//  '<S325>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Tracking Mode Sum/Passthrough'
//  '<S326>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Tsamp - Integral/Passthrough'
//  '<S327>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/Tsamp - Ngain/Passthrough'
//  '<S328>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/postSat Signal/Forward_Path'
//  '<S329>' : 'multiModeQuad_ROS/Subsystem1/PID angular pitch/preSat Signal/Forward_Path'
//  '<S330>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Anti-windup'
//  '<S331>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/D Gain'
//  '<S332>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Filter'
//  '<S333>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Filter ICs'
//  '<S334>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/I Gain'
//  '<S335>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Ideal P Gain'
//  '<S336>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Ideal P Gain Fdbk'
//  '<S337>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Integrator'
//  '<S338>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Integrator ICs'
//  '<S339>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/N Copy'
//  '<S340>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/N Gain'
//  '<S341>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/P Copy'
//  '<S342>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Parallel P Gain'
//  '<S343>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Reset Signal'
//  '<S344>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Saturation'
//  '<S345>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Saturation Fdbk'
//  '<S346>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Sum'
//  '<S347>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Sum Fdbk'
//  '<S348>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Tracking Mode'
//  '<S349>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Tracking Mode Sum'
//  '<S350>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Tsamp - Integral'
//  '<S351>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Tsamp - Ngain'
//  '<S352>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/postSat Signal'
//  '<S353>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/preSat Signal'
//  '<S354>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Anti-windup/Passthrough'
//  '<S355>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/D Gain/Internal Parameters'
//  '<S356>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Filter/Disc. Forward Euler Filter'
//  '<S357>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Filter ICs/Internal IC - Filter'
//  '<S358>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/I Gain/Internal Parameters'
//  '<S359>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Ideal P Gain/Internal Parameters'
//  '<S360>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Ideal P Gain Fdbk/Disabled'
//  '<S361>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Integrator/Discrete'
//  '<S362>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Integrator ICs/Internal IC'
//  '<S363>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/N Copy/Disabled'
//  '<S364>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/N Gain/Internal Parameters'
//  '<S365>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/P Copy/Disabled'
//  '<S366>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Parallel P Gain/Passthrough'
//  '<S367>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Reset Signal/Disabled'
//  '<S368>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Saturation/Passthrough'
//  '<S369>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Saturation Fdbk/Disabled'
//  '<S370>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Sum/Sum_PID'
//  '<S371>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Sum Fdbk/Disabled'
//  '<S372>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Tracking Mode/Disabled'
//  '<S373>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Tracking Mode Sum/Passthrough'
//  '<S374>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Tsamp - Integral/Passthrough'
//  '<S375>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/Tsamp - Ngain/Passthrough'
//  '<S376>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/postSat Signal/Forward_Path'
//  '<S377>' : 'multiModeQuad_ROS/Subsystem1/PID angular roll/preSat Signal/Forward_Path'

#endif                                 // RTW_HEADER_multiModeQuad_ROS_h_

//
// File trailer for generated code.
//
// [EOF]
//
