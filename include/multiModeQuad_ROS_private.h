//
// File: multiModeQuad_ROS_private.h
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
#ifndef RTW_HEADER_multiModeQuad_ROS_private_h_
#define RTW_HEADER_multiModeQuad_ROS_private_h_
#include "rtwtypes.h"
#include "multiModeQuad_ROS.h"

// Private macros used by the generated code to access rtModel
#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

extern int32_T div_nzp_s32(int32_T numerator, int32_T denominator);
extern void multiModeQuad_ROS_timetosecnsec(real_T rtu_time, real_T *rty_sec,
  real_T *rty_nsec);

// private model entry point functions
extern void multiModeQuad_ROS_derivatives(void);

#endif                               // RTW_HEADER_multiModeQuad_ROS_private_h_

//
// File trailer for generated code.
//
// [EOF]
//
