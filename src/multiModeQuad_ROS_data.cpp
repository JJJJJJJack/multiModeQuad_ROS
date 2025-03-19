//
// File: multiModeQuad_ROS_data.cpp
//
// Code generated for Simulink model 'multiModeQuad_ROS'.
//
// Model version                  : 1.71
// Simulink Coder version         : 9.6 (R2021b) 14-May-2021
// C/C++ source code generated on : Thu Mar 20 00:46:27 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: Generic->Unspecified (assume 32-bit Generic)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#include "multiModeQuad_ROS.h"
#include "multiModeQuad_ROS_private.h"

// Block parameters (default storage)
P_multiModeQuad_ROS_T multiModeQuad_ROS_P = {
  // Mask Parameter: PIDVelocityx_D
  //  Referenced by: '<S116>/Derivative Gain'

  0.01,

  // Mask Parameter: PIDVelocityy_D
  //  Referenced by: '<S164>/Derivative Gain'

  0.01,

  // Mask Parameter: PIDangularroll_D
  //  Referenced by: '<S356>/Derivative Gain'

  0.05,

  // Mask Parameter: PIDangularpitch_D
  //  Referenced by: '<S308>/Derivative Gain'

  0.05,

  // Mask Parameter: PIDangulayaw_D
  //  Referenced by: '<S260>/Derivative Gain'

  0.04,

  // Mask Parameter: PIDVelocityz_D
  //  Referenced by: '<S212>/Derivative Gain'

  2.0,

  // Mask Parameter: PIDangulayaw_I
  //  Referenced by: '<S263>/Integral Gain'

  0.05,

  // Mask Parameter: PIDangularpitch_I
  //  Referenced by: '<S311>/Integral Gain'

  0.05,

  // Mask Parameter: PIDangularroll_I
  //  Referenced by: '<S359>/Integral Gain'

  0.05,

  // Mask Parameter: PIDVelocityx_I
  //  Referenced by: '<S119>/Integral Gain'

  0.0,

  // Mask Parameter: PIDVelocityy_I
  //  Referenced by: '<S167>/Integral Gain'

  0.02,

  // Mask Parameter: PIDVelocityz_I
  //  Referenced by: '<S215>/Integral Gain'

  0.2,

  // Mask Parameter: PIDVelocityx_InitialConditionFo
  //  Referenced by: '<S117>/Filter'

  0.0,

  // Mask Parameter: PIDVelocityy_InitialConditionFo
  //  Referenced by: '<S165>/Filter'

  0.0,

  // Mask Parameter: PIDangularroll_InitialCondition
  //  Referenced by: '<S357>/Filter'

  0.0,

  // Mask Parameter: PIDangularpitch_InitialConditio
  //  Referenced by: '<S309>/Filter'

  0.0,

  // Mask Parameter: PIDangulayaw_InitialConditionFo
  //  Referenced by: '<S261>/Filter'

  0.0,

  // Mask Parameter: PIDVelocityz_InitialConditionFo
  //  Referenced by: '<S213>/Filter'

  0.0,

  // Mask Parameter: PIDVelocityx_InitialCondition_o
  //  Referenced by: '<S122>/Integrator'

  0.0,

  // Mask Parameter: PIDVelocityy_InitialCondition_i
  //  Referenced by: '<S170>/Integrator'

  0.0,

  // Mask Parameter: PIDangularroll_InitialConditi_o
  //  Referenced by: '<S362>/Integrator'

  0.0,

  // Mask Parameter: PIDangularpitch_InitialCondit_b
  //  Referenced by: '<S314>/Integrator'

  0.0,

  // Mask Parameter: PIDangulayaw_InitialCondition_m
  //  Referenced by: '<S266>/Integrator'

  0.0,

  // Mask Parameter: PIDVelocityz_InitialCondition_g
  //  Referenced by: '<S218>/Integrator'

  0.0,

  // Mask Parameter: PIDVelocityx_N
  //  Referenced by: '<S125>/Filter Coefficient'

  100.0,

  // Mask Parameter: PIDVelocityy_N
  //  Referenced by: '<S173>/Filter Coefficient'

  100.0,

  // Mask Parameter: PIDangularroll_N
  //  Referenced by: '<S365>/Filter Coefficient'

  300.0,

  // Mask Parameter: PIDangularpitch_N
  //  Referenced by: '<S317>/Filter Coefficient'

  300.0,

  // Mask Parameter: PIDangulayaw_N
  //  Referenced by: '<S269>/Filter Coefficient'

  300.0,

  // Mask Parameter: PIDVelocityz_N
  //  Referenced by: '<S221>/Filter Coefficient'

  100.0,

  // Mask Parameter: PIDVelocityx_P
  //  Referenced by: '<S127>/Proportional Gain'

  3.0,

  // Mask Parameter: PIDVelocityy_P
  //  Referenced by: '<S175>/Proportional Gain'

  3.0,

  // Mask Parameter: PIDangularroll_P
  //  Referenced by: '<S360>/Proportional Gain'

  0.6,

  // Mask Parameter: PIDangularpitch_P
  //  Referenced by: '<S312>/Proportional Gain'

  0.6,

  // Mask Parameter: PIDangulayaw_P
  //  Referenced by: '<S264>/Proportional Gain'

  1.0,

  // Mask Parameter: PIDVelocityz_P
  //  Referenced by: '<S223>/Proportional Gain'

  3.0,

  // Mask Parameter: uDOFEulerAngles2_Vm_0
  //  Referenced by: '<S18>/ub,vb,wb'

  { 0.0, 0.0, 0.0 },

  // Mask Parameter: DirectionCosineMatrixtoQuaterni
  //  Referenced by:
  //    '<S55>/Constant'
  //    '<S80>/Constant'
  //    '<S81>/Constant'

  1.0,

  // Mask Parameter: uDOFEulerAngles2_eul_0
  //  Referenced by: '<S35>/phi theta psi'

  { 0.0, 0.0, 0.0 },

  // Mask Parameter: uDOFEulerAngles2_inertia
  //  Referenced by: '<S37>/Constant1'

  { 0.015, 0.0, 0.0, 0.0, 0.015, 0.0, 0.0, 0.0, 0.03 },

  // Mask Parameter: uDOFEulerAngles2_mass_0
  //  Referenced by: '<S37>/Constant'

  2.0,

  // Mask Parameter: uDOFEulerAngles2_pm_0
  //  Referenced by: '<S18>/p,q,r '

  { 0.0, 0.0, 0.0 },

  // Mask Parameter: DirectionCosineMatrixtoQuater_h
  //  Referenced by:
  //    '<S88>/Constant'
  //    '<S90>/Constant'

  4.4408920985006262E-16,

  // Mask Parameter: uDOFEulerAngles2_xme_0
  //  Referenced by: '<S18>/xe,ye,ze'

  { 0.0, 0.0, 0.0 },

  // Computed Parameter: Constant_Value
  //  Referenced by: '<S2>/Constant'

  {
    {
      0U,                              // Seq

      {
        0.0,                           // Sec
        0.0                            // Nsec
      },                               // Stamp

      {
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U }
      ,                                // FrameId

      {
        0U,                            // CurrentLength
        0U                             // ReceivedLength
      }                                // FrameId_SL_Info
    },                                 // Header

    {
      0.0,                             // X
      0.0,                             // Y
      0.0,                             // Z
      0.0                              // W
    },                                 // Orientation

    {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
    ,                                  // OrientationCovariance

    {
      0.0,                             // X
      0.0,                             // Y
      0.0                              // Z
    },                                 // AngularVelocity

    {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
    ,                                  // AngularVelocityCovariance

    {
      0.0,                             // X
      0.0,                             // Y
      0.0                              // Z
    },                                 // LinearAcceleration

    {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
    // LinearAccelerationCovariance
  },

  // Computed Parameter: Constant_Value_f
  //  Referenced by: '<S1>/Constant'

  {
    {
      0U,                              // Seq

      {
        0.0,                           // Sec
        0.0                            // Nsec
      },                               // Stamp

      {
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U }
      ,                                // FrameId

      {
        0U,                            // CurrentLength
        0U                             // ReceivedLength
      }                                // FrameId_SL_Info
    },                                 // Header

    {
      {
        0.0,                           // X
        0.0,                           // Y
        0.0                            // Z
      },                               // Position

      {
        0.0,                           // X
        0.0,                           // Y
        0.0,                           // Z
        0.0                            // W
      }                                // Orientation
    }                                  // Pose
  },

  // Computed Parameter: Out1_Y0
  //  Referenced by: '<S14>/Out1'

  {
    {
      0U,                              // Seq

      {
        0.0,                           // Sec
        0.0                            // Nsec
      },                               // Stamp

      {
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U }
      ,                                // FrameId

      {
        0U,                            // CurrentLength
        0U                             // ReceivedLength
      }                                // FrameId_SL_Info
    },                                 // Header

    {
      {
        0.0,                           // X
        0.0,                           // Y
        0.0                            // Z
      },                               // Position

      {
        0.0,                           // X
        0.0,                           // Y
        0.0,                           // Z
        0.0                            // W
      }                                // Orientation
    }                                  // Pose
  },

  // Computed Parameter: Constant_Value_a
  //  Referenced by: '<S6>/Constant'

  {
    {
      0U,                              // Seq

      {
        0.0,                           // Sec
        0.0                            // Nsec
      },                               // Stamp

      {
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U }
      ,                                // FrameId

      {
        0U,                            // CurrentLength
        0U                             // ReceivedLength
      }                                // FrameId_SL_Info
    },                                 // Header

    {
      {
        0.0,                           // X
        0.0,                           // Y
        0.0                            // Z
      },                               // Position

      {
        0.0,                           // X
        0.0,                           // Y
        0.0,                           // Z
        0.0                            // W
      }                                // Orientation
    }                                  // Pose
  },

  // Computed Parameter: Out1_Y0_k
  //  Referenced by: '<S15>/Out1'

  {
    {
      0U,                              // Seq

      {
        0.0,                           // Sec
        0.0                            // Nsec
      },                               // Stamp

      {
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U }
      ,                                // FrameId

      {
        0U,                            // CurrentLength
        0U                             // ReceivedLength
      }                                // FrameId_SL_Info
    },                                 // Header

    {
      {
        0.0,                           // X
        0.0,                           // Y
        0.0                            // Z
      },                               // Linear

      {
        0.0,                           // X
        0.0,                           // Y
        0.0                            // Z
      }                                // Angular
    }                                  // Twist
  },

  // Computed Parameter: Constant_Value_c
  //  Referenced by: '<S7>/Constant'

  {
    {
      0U,                              // Seq

      {
        0.0,                           // Sec
        0.0                            // Nsec
      },                               // Stamp

      {
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U }
      ,                                // FrameId

      {
        0U,                            // CurrentLength
        0U                             // ReceivedLength
      }                                // FrameId_SL_Info
    },                                 // Header

    {
      {
        0.0,                           // X
        0.0,                           // Y
        0.0                            // Z
      },                               // Linear

      {
        0.0,                           // X
        0.0,                           // Y
        0.0                            // Z
      }                                // Angular
    }                                  // Twist
  },

  // Computed Parameter: Out1_Y0_g
  //  Referenced by: '<S17>/Out1'

  {
    {
      0.0,                             // X
      0.0,                             // Y
      0.0                              // Z
    },                                 // Linear

    {
      0.0,                             // X
      0.0,                             // Y
      0.0                              // Z
    }                                  // Angular
  },

  // Computed Parameter: Constant_Value_fz
  //  Referenced by: '<S9>/Constant'

  {
    {
      0.0,                             // X
      0.0,                             // Y
      0.0                              // Z
    },                                 // Linear

    {
      0.0,                             // X
      0.0,                             // Y
      0.0                              // Z
    }                                  // Angular
  },

  // Computed Parameter: Out1_Y0_kv
  //  Referenced by: '<S16>/Out1'

  {
    0.0F                               // Data
  },

  // Computed Parameter: Constant_Value_p
  //  Referenced by: '<S8>/Constant'

  {
    0.0F                               // Data
  },

  // Computed Parameter: Out1_Y0_p
  //  Referenced by: '<S13>/Out1'

  {
    0                                  // Data
  },

  // Computed Parameter: Constant_Value_e
  //  Referenced by: '<S3>/Constant'

  {
    0                                  // Data
  },

  // Expression: 1
  //  Referenced by: '<S54>/Constant'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S54>/Gain'

  0.5,

  // Expression: 2
  //  Referenced by: '<S54>/Gain1'

  2.0,

  // Expression: 1
  //  Referenced by: '<S70>/Constant'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S69>/Constant1'

  0.5,

  // Expression: [0 1]
  //  Referenced by: '<S69>/Constant2'

  { 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S58>/Gain1'

  1.0,

  // Expression: 1
  //  Referenced by: '<S58>/Gain3'

  1.0,

  // Expression: 1
  //  Referenced by: '<S58>/Gain4'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S58>/Gain'

  0.5,

  // Expression: 1
  //  Referenced by: '<S75>/Constant'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S74>/Constant1'

  0.5,

  // Expression: [0 1]
  //  Referenced by: '<S74>/Constant2'

  { 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S59>/Gain1'

  1.0,

  // Expression: 1
  //  Referenced by: '<S59>/Gain2'

  1.0,

  // Expression: 1
  //  Referenced by: '<S59>/Gain3'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S59>/Gain'

  0.5,

  // Expression: 1
  //  Referenced by: '<S65>/Constant'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S64>/Constant1'

  0.5,

  // Expression: [0 1]
  //  Referenced by: '<S64>/Constant2'

  { 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S57>/Gain1'

  1.0,

  // Expression: 1
  //  Referenced by: '<S57>/Gain2'

  1.0,

  // Expression: 1
  //  Referenced by: '<S57>/Gain3'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S57>/Gain'

  0.5,

  // Expression: 0
  //  Referenced by: '<S81>/Constant1'

  0.0,

  // Expression: 0
  //  Referenced by: '<S80>/Constant1'

  0.0,

  // Expression: -eye(3)
  //  Referenced by: '<S82>/Bias1'

  { -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0 },

  // Expression: -1
  //  Referenced by: '<S83>/Bias'

  -1.0,

  // Expression: 0.5
  //  Referenced by: '<S32>/1//2'

  0.5,

  // Expression: zeros(3)
  //  Referenced by: '<S37>/Constant2'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pi
  //  Referenced by: '<S35>/phi theta psi'

  3.1415926535897931,

  // Expression: -pi
  //  Referenced by: '<S35>/phi theta psi'

  -3.1415926535897931,

  // Expression: -1
  //  Referenced by: '<S33>/Gain'

  -1.0,

  // Computed Parameter: TransferFcn2_A
  //  Referenced by: '<S10>/Transfer Fcn2'

  -24.390243902439025,

  // Computed Parameter: TransferFcn2_C
  //  Referenced by: '<S10>/Transfer Fcn2'

  24.390243902439025,

  // Expression: 10
  //  Referenced by: '<S10>/Saturation'

  10.0,

  // Expression: 0
  //  Referenced by: '<S10>/Saturation'

  0.0,

  // Computed Parameter: TransferFcn1_A
  //  Referenced by: '<S10>/Transfer Fcn1'

  -24.390243902439025,

  // Computed Parameter: TransferFcn1_C
  //  Referenced by: '<S10>/Transfer Fcn1'

  24.390243902439025,

  // Expression: 10
  //  Referenced by: '<S10>/Saturation1'

  10.0,

  // Expression: 0
  //  Referenced by: '<S10>/Saturation1'

  0.0,

  // Computed Parameter: TransferFcn3_A
  //  Referenced by: '<S10>/Transfer Fcn3'

  -24.390243902439025,

  // Computed Parameter: TransferFcn3_C
  //  Referenced by: '<S10>/Transfer Fcn3'

  24.390243902439025,

  // Expression: 10
  //  Referenced by: '<S10>/Saturation2'

  10.0,

  // Expression: 0
  //  Referenced by: '<S10>/Saturation2'

  0.0,

  // Computed Parameter: TransferFcn4_A
  //  Referenced by: '<S10>/Transfer Fcn4'

  -24.390243902439025,

  // Computed Parameter: TransferFcn4_C
  //  Referenced by: '<S10>/Transfer Fcn4'

  24.390243902439025,

  // Expression: 10
  //  Referenced by: '<S10>/Saturation3'

  10.0,

  // Expression: 0
  //  Referenced by: '<S10>/Saturation3'

  0.0,

  // Expression: 0
  //  Referenced by: '<S21>/Constant4'

  0.0,

  // Expression: 0
  //  Referenced by: '<S21>/Constant5'

  0.0,

  // Expression: 2
  //  Referenced by: '<S21>/Mass'

  2.0,

  // Expression: 9.8
  //  Referenced by: '<S21>/g'

  9.8,

  // Expression: -1
  //  Referenced by: '<S34>/Gain'

  -1.0,

  // Expression: [1 0 0 0]
  //  Referenced by: '<S20>/Merge'

  { 1.0, 0.0, 0.0, 0.0 },

  // Computed Parameter: Filter_gainval
  //  Referenced by: '<S357>/Filter'

  0.005,

  // Computed Parameter: Integrator_gainval
  //  Referenced by: '<S362>/Integrator'

  0.005,

  // Computed Parameter: Integrator_gainval_e
  //  Referenced by: '<S314>/Integrator'

  0.005,

  // Computed Parameter: Filter_gainval_j
  //  Referenced by: '<S309>/Filter'

  0.005,

  // Computed Parameter: Integrator_gainval_a
  //  Referenced by: '<S266>/Integrator'

  0.005,

  // Computed Parameter: Filter_gainval_c
  //  Referenced by: '<S261>/Filter'

  0.005,

  // Expression: 3.13
  //  Referenced by: '<S10>/Hover throttle'

  3.13,

  // Expression: inf
  //  Referenced by: '<S10>/Saturation4'

  0.0,

  // Expression: 0
  //  Referenced by: '<S10>/Saturation4'

  0.0,

  // Expression: 0.5
  //  Referenced by: '<S10>/Komega'

  0.5,

  // Expression: inf
  //  Referenced by: '<S10>/Saturation5'

  0.0,

  // Expression: 0
  //  Referenced by: '<S10>/Saturation5'

  0.0,

  // Expression: 0.5
  //  Referenced by: '<S10>/Komega1'

  0.5,

  // Expression: inf
  //  Referenced by: '<S10>/Saturation6'

  0.0,

  // Expression: 0
  //  Referenced by: '<S10>/Saturation6'

  0.0,

  // Expression: 0.5
  //  Referenced by: '<S10>/Komega2'

  0.5,

  // Expression: inf
  //  Referenced by: '<S10>/Saturation7'

  0.0,

  // Expression: 0
  //  Referenced by: '<S10>/Saturation7'

  0.0,

  // Expression: 0.5
  //  Referenced by: '<S10>/Komega3'

  0.5,

  // Computed Parameter: Gain_Gain_dx
  //  Referenced by: '<S10>/Gain'

  10.0F,

  // Computed Parameter: Switch3_Threshold
  //  Referenced by: '<S10>/Switch3'

  1,

  // Computed Parameter: Switch_Threshold
  //  Referenced by: '<S10>/Switch'

  1,

  // Computed Parameter: Switch1_Threshold
  //  Referenced by: '<S10>/Switch1'

  2,

  // Computed Parameter: Switch2_Threshold
  //  Referenced by: '<S10>/Switch2'

  1
};

//
// File trailer for generated code.
//
// [EOF]
//
