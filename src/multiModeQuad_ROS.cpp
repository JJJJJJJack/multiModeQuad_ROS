//
// File: multiModeQuad_ROS.cpp
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

// Block signals (default storage)
B_multiModeQuad_ROS_T multiModeQuad_ROS_B;

// Continuous states
X_multiModeQuad_ROS_T multiModeQuad_ROS_X;

// Block states (default storage)
DW_multiModeQuad_ROS_T multiModeQuad_ROS_DW;

// Periodic continuous states
PeriodicIndX_multiModeQuad_RO_T multiModeQuad_ROS_PeriodicIndX;
PeriodicRngX_multiModeQuad_RO_T multiModeQuad_ROS_PeriodicRngX;

// Real-time model
RT_MODEL_multiModeQuad_ROS_T multiModeQuad_ROS_M_ = RT_MODEL_multiModeQuad_ROS_T
  ();
RT_MODEL_multiModeQuad_ROS_T *const multiModeQuad_ROS_M = &multiModeQuad_ROS_M_;

// Forward declaration for local functions
static boolean_T multiModeQuad_ROS_anyNonFinite(const real_T x[16]);
static real_T multiModeQuad_ROS_rt_hypotd_snf(real_T u0, real_T u1);
static void multiModeQuad_ROS_xzggbal(creal_T A[16], int32_T *ilo, int32_T *ihi,
  int32_T rscale[4]);
static void multiModeQuad_ROS_sqrt(creal_T *x);
static void multiModeQuad_ROS_xzlartg_h(const creal_T f, const creal_T g, real_T
  *cs, creal_T *sn);
static void multiModeQuad_ROS_xzlartg(const creal_T f, const creal_T g, real_T
  *cs, creal_T *sn, creal_T *r);
static void multiModeQuad_ROS_xzhgeqz(creal_T A[16], int32_T ilo, int32_T ihi,
  creal_T Z[16], int32_T *info, creal_T alpha1[4], creal_T beta1[4]);
static void multiModeQuad_ROS_xztgevc(const creal_T A[16], creal_T V[16]);
static void multiModeQuad_ROS_eigStandard(const real_T A[16], creal_T V[16],
  creal_T D[4]);
static real_T multiModeQuad_ROS_xnrm2(int32_T n, const real_T x[16], int32_T ix0);
static void multiModeQuad_ROS_xzlarf(int32_T m, int32_T n, int32_T iv0, real_T
  tau, real_T C[16], int32_T ic0, real_T work[4]);
static real_T multiModeQuad_ROS_xnrm2_l(int32_T n, const real_T x[3]);
static real_T multiModeQuad_ROS_xzlarfg(int32_T n, real_T *alpha1, real_T x[3]);
static void multiModeQuad_ROS_xdlanv2(real_T *a, real_T *b, real_T *c, real_T *d,
  real_T *rt1r, real_T *rt1i, real_T *rt2r, real_T *rt2i, real_T *cs, real_T *sn);
static void multiModeQuad_ROS_xrot(int32_T n, real_T x[16], int32_T ix0, int32_T
  iy0, real_T c, real_T s);
static void multiModeQuad_ROS_xrot_i(real_T x[16], int32_T ix0, int32_T iy0,
  real_T c, real_T s);
static int32_T multiModeQuad_ROS_xhseqr(real_T h[16], real_T z[16]);
static void multiModeQuad_ROS_schur(const real_T A[16], real_T V[16], real_T T
  [16]);
static void multiModeQuad_ROS_quat2axang(real_T q[4], real_T axang[4]);
static void rt_mrdivide_U1d1x3_U2d_9vOrDY_i(const real_T u0[3], const real_T u1
  [9], real_T y[3]);
int32_T div_nzp_s32(int32_T numerator, int32_T denominator)
{
  uint32_T tempAbsQuotient;
  tempAbsQuotient = (numerator < 0 ? ~static_cast<uint32_T>(numerator) + 1U :
                     static_cast<uint32_T>(numerator)) / (denominator < 0 ? ~
    static_cast<uint32_T>(denominator) + 1U : static_cast<uint32_T>(denominator));
  return (numerator < 0) != (denominator < 0) ? -static_cast<int32_T>
    (tempAbsQuotient) : static_cast<int32_T>(tempAbsQuotient);
}

// State reduction function
void local_stateReduction(real_T* x, int_T* p, int_T n, real_T* r)
{
  int_T i, j;
  for (i = 0, j = 0; i < n; ++i, ++j) {
    int_T k = p[i];
    real_T lb = r[j++];
    real_T xk = x[k]-lb;
    real_T rk = r[j]-lb;
    int_T q = (int_T) floor(xk/rk);
    if (q) {
      x[k] = xk-q*rk+lb;
    }
  }
}

//
// This function updates continuous states using the ODE3 fixed-step
// solver algorithm
//
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  // Solver Matrices
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = static_cast<ODE3_IntgData *>(rtsiGetSolverData(si));
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 22;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  // Save the state values at time t in y, we'll use x as ynew.
  (void) memcpy(y, x,
                static_cast<uint_T>(nXc)*sizeof(real_T));

  // Assumes that rtsiSetT and ModelOutputs are up-to-date
  // f0 = f(t,y)
  rtsiSetdX(si, f0);
  multiModeQuad_ROS_derivatives();

  // f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*));
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  multiModeQuad_ROS_step();
  multiModeQuad_ROS_derivatives();

  // f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*));
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  multiModeQuad_ROS_step();
  multiModeQuad_ROS_derivatives();

  // tnew = t + hA(3);
  // ynew = y + f*hB(:,3);
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  local_stateReduction(rtsiGetContStates(si), rtsiGetPeriodicContStateIndices(si),
                       3,
                       rtsiGetPeriodicContStateRanges(si));
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

//
// Output and update for atomic system:
//    '<Root>/time to sec & nsec'
//    '<Root>/time to sec & nsec1'
//
void multiModeQuad_ROS_timetosecnsec(real_T rtu_time, real_T *rty_sec, real_T
  *rty_nsec)
{
  real_T sec;
  sec = floor(rtu_time);
  *rty_nsec = (rtu_time - sec) * 1.0E+9;
  *rty_sec = sec;
}

// Function for MATLAB Function: '<S10>/Attitude control'
static boolean_T multiModeQuad_ROS_anyNonFinite(const real_T x[16])
{
  boolean_T b_p;
  b_p = true;
  for (int32_T k = 0; k < 16; k++) {
    real_T x_0;
    x_0 = x[k];
    if (b_p && (rtIsInf(x_0) || rtIsNaN(x_0))) {
      b_p = false;
    }
  }

  return !b_p;
}

static real_T multiModeQuad_ROS_rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  multiModeQuad_ROS_B.a = fabs(u0);
  y = fabs(u1);
  if (multiModeQuad_ROS_B.a < y) {
    multiModeQuad_ROS_B.a /= y;
    y *= sqrt(multiModeQuad_ROS_B.a * multiModeQuad_ROS_B.a + 1.0);
  } else if (multiModeQuad_ROS_B.a > y) {
    y /= multiModeQuad_ROS_B.a;
    y = sqrt(y * y + 1.0) * multiModeQuad_ROS_B.a;
  } else if (!rtIsNaN(y)) {
    y = multiModeQuad_ROS_B.a * 1.4142135623730951;
  }

  return y;
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_xzggbal(creal_T A[16], int32_T *ilo, int32_T *ihi,
  int32_T rscale[4])
{
  real_T atmp_im;
  real_T atmp_re;
  int32_T atmp_re_tmp_tmp;
  int32_T exitg2;
  int32_T i;
  int32_T ii;
  int32_T j;
  int32_T jj;
  int32_T nzcount;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T found;
  rscale[0] = 1;
  rscale[1] = 1;
  rscale[2] = 1;
  rscale[3] = 1;
  *ilo = 1;
  *ihi = 4;
  do {
    exitg2 = 0;
    i = 0;
    j = 0;
    found = false;
    ii = *ihi;
    exitg3 = false;
    while ((!exitg3) && (ii > 0)) {
      nzcount = 0;
      i = ii;
      j = *ihi;
      jj = 0;
      exitg4 = false;
      while ((!exitg4) && (jj <= *ihi - 1)) {
        atmp_re_tmp_tmp = ((jj << 2) + ii) - 1;
        if ((A[atmp_re_tmp_tmp].re != 0.0) || (A[atmp_re_tmp_tmp].im != 0.0) ||
            (jj + 1 == ii)) {
          if (nzcount == 0) {
            j = jj + 1;
            nzcount = 1;
            jj++;
          } else {
            nzcount = 2;
            exitg4 = true;
          }
        } else {
          jj++;
        }
      }

      if (nzcount < 2) {
        found = true;
        exitg3 = true;
      } else {
        ii--;
      }
    }

    if (!found) {
      exitg2 = 2;
    } else {
      if (i != *ihi) {
        atmp_re = A[i - 1].re;
        atmp_im = A[i - 1].im;
        A[i - 1] = A[*ihi - 1];
        A[*ihi - 1].re = atmp_re;
        A[*ihi - 1].im = atmp_im;
        atmp_re = A[i + 3].re;
        atmp_im = A[i + 3].im;
        A[i + 3] = A[*ihi + 3];
        A[*ihi + 3].re = atmp_re;
        A[*ihi + 3].im = atmp_im;
        atmp_re = A[i + 7].re;
        atmp_im = A[i + 7].im;
        A[i + 7] = A[*ihi + 7];
        A[*ihi + 7].re = atmp_re;
        A[*ihi + 7].im = atmp_im;
        atmp_re = A[i + 11].re;
        atmp_im = A[i + 11].im;
        A[i + 11] = A[*ihi + 11];
        A[*ihi + 11].re = atmp_re;
        A[*ihi + 11].im = atmp_im;
      }

      if (j != *ihi) {
        for (ii = 0; ii < *ihi; ii++) {
          i = ((j - 1) << 2) + ii;
          atmp_re = A[i].re;
          atmp_im = A[i].im;
          atmp_re_tmp_tmp = ((*ihi - 1) << 2) + ii;
          A[i] = A[atmp_re_tmp_tmp];
          A[atmp_re_tmp_tmp].re = atmp_re;
          A[atmp_re_tmp_tmp].im = atmp_im;
        }
      }

      rscale[*ihi - 1] = j;
      (*ihi)--;
      if (*ihi == 1) {
        rscale[0] = 1;
        exitg2 = 1;
      }
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
    int32_T exitg1;
    do {
      exitg1 = 0;
      ii = 0;
      j = 0;
      found = false;
      i = *ilo;
      exitg3 = false;
      while ((!exitg3) && (i <= *ihi)) {
        nzcount = 0;
        ii = *ihi;
        j = i;
        jj = *ilo;
        exitg4 = false;
        while ((!exitg4) && (jj <= *ihi)) {
          atmp_re_tmp_tmp = (((i - 1) << 2) + jj) - 1;
          if ((A[atmp_re_tmp_tmp].re != 0.0) || (A[atmp_re_tmp_tmp].im != 0.0) ||
              (jj == i)) {
            if (nzcount == 0) {
              ii = jj;
              nzcount = 1;
              jj++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            jj++;
          }
        }

        if (nzcount < 2) {
          found = true;
          exitg3 = true;
        } else {
          i++;
        }
      }

      if (!found) {
        exitg1 = 1;
      } else {
        if (ii != *ilo) {
          for (nzcount = *ilo - 1; nzcount + 1 < 5; nzcount++) {
            atmp_re_tmp_tmp = nzcount << 2;
            i = (atmp_re_tmp_tmp + ii) - 1;
            atmp_re = A[i].re;
            atmp_im = A[i].im;
            atmp_re_tmp_tmp = (atmp_re_tmp_tmp + *ilo) - 1;
            A[i] = A[atmp_re_tmp_tmp];
            A[atmp_re_tmp_tmp].re = atmp_re;
            A[atmp_re_tmp_tmp].im = atmp_im;
          }
        }

        if (j != *ilo) {
          for (ii = 0; ii < *ihi; ii++) {
            i = ((j - 1) << 2) + ii;
            atmp_re = A[i].re;
            atmp_im = A[i].im;
            atmp_re_tmp_tmp = ((*ilo - 1) << 2) + ii;
            A[i] = A[atmp_re_tmp_tmp];
            A[atmp_re_tmp_tmp].re = atmp_re;
            A[atmp_re_tmp_tmp].im = atmp_im;
          }
        }

        rscale[*ilo - 1] = j;
        (*ilo)++;
        if (*ilo == *ihi) {
          rscale[*ilo - 1] = *ilo;
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_sqrt(creal_T *x)
{
  real_T absxi;
  if (x->im == 0.0) {
    if (x->re < 0.0) {
      multiModeQuad_ROS_B.absxr = 0.0;
      absxi = sqrt(-x->re);
    } else {
      multiModeQuad_ROS_B.absxr = sqrt(x->re);
      absxi = 0.0;
    }
  } else if (x->re == 0.0) {
    if (x->im < 0.0) {
      multiModeQuad_ROS_B.absxr = sqrt(-x->im / 2.0);
      absxi = -multiModeQuad_ROS_B.absxr;
    } else {
      multiModeQuad_ROS_B.absxr = sqrt(x->im / 2.0);
      absxi = multiModeQuad_ROS_B.absxr;
    }
  } else if (rtIsNaN(x->re)) {
    multiModeQuad_ROS_B.absxr = x->re;
    absxi = x->re;
  } else if (rtIsNaN(x->im)) {
    multiModeQuad_ROS_B.absxr = x->im;
    absxi = x->im;
  } else if (rtIsInf(x->im)) {
    multiModeQuad_ROS_B.absxr = fabs(x->im);
    absxi = x->im;
  } else if (rtIsInf(x->re)) {
    if (x->re < 0.0) {
      multiModeQuad_ROS_B.absxr = 0.0;
      absxi = x->im * -x->re;
    } else {
      multiModeQuad_ROS_B.absxr = x->re;
      absxi = 0.0;
    }
  } else {
    multiModeQuad_ROS_B.absxr = fabs(x->re);
    absxi = fabs(x->im);
    if ((multiModeQuad_ROS_B.absxr > 4.4942328371557893E+307) || (absxi >
         4.4942328371557893E+307)) {
      multiModeQuad_ROS_B.absxr *= 0.5;
      absxi = multiModeQuad_ROS_rt_hypotd_snf(multiModeQuad_ROS_B.absxr, absxi *
        0.5);
      if (absxi > multiModeQuad_ROS_B.absxr) {
        multiModeQuad_ROS_B.absxr = sqrt(multiModeQuad_ROS_B.absxr / absxi + 1.0)
          * sqrt(absxi);
      } else {
        multiModeQuad_ROS_B.absxr = sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      multiModeQuad_ROS_B.absxr = sqrt((multiModeQuad_ROS_rt_hypotd_snf
        (multiModeQuad_ROS_B.absxr, absxi) + multiModeQuad_ROS_B.absxr) * 0.5);
    }

    if (x->re > 0.0) {
      absxi = x->im / multiModeQuad_ROS_B.absxr * 0.5;
    } else {
      if (x->im < 0.0) {
        absxi = -multiModeQuad_ROS_B.absxr;
      } else {
        absxi = multiModeQuad_ROS_B.absxr;
      }

      multiModeQuad_ROS_B.absxr = x->im / absxi * 0.5;
    }
  }

  x->re = multiModeQuad_ROS_B.absxr;
  x->im = absxi;
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_xzlartg_h(const creal_T f, const creal_T g, real_T
  *cs, creal_T *sn)
{
  real_T gs_im;
  boolean_T guard1 = false;
  multiModeQuad_ROS_B.d = fabs(f.re);
  multiModeQuad_ROS_B.scale_dh = multiModeQuad_ROS_B.d;
  multiModeQuad_ROS_B.f2s_l = fabs(f.im);
  if (multiModeQuad_ROS_B.f2s_l > multiModeQuad_ROS_B.d) {
    multiModeQuad_ROS_B.scale_dh = multiModeQuad_ROS_B.f2s_l;
  }

  multiModeQuad_ROS_B.gs_re_b = fabs(g.re);
  gs_im = fabs(g.im);
  if (gs_im > multiModeQuad_ROS_B.gs_re_b) {
    multiModeQuad_ROS_B.gs_re_b = gs_im;
  }

  if (multiModeQuad_ROS_B.gs_re_b > multiModeQuad_ROS_B.scale_dh) {
    multiModeQuad_ROS_B.scale_dh = multiModeQuad_ROS_B.gs_re_b;
  }

  multiModeQuad_ROS_B.fs_re_l = f.re;
  multiModeQuad_ROS_B.fs_im_o = f.im;
  multiModeQuad_ROS_B.gs_re_b = g.re;
  gs_im = g.im;
  guard1 = false;
  if (multiModeQuad_ROS_B.scale_dh >= 7.4428285367870146E+137) {
    do {
      multiModeQuad_ROS_B.fs_re_l *= 1.3435752215134178E-138;
      multiModeQuad_ROS_B.fs_im_o *= 1.3435752215134178E-138;
      multiModeQuad_ROS_B.gs_re_b *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      multiModeQuad_ROS_B.scale_dh *= 1.3435752215134178E-138;
    } while (!(multiModeQuad_ROS_B.scale_dh < 7.4428285367870146E+137));

    guard1 = true;
  } else if (multiModeQuad_ROS_B.scale_dh <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
    } else {
      do {
        multiModeQuad_ROS_B.fs_re_l *= 7.4428285367870146E+137;
        multiModeQuad_ROS_B.fs_im_o *= 7.4428285367870146E+137;
        multiModeQuad_ROS_B.gs_re_b *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        multiModeQuad_ROS_B.scale_dh *= 7.4428285367870146E+137;
      } while (!(multiModeQuad_ROS_B.scale_dh > 1.3435752215134178E-138));

      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    multiModeQuad_ROS_B.scale_dh = multiModeQuad_ROS_B.fs_re_l *
      multiModeQuad_ROS_B.fs_re_l + multiModeQuad_ROS_B.fs_im_o *
      multiModeQuad_ROS_B.fs_im_o;
    multiModeQuad_ROS_B.g2_g = multiModeQuad_ROS_B.gs_re_b *
      multiModeQuad_ROS_B.gs_re_b + gs_im * gs_im;
    multiModeQuad_ROS_B.x_d = multiModeQuad_ROS_B.g2_g;
    if (1.0 > multiModeQuad_ROS_B.g2_g) {
      multiModeQuad_ROS_B.x_d = 1.0;
    }

    if (multiModeQuad_ROS_B.scale_dh <= multiModeQuad_ROS_B.x_d *
        2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        multiModeQuad_ROS_B.d = multiModeQuad_ROS_rt_hypotd_snf
          (multiModeQuad_ROS_B.gs_re_b, gs_im);
        sn->re = multiModeQuad_ROS_B.gs_re_b / multiModeQuad_ROS_B.d;
        sn->im = -gs_im / multiModeQuad_ROS_B.d;
      } else {
        multiModeQuad_ROS_B.scale_dh = sqrt(multiModeQuad_ROS_B.g2_g);
        *cs = multiModeQuad_ROS_rt_hypotd_snf(multiModeQuad_ROS_B.fs_re_l,
          multiModeQuad_ROS_B.fs_im_o) / multiModeQuad_ROS_B.scale_dh;
        if (multiModeQuad_ROS_B.f2s_l > multiModeQuad_ROS_B.d) {
          multiModeQuad_ROS_B.d = multiModeQuad_ROS_B.f2s_l;
        }

        if (multiModeQuad_ROS_B.d > 1.0) {
          multiModeQuad_ROS_B.d = multiModeQuad_ROS_rt_hypotd_snf(f.re, f.im);
          multiModeQuad_ROS_B.fs_re_l = f.re / multiModeQuad_ROS_B.d;
          multiModeQuad_ROS_B.fs_im_o = f.im / multiModeQuad_ROS_B.d;
        } else {
          multiModeQuad_ROS_B.fs_re_l = 7.4428285367870146E+137 * f.re;
          multiModeQuad_ROS_B.f2s_l = 7.4428285367870146E+137 * f.im;
          multiModeQuad_ROS_B.d = multiModeQuad_ROS_rt_hypotd_snf
            (multiModeQuad_ROS_B.fs_re_l, multiModeQuad_ROS_B.f2s_l);
          multiModeQuad_ROS_B.fs_re_l /= multiModeQuad_ROS_B.d;
          multiModeQuad_ROS_B.fs_im_o = multiModeQuad_ROS_B.f2s_l /
            multiModeQuad_ROS_B.d;
        }

        multiModeQuad_ROS_B.gs_re_b /= multiModeQuad_ROS_B.scale_dh;
        gs_im = -gs_im / multiModeQuad_ROS_B.scale_dh;
        sn->re = multiModeQuad_ROS_B.fs_re_l * multiModeQuad_ROS_B.gs_re_b -
          multiModeQuad_ROS_B.fs_im_o * gs_im;
        sn->im = multiModeQuad_ROS_B.fs_re_l * gs_im +
          multiModeQuad_ROS_B.fs_im_o * multiModeQuad_ROS_B.gs_re_b;
      }
    } else {
      multiModeQuad_ROS_B.f2s_l = sqrt(multiModeQuad_ROS_B.g2_g /
        multiModeQuad_ROS_B.scale_dh + 1.0);
      *cs = 1.0 / multiModeQuad_ROS_B.f2s_l;
      multiModeQuad_ROS_B.d = multiModeQuad_ROS_B.scale_dh +
        multiModeQuad_ROS_B.g2_g;
      multiModeQuad_ROS_B.fs_re_l = multiModeQuad_ROS_B.f2s_l *
        multiModeQuad_ROS_B.fs_re_l / multiModeQuad_ROS_B.d;
      multiModeQuad_ROS_B.fs_im_o = multiModeQuad_ROS_B.f2s_l *
        multiModeQuad_ROS_B.fs_im_o / multiModeQuad_ROS_B.d;
      sn->re = multiModeQuad_ROS_B.fs_re_l * multiModeQuad_ROS_B.gs_re_b -
        multiModeQuad_ROS_B.fs_im_o * -gs_im;
      sn->im = multiModeQuad_ROS_B.fs_re_l * -gs_im +
        multiModeQuad_ROS_B.fs_im_o * multiModeQuad_ROS_B.gs_re_b;
    }
  }
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_xzlartg(const creal_T f, const creal_T g, real_T
  *cs, creal_T *sn, creal_T *r)
{
  int32_T count;
  int32_T rescaledir;
  boolean_T guard1 = false;
  multiModeQuad_ROS_B.f2s = fabs(f.re);
  multiModeQuad_ROS_B.scale_n = multiModeQuad_ROS_B.f2s;
  multiModeQuad_ROS_B.di = fabs(f.im);
  if (multiModeQuad_ROS_B.di > multiModeQuad_ROS_B.f2s) {
    multiModeQuad_ROS_B.scale_n = multiModeQuad_ROS_B.di;
  }

  multiModeQuad_ROS_B.gs_re = fabs(g.re);
  multiModeQuad_ROS_B.gs_im = fabs(g.im);
  if (multiModeQuad_ROS_B.gs_im > multiModeQuad_ROS_B.gs_re) {
    multiModeQuad_ROS_B.gs_re = multiModeQuad_ROS_B.gs_im;
  }

  if (multiModeQuad_ROS_B.gs_re > multiModeQuad_ROS_B.scale_n) {
    multiModeQuad_ROS_B.scale_n = multiModeQuad_ROS_B.gs_re;
  }

  multiModeQuad_ROS_B.fs_re = f.re;
  multiModeQuad_ROS_B.fs_im = f.im;
  multiModeQuad_ROS_B.gs_re = g.re;
  multiModeQuad_ROS_B.gs_im = g.im;
  count = -1;
  rescaledir = 0;
  guard1 = false;
  if (multiModeQuad_ROS_B.scale_n >= 7.4428285367870146E+137) {
    do {
      count++;
      multiModeQuad_ROS_B.fs_re *= 1.3435752215134178E-138;
      multiModeQuad_ROS_B.fs_im *= 1.3435752215134178E-138;
      multiModeQuad_ROS_B.gs_re *= 1.3435752215134178E-138;
      multiModeQuad_ROS_B.gs_im *= 1.3435752215134178E-138;
      multiModeQuad_ROS_B.scale_n *= 1.3435752215134178E-138;
    } while (!(multiModeQuad_ROS_B.scale_n < 7.4428285367870146E+137));

    rescaledir = 1;
    guard1 = true;
  } else if (multiModeQuad_ROS_B.scale_n <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
      *r = f;
    } else {
      do {
        count++;
        multiModeQuad_ROS_B.fs_re *= 7.4428285367870146E+137;
        multiModeQuad_ROS_B.fs_im *= 7.4428285367870146E+137;
        multiModeQuad_ROS_B.gs_re *= 7.4428285367870146E+137;
        multiModeQuad_ROS_B.gs_im *= 7.4428285367870146E+137;
        multiModeQuad_ROS_B.scale_n *= 7.4428285367870146E+137;
      } while (!(multiModeQuad_ROS_B.scale_n > 1.3435752215134178E-138));

      rescaledir = -1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    multiModeQuad_ROS_B.scale_n = multiModeQuad_ROS_B.fs_re *
      multiModeQuad_ROS_B.fs_re + multiModeQuad_ROS_B.fs_im *
      multiModeQuad_ROS_B.fs_im;
    multiModeQuad_ROS_B.g2 = multiModeQuad_ROS_B.gs_re *
      multiModeQuad_ROS_B.gs_re + multiModeQuad_ROS_B.gs_im *
      multiModeQuad_ROS_B.gs_im;
    multiModeQuad_ROS_B.x = multiModeQuad_ROS_B.g2;
    if (1.0 > multiModeQuad_ROS_B.g2) {
      multiModeQuad_ROS_B.x = 1.0;
    }

    if (multiModeQuad_ROS_B.scale_n <= multiModeQuad_ROS_B.x *
        2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        r->re = multiModeQuad_ROS_rt_hypotd_snf(g.re, g.im);
        r->im = 0.0;
        multiModeQuad_ROS_B.f2s = multiModeQuad_ROS_rt_hypotd_snf
          (multiModeQuad_ROS_B.gs_re, multiModeQuad_ROS_B.gs_im);
        sn->re = multiModeQuad_ROS_B.gs_re / multiModeQuad_ROS_B.f2s;
        sn->im = -multiModeQuad_ROS_B.gs_im / multiModeQuad_ROS_B.f2s;
      } else {
        multiModeQuad_ROS_B.scale_n = sqrt(multiModeQuad_ROS_B.g2);
        *cs = multiModeQuad_ROS_rt_hypotd_snf(multiModeQuad_ROS_B.fs_re,
          multiModeQuad_ROS_B.fs_im) / multiModeQuad_ROS_B.scale_n;
        if (multiModeQuad_ROS_B.di > multiModeQuad_ROS_B.f2s) {
          multiModeQuad_ROS_B.f2s = multiModeQuad_ROS_B.di;
        }

        if (multiModeQuad_ROS_B.f2s > 1.0) {
          multiModeQuad_ROS_B.f2s = multiModeQuad_ROS_rt_hypotd_snf(f.re, f.im);
          multiModeQuad_ROS_B.fs_re = f.re / multiModeQuad_ROS_B.f2s;
          multiModeQuad_ROS_B.fs_im = f.im / multiModeQuad_ROS_B.f2s;
        } else {
          multiModeQuad_ROS_B.fs_re = 7.4428285367870146E+137 * f.re;
          multiModeQuad_ROS_B.di = 7.4428285367870146E+137 * f.im;
          multiModeQuad_ROS_B.f2s = multiModeQuad_ROS_rt_hypotd_snf
            (multiModeQuad_ROS_B.fs_re, multiModeQuad_ROS_B.di);
          multiModeQuad_ROS_B.fs_re /= multiModeQuad_ROS_B.f2s;
          multiModeQuad_ROS_B.fs_im = multiModeQuad_ROS_B.di /
            multiModeQuad_ROS_B.f2s;
        }

        multiModeQuad_ROS_B.gs_re /= multiModeQuad_ROS_B.scale_n;
        multiModeQuad_ROS_B.gs_im = -multiModeQuad_ROS_B.gs_im /
          multiModeQuad_ROS_B.scale_n;
        sn->re = multiModeQuad_ROS_B.fs_re * multiModeQuad_ROS_B.gs_re -
          multiModeQuad_ROS_B.fs_im * multiModeQuad_ROS_B.gs_im;
        sn->im = multiModeQuad_ROS_B.fs_re * multiModeQuad_ROS_B.gs_im +
          multiModeQuad_ROS_B.fs_im * multiModeQuad_ROS_B.gs_re;
        r->re = (sn->re * g.re - sn->im * g.im) + *cs * f.re;
        r->im = (sn->re * g.im + sn->im * g.re) + *cs * f.im;
      }
    } else {
      multiModeQuad_ROS_B.f2s = sqrt(multiModeQuad_ROS_B.g2 /
        multiModeQuad_ROS_B.scale_n + 1.0);
      r->re = multiModeQuad_ROS_B.f2s * multiModeQuad_ROS_B.fs_re;
      r->im = multiModeQuad_ROS_B.f2s * multiModeQuad_ROS_B.fs_im;
      *cs = 1.0 / multiModeQuad_ROS_B.f2s;
      multiModeQuad_ROS_B.f2s = multiModeQuad_ROS_B.scale_n +
        multiModeQuad_ROS_B.g2;
      multiModeQuad_ROS_B.fs_re = r->re / multiModeQuad_ROS_B.f2s;
      multiModeQuad_ROS_B.f2s = r->im / multiModeQuad_ROS_B.f2s;
      sn->re = multiModeQuad_ROS_B.fs_re * multiModeQuad_ROS_B.gs_re -
        multiModeQuad_ROS_B.f2s * -multiModeQuad_ROS_B.gs_im;
      sn->im = multiModeQuad_ROS_B.fs_re * -multiModeQuad_ROS_B.gs_im +
        multiModeQuad_ROS_B.f2s * multiModeQuad_ROS_B.gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 7.4428285367870146E+137;
          r->im *= 7.4428285367870146E+137;
        }
      } else if (rescaledir < 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 1.3435752215134178E-138;
          r->im *= 1.3435752215134178E-138;
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_xzhgeqz(creal_T A[16], int32_T ilo, int32_T ihi,
  creal_T Z[16], int32_T *info, creal_T alpha1[4], creal_T beta1[4])
{
  int32_T absxk_tmp;
  int32_T col;
  int32_T iiter;
  int32_T ilastm1;
  int32_T jp1;
  int32_T nm1;
  boolean_T failed;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  *info = 0;
  alpha1[0].re = 0.0;
  alpha1[0].im = 0.0;
  beta1[0].re = 1.0;
  beta1[0].im = 0.0;
  alpha1[1].re = 0.0;
  alpha1[1].im = 0.0;
  beta1[1].re = 1.0;
  beta1[1].im = 0.0;
  alpha1[2].re = 0.0;
  alpha1[2].im = 0.0;
  beta1[2].re = 1.0;
  beta1[2].im = 0.0;
  alpha1[3].re = 0.0;
  alpha1[3].im = 0.0;
  beta1[3].re = 1.0;
  beta1[3].im = 0.0;
  multiModeQuad_ROS_B.eshift_re = 0.0;
  multiModeQuad_ROS_B.eshift_im = 0.0;
  multiModeQuad_ROS_B.ctemp.re = 0.0;
  multiModeQuad_ROS_B.ctemp.im = 0.0;
  multiModeQuad_ROS_B.anorm = 0.0;
  if (ilo <= ihi) {
    multiModeQuad_ROS_B.scale = 3.3121686421112381E-170;
    multiModeQuad_ROS_B.ssq = 0.0;
    nm1 = ihi - ilo;
    multiModeQuad_ROS_B.ifirst = -1;
    while (multiModeQuad_ROS_B.ifirst + 1 <= nm1) {
      multiModeQuad_ROS_B.colscale = 3.3121686421112381E-170;
      multiModeQuad_ROS_B.anorm = 0.0;
      col = ilo + multiModeQuad_ROS_B.ifirst;
      if (multiModeQuad_ROS_B.ifirst + 2 <= nm1) {
        ilastm1 = multiModeQuad_ROS_B.ifirst + 2;
      } else {
        ilastm1 = nm1;
      }

      ilastm1 += ilo;
      for (iiter = ilo; iiter <= ilastm1; iiter++) {
        absxk_tmp = ((col << 2) + iiter) - 1;
        multiModeQuad_ROS_B.absxk = fabs(A[absxk_tmp].re);
        if (multiModeQuad_ROS_B.absxk > multiModeQuad_ROS_B.colscale) {
          multiModeQuad_ROS_B.t = multiModeQuad_ROS_B.colscale /
            multiModeQuad_ROS_B.absxk;
          multiModeQuad_ROS_B.anorm = multiModeQuad_ROS_B.anorm *
            multiModeQuad_ROS_B.t * multiModeQuad_ROS_B.t + 1.0;
          multiModeQuad_ROS_B.colscale = multiModeQuad_ROS_B.absxk;
        } else {
          multiModeQuad_ROS_B.t = multiModeQuad_ROS_B.absxk /
            multiModeQuad_ROS_B.colscale;
          multiModeQuad_ROS_B.anorm += multiModeQuad_ROS_B.t *
            multiModeQuad_ROS_B.t;
        }

        multiModeQuad_ROS_B.absxk = fabs(A[absxk_tmp].im);
        if (multiModeQuad_ROS_B.absxk > multiModeQuad_ROS_B.colscale) {
          multiModeQuad_ROS_B.t = multiModeQuad_ROS_B.colscale /
            multiModeQuad_ROS_B.absxk;
          multiModeQuad_ROS_B.anorm = multiModeQuad_ROS_B.anorm *
            multiModeQuad_ROS_B.t * multiModeQuad_ROS_B.t + 1.0;
          multiModeQuad_ROS_B.colscale = multiModeQuad_ROS_B.absxk;
        } else {
          multiModeQuad_ROS_B.t = multiModeQuad_ROS_B.absxk /
            multiModeQuad_ROS_B.colscale;
          multiModeQuad_ROS_B.anorm += multiModeQuad_ROS_B.t *
            multiModeQuad_ROS_B.t;
        }
      }

      if (multiModeQuad_ROS_B.scale >= multiModeQuad_ROS_B.colscale) {
        multiModeQuad_ROS_B.colscale /= multiModeQuad_ROS_B.scale;
        multiModeQuad_ROS_B.ssq += multiModeQuad_ROS_B.colscale *
          multiModeQuad_ROS_B.colscale * multiModeQuad_ROS_B.anorm;
      } else {
        multiModeQuad_ROS_B.scale /= multiModeQuad_ROS_B.colscale;
        multiModeQuad_ROS_B.ssq = multiModeQuad_ROS_B.scale *
          multiModeQuad_ROS_B.scale * multiModeQuad_ROS_B.ssq +
          multiModeQuad_ROS_B.anorm;
        multiModeQuad_ROS_B.scale = multiModeQuad_ROS_B.colscale;
      }

      multiModeQuad_ROS_B.ifirst++;
    }

    multiModeQuad_ROS_B.anorm = multiModeQuad_ROS_B.scale * sqrt
      (multiModeQuad_ROS_B.ssq);
  }

  multiModeQuad_ROS_B.ssq = 2.2250738585072014E-308;
  multiModeQuad_ROS_B.scale = 2.2204460492503131E-16 * multiModeQuad_ROS_B.anorm;
  if (multiModeQuad_ROS_B.scale > 2.2250738585072014E-308) {
    multiModeQuad_ROS_B.ssq = multiModeQuad_ROS_B.scale;
  }

  multiModeQuad_ROS_B.scale = 2.2250738585072014E-308;
  if (multiModeQuad_ROS_B.anorm > 2.2250738585072014E-308) {
    multiModeQuad_ROS_B.scale = multiModeQuad_ROS_B.anorm;
  }

  multiModeQuad_ROS_B.anorm = 1.0 / multiModeQuad_ROS_B.scale;
  failed = true;
  for (nm1 = ihi; nm1 + 1 < 5; nm1++) {
    alpha1[nm1] = A[(nm1 << 2) + nm1];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    boolean_T goto60;
    boolean_T goto70;
    boolean_T goto90;
    multiModeQuad_ROS_B.ifirst = ilo;
    col = ilo;
    nm1 = ihi - 1;
    ilastm1 = ihi - 2;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    absxk_tmp = 0;
    int32_T exitg1;
    do {
      exitg1 = 0;
      if (absxk_tmp <= ((ihi - ilo) + 1) * 30 - 1) {
        boolean_T exitg2;
        if (nm1 + 1 == ilo) {
          goto60 = true;
        } else {
          jp1 = (ilastm1 << 2) + nm1;
          if (fabs(A[jp1].re) + fabs(A[jp1].im) <= multiModeQuad_ROS_B.ssq) {
            A[jp1].re = 0.0;
            A[jp1].im = 0.0;
            goto60 = true;
          } else {
            boolean_T guard3 = false;
            multiModeQuad_ROS_B.j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (multiModeQuad_ROS_B.j + 1 >= ilo)) {
              if (multiModeQuad_ROS_B.j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else {
                jp1 = ((multiModeQuad_ROS_B.j - 1) << 2) + multiModeQuad_ROS_B.j;
                if (fabs(A[jp1].re) + fabs(A[jp1].im) <= multiModeQuad_ROS_B.ssq)
                {
                  A[jp1].re = 0.0;
                  A[jp1].im = 0.0;
                  guard3 = true;
                  exitg2 = true;
                } else {
                  multiModeQuad_ROS_B.j--;
                  guard3 = false;
                }
              }
            }

            if (guard3) {
              multiModeQuad_ROS_B.ifirst = multiModeQuad_ROS_B.j + 1;
              goto70 = true;
            }
          }
        }

        if ((!goto60) && (!goto70)) {
          alpha1[0].re = (rtNaN);
          alpha1[0].im = 0.0;
          beta1[0].re = (rtNaN);
          beta1[0].im = 0.0;
          alpha1[1].re = (rtNaN);
          alpha1[1].im = 0.0;
          beta1[1].re = (rtNaN);
          beta1[1].im = 0.0;
          alpha1[2].re = (rtNaN);
          alpha1[2].im = 0.0;
          beta1[2].re = (rtNaN);
          beta1[2].im = 0.0;
          alpha1[3].re = (rtNaN);
          alpha1[3].im = 0.0;
          beta1[3].re = (rtNaN);
          beta1[3].im = 0.0;
          for (jp1 = 0; jp1 < 16; jp1++) {
            Z[jp1].re = (rtNaN);
            Z[jp1].im = 0.0;
          }

          *info = 1;
          exitg1 = 1;
        } else if (goto60) {
          goto60 = false;
          alpha1[nm1] = A[(nm1 << 2) + nm1];
          nm1 = ilastm1;
          ilastm1--;
          if (nm1 + 1 < ilo) {
            failed = false;
            guard2 = true;
            exitg1 = 1;
          } else {
            iiter = 0;
            multiModeQuad_ROS_B.eshift_re = 0.0;
            multiModeQuad_ROS_B.eshift_im = 0.0;
            absxk_tmp++;
          }
        } else {
          if (goto70) {
            int32_T ctemp_tmp;
            int32_T ctemp_tmp_tmp;
            goto70 = false;
            iiter++;
            if (iiter - div_nzp_s32(iiter, 10) * 10 != 0) {
              multiModeQuad_ROS_B.j = (ilastm1 << 2) + ilastm1;
              multiModeQuad_ROS_B.ar = A[multiModeQuad_ROS_B.j].re *
                multiModeQuad_ROS_B.anorm;
              multiModeQuad_ROS_B.ai = A[multiModeQuad_ROS_B.j].im *
                multiModeQuad_ROS_B.anorm;
              if (multiModeQuad_ROS_B.ai == 0.0) {
                multiModeQuad_ROS_B.shift.re = multiModeQuad_ROS_B.ar / 0.5;
                multiModeQuad_ROS_B.shift.im = 0.0;
              } else if (multiModeQuad_ROS_B.ar == 0.0) {
                multiModeQuad_ROS_B.shift.re = 0.0;
                multiModeQuad_ROS_B.shift.im = multiModeQuad_ROS_B.ai / 0.5;
              } else {
                multiModeQuad_ROS_B.shift.re = multiModeQuad_ROS_B.ar / 0.5;
                multiModeQuad_ROS_B.shift.im = multiModeQuad_ROS_B.ai / 0.5;
              }

              multiModeQuad_ROS_B.j = (nm1 << 2) + nm1;
              multiModeQuad_ROS_B.ar = A[multiModeQuad_ROS_B.j].re *
                multiModeQuad_ROS_B.anorm;
              multiModeQuad_ROS_B.ai = A[multiModeQuad_ROS_B.j].im *
                multiModeQuad_ROS_B.anorm;
              if (multiModeQuad_ROS_B.ai == 0.0) {
                multiModeQuad_ROS_B.absxk = multiModeQuad_ROS_B.ar / 0.5;
                multiModeQuad_ROS_B.t = 0.0;
              } else if (multiModeQuad_ROS_B.ar == 0.0) {
                multiModeQuad_ROS_B.absxk = 0.0;
                multiModeQuad_ROS_B.t = multiModeQuad_ROS_B.ai / 0.5;
              } else {
                multiModeQuad_ROS_B.absxk = multiModeQuad_ROS_B.ar / 0.5;
                multiModeQuad_ROS_B.t = multiModeQuad_ROS_B.ai / 0.5;
              }

              multiModeQuad_ROS_B.t1_re = (multiModeQuad_ROS_B.shift.re +
                multiModeQuad_ROS_B.absxk) * 0.5;
              multiModeQuad_ROS_B.t1_im = (multiModeQuad_ROS_B.shift.im +
                multiModeQuad_ROS_B.t) * 0.5;
              multiModeQuad_ROS_B.j = (nm1 << 2) + ilastm1;
              multiModeQuad_ROS_B.ar = A[multiModeQuad_ROS_B.j].re *
                multiModeQuad_ROS_B.anorm;
              multiModeQuad_ROS_B.ai = A[multiModeQuad_ROS_B.j].im *
                multiModeQuad_ROS_B.anorm;
              if (multiModeQuad_ROS_B.ai == 0.0) {
                multiModeQuad_ROS_B.scale = multiModeQuad_ROS_B.ar / 0.5;
                multiModeQuad_ROS_B.colscale = 0.0;
              } else if (multiModeQuad_ROS_B.ar == 0.0) {
                multiModeQuad_ROS_B.scale = 0.0;
                multiModeQuad_ROS_B.colscale = multiModeQuad_ROS_B.ai / 0.5;
              } else {
                multiModeQuad_ROS_B.scale = multiModeQuad_ROS_B.ar / 0.5;
                multiModeQuad_ROS_B.colscale = multiModeQuad_ROS_B.ai / 0.5;
              }

              multiModeQuad_ROS_B.j = (ilastm1 << 2) + nm1;
              multiModeQuad_ROS_B.ar = A[multiModeQuad_ROS_B.j].re *
                multiModeQuad_ROS_B.anorm;
              multiModeQuad_ROS_B.ai = A[multiModeQuad_ROS_B.j].im *
                multiModeQuad_ROS_B.anorm;
              if (multiModeQuad_ROS_B.ai == 0.0) {
                multiModeQuad_ROS_B.ar /= 0.5;
                multiModeQuad_ROS_B.ai = 0.0;
              } else if (multiModeQuad_ROS_B.ar == 0.0) {
                multiModeQuad_ROS_B.ar = 0.0;
                multiModeQuad_ROS_B.ai /= 0.5;
              } else {
                multiModeQuad_ROS_B.ar /= 0.5;
                multiModeQuad_ROS_B.ai /= 0.5;
              }

              multiModeQuad_ROS_B.shift_im = multiModeQuad_ROS_B.shift.re *
                multiModeQuad_ROS_B.t + multiModeQuad_ROS_B.shift.im *
                multiModeQuad_ROS_B.absxk;
              multiModeQuad_ROS_B.shift.re = ((multiModeQuad_ROS_B.t1_re *
                multiModeQuad_ROS_B.t1_re - multiModeQuad_ROS_B.t1_im *
                multiModeQuad_ROS_B.t1_im) + (multiModeQuad_ROS_B.scale *
                multiModeQuad_ROS_B.ar - multiModeQuad_ROS_B.colscale *
                multiModeQuad_ROS_B.ai)) - (multiModeQuad_ROS_B.shift.re *
                multiModeQuad_ROS_B.absxk - multiModeQuad_ROS_B.shift.im *
                multiModeQuad_ROS_B.t);
              multiModeQuad_ROS_B.shift_tmp = multiModeQuad_ROS_B.t1_re *
                multiModeQuad_ROS_B.t1_im;
              multiModeQuad_ROS_B.shift.im = ((multiModeQuad_ROS_B.scale *
                multiModeQuad_ROS_B.ai + multiModeQuad_ROS_B.colscale *
                multiModeQuad_ROS_B.ar) + (multiModeQuad_ROS_B.shift_tmp +
                multiModeQuad_ROS_B.shift_tmp)) - multiModeQuad_ROS_B.shift_im;
              multiModeQuad_ROS_sqrt(&multiModeQuad_ROS_B.shift);
              if ((multiModeQuad_ROS_B.t1_re - multiModeQuad_ROS_B.absxk) *
                  multiModeQuad_ROS_B.shift.re + (multiModeQuad_ROS_B.t1_im -
                   multiModeQuad_ROS_B.t) * multiModeQuad_ROS_B.shift.im <= 0.0)
              {
                multiModeQuad_ROS_B.shift.re += multiModeQuad_ROS_B.t1_re;
                multiModeQuad_ROS_B.shift.im += multiModeQuad_ROS_B.t1_im;
              } else {
                multiModeQuad_ROS_B.shift.re = multiModeQuad_ROS_B.t1_re -
                  multiModeQuad_ROS_B.shift.re;
                multiModeQuad_ROS_B.shift.im = multiModeQuad_ROS_B.t1_im -
                  multiModeQuad_ROS_B.shift.im;
              }
            } else {
              multiModeQuad_ROS_B.j = (ilastm1 << 2) + nm1;
              multiModeQuad_ROS_B.ar = A[multiModeQuad_ROS_B.j].re *
                multiModeQuad_ROS_B.anorm;
              multiModeQuad_ROS_B.ai = A[multiModeQuad_ROS_B.j].im *
                multiModeQuad_ROS_B.anorm;
              if (multiModeQuad_ROS_B.ai == 0.0) {
                multiModeQuad_ROS_B.scale = multiModeQuad_ROS_B.ar / 0.5;
                multiModeQuad_ROS_B.colscale = 0.0;
              } else if (multiModeQuad_ROS_B.ar == 0.0) {
                multiModeQuad_ROS_B.scale = 0.0;
                multiModeQuad_ROS_B.colscale = multiModeQuad_ROS_B.ai / 0.5;
              } else {
                multiModeQuad_ROS_B.scale = multiModeQuad_ROS_B.ar / 0.5;
                multiModeQuad_ROS_B.colscale = multiModeQuad_ROS_B.ai / 0.5;
              }

              multiModeQuad_ROS_B.eshift_re += multiModeQuad_ROS_B.scale;
              multiModeQuad_ROS_B.eshift_im += multiModeQuad_ROS_B.colscale;
              multiModeQuad_ROS_B.shift.re = multiModeQuad_ROS_B.eshift_re;
              multiModeQuad_ROS_B.shift.im = multiModeQuad_ROS_B.eshift_im;
            }

            multiModeQuad_ROS_B.j = ilastm1;
            jp1 = ilastm1 + 1;
            exitg2 = false;
            while ((!exitg2) && (multiModeQuad_ROS_B.j + 1 >
                                 multiModeQuad_ROS_B.ifirst)) {
              col = multiModeQuad_ROS_B.j + 1;
              ctemp_tmp_tmp = multiModeQuad_ROS_B.j << 2;
              ctemp_tmp = ctemp_tmp_tmp + multiModeQuad_ROS_B.j;
              multiModeQuad_ROS_B.ctemp.re = A[ctemp_tmp].re *
                multiModeQuad_ROS_B.anorm - multiModeQuad_ROS_B.shift.re * 0.5;
              multiModeQuad_ROS_B.ctemp.im = A[ctemp_tmp].im *
                multiModeQuad_ROS_B.anorm - multiModeQuad_ROS_B.shift.im * 0.5;
              multiModeQuad_ROS_B.scale = fabs(multiModeQuad_ROS_B.ctemp.re) +
                fabs(multiModeQuad_ROS_B.ctemp.im);
              jp1 += ctemp_tmp_tmp;
              multiModeQuad_ROS_B.colscale = (fabs(A[jp1].re) + fabs(A[jp1].im))
                * multiModeQuad_ROS_B.anorm;
              multiModeQuad_ROS_B.absxk = multiModeQuad_ROS_B.scale;
              if (multiModeQuad_ROS_B.colscale > multiModeQuad_ROS_B.scale) {
                multiModeQuad_ROS_B.absxk = multiModeQuad_ROS_B.colscale;
              }

              if ((multiModeQuad_ROS_B.absxk < 1.0) &&
                  (multiModeQuad_ROS_B.absxk != 0.0)) {
                multiModeQuad_ROS_B.scale /= multiModeQuad_ROS_B.absxk;
                multiModeQuad_ROS_B.colscale /= multiModeQuad_ROS_B.absxk;
              }

              jp1 = ((multiModeQuad_ROS_B.j - 1) << 2) + multiModeQuad_ROS_B.j;
              if ((fabs(A[jp1].re) + fabs(A[jp1].im)) *
                  multiModeQuad_ROS_B.colscale <= multiModeQuad_ROS_B.scale *
                  multiModeQuad_ROS_B.ssq) {
                goto90 = true;
                exitg2 = true;
              } else {
                jp1 = multiModeQuad_ROS_B.j;
                multiModeQuad_ROS_B.j--;
              }
            }

            if (!goto90) {
              col = multiModeQuad_ROS_B.ifirst;
              ctemp_tmp = (((multiModeQuad_ROS_B.ifirst - 1) << 2) +
                           multiModeQuad_ROS_B.ifirst) - 1;
              multiModeQuad_ROS_B.ctemp.re = A[ctemp_tmp].re *
                multiModeQuad_ROS_B.anorm - multiModeQuad_ROS_B.shift.re * 0.5;
              multiModeQuad_ROS_B.ctemp.im = A[ctemp_tmp].im *
                multiModeQuad_ROS_B.anorm - multiModeQuad_ROS_B.shift.im * 0.5;
            }

            goto90 = false;
            multiModeQuad_ROS_B.j = ((col - 1) << 2) + col;
            multiModeQuad_ROS_B.ascale.re = A[multiModeQuad_ROS_B.j].re *
              multiModeQuad_ROS_B.anorm;
            multiModeQuad_ROS_B.ascale.im = A[multiModeQuad_ROS_B.j].im *
              multiModeQuad_ROS_B.anorm;
            multiModeQuad_ROS_xzlartg_h(multiModeQuad_ROS_B.ctemp,
              multiModeQuad_ROS_B.ascale, &multiModeQuad_ROS_B.scale,
              &multiModeQuad_ROS_B.shift);
            multiModeQuad_ROS_B.j = col;
            jp1 = col - 2;
            while (multiModeQuad_ROS_B.j < nm1 + 1) {
              int32_T ad22_re_tmp;
              if (multiModeQuad_ROS_B.j > col) {
                multiModeQuad_ROS_xzlartg(A[(multiModeQuad_ROS_B.j + (jp1 << 2))
                  - 1], A[multiModeQuad_ROS_B.j + (jp1 << 2)],
                  &multiModeQuad_ROS_B.scale, &multiModeQuad_ROS_B.shift, &A
                  [(multiModeQuad_ROS_B.j + (jp1 << 2)) - 1]);
                jp1 = (jp1 << 2) + multiModeQuad_ROS_B.j;
                A[jp1].re = 0.0;
                A[jp1].im = 0.0;
              }

              for (ctemp_tmp_tmp = multiModeQuad_ROS_B.j - 1; ctemp_tmp_tmp + 1 <
                   5; ctemp_tmp_tmp++) {
                jp1 = (ctemp_tmp_tmp << 2) + multiModeQuad_ROS_B.j;
                multiModeQuad_ROS_B.ar = A[jp1].im;
                multiModeQuad_ROS_B.absxk = A[jp1].re;
                multiModeQuad_ROS_B.t = A[jp1 - 1].re;
                multiModeQuad_ROS_B.t1_re = A[jp1 - 1].im;
                A[jp1].re = multiModeQuad_ROS_B.absxk *
                  multiModeQuad_ROS_B.scale - (multiModeQuad_ROS_B.t *
                  multiModeQuad_ROS_B.shift.re + multiModeQuad_ROS_B.t1_re *
                  multiModeQuad_ROS_B.shift.im);
                A[jp1].im = A[jp1].im * multiModeQuad_ROS_B.scale -
                  (multiModeQuad_ROS_B.shift.re * multiModeQuad_ROS_B.t1_re -
                   multiModeQuad_ROS_B.shift.im * multiModeQuad_ROS_B.t);
                A[jp1 - 1].re = (multiModeQuad_ROS_B.absxk *
                                 multiModeQuad_ROS_B.shift.re -
                                 multiModeQuad_ROS_B.ar *
                                 multiModeQuad_ROS_B.shift.im) +
                  multiModeQuad_ROS_B.t * multiModeQuad_ROS_B.scale;
                A[jp1 - 1].im = (multiModeQuad_ROS_B.ar *
                                 multiModeQuad_ROS_B.shift.re +
                                 multiModeQuad_ROS_B.absxk *
                                 multiModeQuad_ROS_B.shift.im) +
                  multiModeQuad_ROS_B.t1_re * multiModeQuad_ROS_B.scale;
              }

              multiModeQuad_ROS_B.shift.re = -multiModeQuad_ROS_B.shift.re;
              multiModeQuad_ROS_B.shift.im = -multiModeQuad_ROS_B.shift.im;
              ctemp_tmp_tmp = multiModeQuad_ROS_B.j;
              if (nm1 + 1 < multiModeQuad_ROS_B.j + 2) {
                ctemp_tmp_tmp = nm1 - 1;
              }

              for (ctemp_tmp = 0; ctemp_tmp < ctemp_tmp_tmp + 2; ctemp_tmp++) {
                jp1 = ((multiModeQuad_ROS_B.j - 1) << 2) + ctemp_tmp;
                multiModeQuad_ROS_B.ar = A[jp1].im;
                multiModeQuad_ROS_B.absxk = A[jp1].re;
                ad22_re_tmp = (multiModeQuad_ROS_B.j << 2) + ctemp_tmp;
                multiModeQuad_ROS_B.t = A[ad22_re_tmp].re;
                multiModeQuad_ROS_B.t1_re = A[ad22_re_tmp].im;
                A[jp1].re = multiModeQuad_ROS_B.absxk *
                  multiModeQuad_ROS_B.scale - (multiModeQuad_ROS_B.t *
                  multiModeQuad_ROS_B.shift.re + multiModeQuad_ROS_B.t1_re *
                  multiModeQuad_ROS_B.shift.im);
                A[jp1].im = A[jp1].im * multiModeQuad_ROS_B.scale -
                  (multiModeQuad_ROS_B.shift.re * multiModeQuad_ROS_B.t1_re -
                   multiModeQuad_ROS_B.shift.im * multiModeQuad_ROS_B.t);
                A[ad22_re_tmp].re = (multiModeQuad_ROS_B.absxk *
                                     multiModeQuad_ROS_B.shift.re -
                                     multiModeQuad_ROS_B.ar *
                                     multiModeQuad_ROS_B.shift.im) +
                  multiModeQuad_ROS_B.t * multiModeQuad_ROS_B.scale;
                A[ad22_re_tmp].im = (multiModeQuad_ROS_B.ar *
                                     multiModeQuad_ROS_B.shift.re +
                                     multiModeQuad_ROS_B.absxk *
                                     multiModeQuad_ROS_B.shift.im) +
                  multiModeQuad_ROS_B.t1_re * multiModeQuad_ROS_B.scale;
              }

              jp1 = (multiModeQuad_ROS_B.j - 1) << 2;
              multiModeQuad_ROS_B.ar = Z[jp1].im;
              multiModeQuad_ROS_B.absxk = Z[jp1].re;
              ad22_re_tmp = multiModeQuad_ROS_B.j << 2;
              multiModeQuad_ROS_B.t = Z[ad22_re_tmp].re;
              multiModeQuad_ROS_B.t1_re = Z[ad22_re_tmp].im;
              Z[jp1].re = multiModeQuad_ROS_B.absxk * multiModeQuad_ROS_B.scale
                - (multiModeQuad_ROS_B.t * multiModeQuad_ROS_B.shift.re +
                   multiModeQuad_ROS_B.t1_re * multiModeQuad_ROS_B.shift.im);
              Z[jp1].im = Z[jp1].im * multiModeQuad_ROS_B.scale -
                (multiModeQuad_ROS_B.shift.re * multiModeQuad_ROS_B.t1_re -
                 multiModeQuad_ROS_B.shift.im * multiModeQuad_ROS_B.t);
              Z[ad22_re_tmp].re = (multiModeQuad_ROS_B.absxk *
                                   multiModeQuad_ROS_B.shift.re -
                                   multiModeQuad_ROS_B.ar *
                                   multiModeQuad_ROS_B.shift.im) +
                multiModeQuad_ROS_B.t * multiModeQuad_ROS_B.scale;
              Z[ad22_re_tmp].im = (multiModeQuad_ROS_B.ar *
                                   multiModeQuad_ROS_B.shift.re +
                                   multiModeQuad_ROS_B.absxk *
                                   multiModeQuad_ROS_B.shift.im) +
                multiModeQuad_ROS_B.t1_re * multiModeQuad_ROS_B.scale;
              multiModeQuad_ROS_B.ar = Z[jp1 + 1].im;
              multiModeQuad_ROS_B.absxk = Z[jp1 + 1].re;
              multiModeQuad_ROS_B.t = Z[ad22_re_tmp + 1].re;
              multiModeQuad_ROS_B.t1_re = Z[ad22_re_tmp + 1].im;
              Z[jp1 + 1].re = multiModeQuad_ROS_B.absxk *
                multiModeQuad_ROS_B.scale - (multiModeQuad_ROS_B.t *
                multiModeQuad_ROS_B.shift.re + multiModeQuad_ROS_B.t1_re *
                multiModeQuad_ROS_B.shift.im);
              Z[jp1 + 1].im = Z[jp1 + 1].im * multiModeQuad_ROS_B.scale -
                (multiModeQuad_ROS_B.shift.re * multiModeQuad_ROS_B.t1_re -
                 multiModeQuad_ROS_B.shift.im * multiModeQuad_ROS_B.t);
              Z[ad22_re_tmp + 1].re = (multiModeQuad_ROS_B.absxk *
                multiModeQuad_ROS_B.shift.re - multiModeQuad_ROS_B.ar *
                multiModeQuad_ROS_B.shift.im) + multiModeQuad_ROS_B.t *
                multiModeQuad_ROS_B.scale;
              Z[ad22_re_tmp + 1].im = (multiModeQuad_ROS_B.ar *
                multiModeQuad_ROS_B.shift.re + multiModeQuad_ROS_B.absxk *
                multiModeQuad_ROS_B.shift.im) + multiModeQuad_ROS_B.t1_re *
                multiModeQuad_ROS_B.scale;
              multiModeQuad_ROS_B.ar = Z[jp1 + 2].im;
              multiModeQuad_ROS_B.absxk = Z[jp1 + 2].re;
              multiModeQuad_ROS_B.t = Z[ad22_re_tmp + 2].re;
              multiModeQuad_ROS_B.t1_re = Z[ad22_re_tmp + 2].im;
              Z[jp1 + 2].re = multiModeQuad_ROS_B.absxk *
                multiModeQuad_ROS_B.scale - (multiModeQuad_ROS_B.t *
                multiModeQuad_ROS_B.shift.re + multiModeQuad_ROS_B.t1_re *
                multiModeQuad_ROS_B.shift.im);
              Z[jp1 + 2].im = Z[jp1 + 2].im * multiModeQuad_ROS_B.scale -
                (multiModeQuad_ROS_B.shift.re * multiModeQuad_ROS_B.t1_re -
                 multiModeQuad_ROS_B.shift.im * multiModeQuad_ROS_B.t);
              Z[ad22_re_tmp + 2].re = (multiModeQuad_ROS_B.absxk *
                multiModeQuad_ROS_B.shift.re - multiModeQuad_ROS_B.ar *
                multiModeQuad_ROS_B.shift.im) + multiModeQuad_ROS_B.t *
                multiModeQuad_ROS_B.scale;
              Z[ad22_re_tmp + 2].im = (multiModeQuad_ROS_B.ar *
                multiModeQuad_ROS_B.shift.re + multiModeQuad_ROS_B.absxk *
                multiModeQuad_ROS_B.shift.im) + multiModeQuad_ROS_B.t1_re *
                multiModeQuad_ROS_B.scale;
              multiModeQuad_ROS_B.ar = Z[jp1 + 3].im;
              multiModeQuad_ROS_B.absxk = Z[jp1 + 3].re;
              multiModeQuad_ROS_B.t = Z[ad22_re_tmp + 3].re;
              multiModeQuad_ROS_B.t1_re = Z[ad22_re_tmp + 3].im;
              Z[jp1 + 3].re = multiModeQuad_ROS_B.absxk *
                multiModeQuad_ROS_B.scale - (multiModeQuad_ROS_B.t *
                multiModeQuad_ROS_B.shift.re + multiModeQuad_ROS_B.t1_re *
                multiModeQuad_ROS_B.shift.im);
              Z[jp1 + 3].im = Z[jp1 + 3].im * multiModeQuad_ROS_B.scale -
                (multiModeQuad_ROS_B.shift.re * multiModeQuad_ROS_B.t1_re -
                 multiModeQuad_ROS_B.shift.im * multiModeQuad_ROS_B.t);
              Z[ad22_re_tmp + 3].re = (multiModeQuad_ROS_B.absxk *
                multiModeQuad_ROS_B.shift.re - multiModeQuad_ROS_B.ar *
                multiModeQuad_ROS_B.shift.im) + multiModeQuad_ROS_B.t *
                multiModeQuad_ROS_B.scale;
              Z[ad22_re_tmp + 3].im = (multiModeQuad_ROS_B.ar *
                multiModeQuad_ROS_B.shift.re + multiModeQuad_ROS_B.absxk *
                multiModeQuad_ROS_B.shift.im) + multiModeQuad_ROS_B.t1_re *
                multiModeQuad_ROS_B.scale;
              jp1 = multiModeQuad_ROS_B.j - 1;
              multiModeQuad_ROS_B.j++;
            }
          }

          absxk_tmp++;
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (failed) {
      *info = nm1 + 1;
      multiModeQuad_ROS_B.ifirst = 0;
      while (multiModeQuad_ROS_B.ifirst <= nm1) {
        alpha1[multiModeQuad_ROS_B.ifirst].re = (rtNaN);
        alpha1[multiModeQuad_ROS_B.ifirst].im = 0.0;
        beta1[multiModeQuad_ROS_B.ifirst].re = (rtNaN);
        beta1[multiModeQuad_ROS_B.ifirst].im = 0.0;
        multiModeQuad_ROS_B.ifirst++;
      }

      for (jp1 = 0; jp1 < 16; jp1++) {
        Z[jp1].re = (rtNaN);
        Z[jp1].im = 0.0;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (nm1 = 0; nm1 <= ilo - 2; nm1++) {
      alpha1[nm1] = A[(nm1 << 2) + nm1];
    }
  }
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_xztgevc(const creal_T A[16], creal_T V[16])
{
  int32_T i;
  multiModeQuad_ROS_B.rworka[0] = 0.0;
  multiModeQuad_ROS_B.rworka[2] = 0.0;
  multiModeQuad_ROS_B.rworka[3] = 0.0;
  multiModeQuad_ROS_B.anorm_l = fabs(A[0].re) + fabs(A[0].im);
  multiModeQuad_ROS_B.rworka[1] = fabs(A[4].re) + fabs(A[4].im);
  multiModeQuad_ROS_B.ascale_j = (fabs(A[5].re) + fabs(A[5].im)) +
    multiModeQuad_ROS_B.rworka[1];
  if (multiModeQuad_ROS_B.ascale_j > multiModeQuad_ROS_B.anorm_l) {
    multiModeQuad_ROS_B.anorm_l = multiModeQuad_ROS_B.ascale_j;
  }

  for (i = 0; i < 2; i++) {
    multiModeQuad_ROS_B.rworka[2] += fabs(A[i + 8].re) + fabs(A[i + 8].im);
  }

  multiModeQuad_ROS_B.ascale_j = (fabs(A[10].re) + fabs(A[10].im)) +
    multiModeQuad_ROS_B.rworka[2];
  if (multiModeQuad_ROS_B.ascale_j > multiModeQuad_ROS_B.anorm_l) {
    multiModeQuad_ROS_B.anorm_l = multiModeQuad_ROS_B.ascale_j;
  }

  for (i = 0; i < 3; i++) {
    multiModeQuad_ROS_B.rworka[3] += fabs(A[i + 12].re) + fabs(A[i + 12].im);
  }

  multiModeQuad_ROS_B.ascale_j = (fabs(A[15].re) + fabs(A[15].im)) +
    multiModeQuad_ROS_B.rworka[3];
  if (multiModeQuad_ROS_B.ascale_j > multiModeQuad_ROS_B.anorm_l) {
    multiModeQuad_ROS_B.anorm_l = multiModeQuad_ROS_B.ascale_j;
  }

  multiModeQuad_ROS_B.ascale_j = multiModeQuad_ROS_B.anorm_l;
  if (2.2250738585072014E-308 > multiModeQuad_ROS_B.anorm_l) {
    multiModeQuad_ROS_B.ascale_j = 2.2250738585072014E-308;
  }

  multiModeQuad_ROS_B.ascale_j = 1.0 / multiModeQuad_ROS_B.ascale_j;
  for (i = 3; i >= 0; i--) {
    real_T work2_idx_3_im;
    int32_T b_x_tmp;
    int32_T b_x_tmp_tmp;
    int32_T d_re_tmp;
    int32_T work2_idx_0_re_tmp;
    boolean_T lscalea;
    boolean_T lscaleb;
    b_x_tmp_tmp = i << 2;
    b_x_tmp = b_x_tmp_tmp + i;
    multiModeQuad_ROS_B.salpha_re = A[b_x_tmp].re;
    multiModeQuad_ROS_B.salpha_im = A[b_x_tmp].im;
    multiModeQuad_ROS_B.temp = (fabs(multiModeQuad_ROS_B.salpha_re) + fabs
      (multiModeQuad_ROS_B.salpha_im)) * multiModeQuad_ROS_B.ascale_j;
    if (1.0 > multiModeQuad_ROS_B.temp) {
      multiModeQuad_ROS_B.temp = 1.0;
    }

    multiModeQuad_ROS_B.temp = 1.0 / multiModeQuad_ROS_B.temp;
    multiModeQuad_ROS_B.salpha_re = multiModeQuad_ROS_B.salpha_re *
      multiModeQuad_ROS_B.temp * multiModeQuad_ROS_B.ascale_j;
    multiModeQuad_ROS_B.salpha_im = multiModeQuad_ROS_B.salpha_im *
      multiModeQuad_ROS_B.temp * multiModeQuad_ROS_B.ascale_j;
    multiModeQuad_ROS_B.acoeff = multiModeQuad_ROS_B.temp *
      multiModeQuad_ROS_B.ascale_j;
    lscalea = ((multiModeQuad_ROS_B.temp >= 2.2250738585072014E-308) &&
               (multiModeQuad_ROS_B.acoeff < 4.0083367200179456E-292));
    multiModeQuad_ROS_B.dmin = fabs(multiModeQuad_ROS_B.salpha_re) + fabs
      (multiModeQuad_ROS_B.salpha_im);
    lscaleb = ((multiModeQuad_ROS_B.dmin >= 2.2250738585072014E-308) &&
               (multiModeQuad_ROS_B.dmin < 4.0083367200179456E-292));
    multiModeQuad_ROS_B.scale_d = 1.0;
    if (lscalea) {
      multiModeQuad_ROS_B.scale_d = multiModeQuad_ROS_B.anorm_l;
      if (2.4948003869184E+291 < multiModeQuad_ROS_B.anorm_l) {
        multiModeQuad_ROS_B.scale_d = 2.4948003869184E+291;
      }

      multiModeQuad_ROS_B.scale_d *= 4.0083367200179456E-292 /
        multiModeQuad_ROS_B.temp;
    }

    if (lscaleb) {
      multiModeQuad_ROS_B.work2_idx_2_im = 4.0083367200179456E-292 /
        multiModeQuad_ROS_B.dmin;
      if (multiModeQuad_ROS_B.work2_idx_2_im > multiModeQuad_ROS_B.scale_d) {
        multiModeQuad_ROS_B.scale_d = multiModeQuad_ROS_B.work2_idx_2_im;
      }
    }

    if (lscalea || lscaleb) {
      multiModeQuad_ROS_B.work2_idx_2_im = multiModeQuad_ROS_B.acoeff;
      if (1.0 > multiModeQuad_ROS_B.acoeff) {
        multiModeQuad_ROS_B.work2_idx_2_im = 1.0;
      }

      if (multiModeQuad_ROS_B.dmin > multiModeQuad_ROS_B.work2_idx_2_im) {
        multiModeQuad_ROS_B.work2_idx_2_im = multiModeQuad_ROS_B.dmin;
      }

      multiModeQuad_ROS_B.dmin = 1.0 / (2.2250738585072014E-308 *
        multiModeQuad_ROS_B.work2_idx_2_im);
      if (multiModeQuad_ROS_B.dmin < multiModeQuad_ROS_B.scale_d) {
        multiModeQuad_ROS_B.scale_d = multiModeQuad_ROS_B.dmin;
      }

      if (lscalea) {
        multiModeQuad_ROS_B.acoeff = multiModeQuad_ROS_B.scale_d *
          multiModeQuad_ROS_B.temp * multiModeQuad_ROS_B.ascale_j;
      } else {
        multiModeQuad_ROS_B.acoeff *= multiModeQuad_ROS_B.scale_d;
      }

      multiModeQuad_ROS_B.salpha_re *= multiModeQuad_ROS_B.scale_d;
      multiModeQuad_ROS_B.salpha_im *= multiModeQuad_ROS_B.scale_d;
    }

    memset(&multiModeQuad_ROS_B.work1[0], 0, sizeof(creal_T) << 2U);
    multiModeQuad_ROS_B.work1[i].re = 1.0;
    multiModeQuad_ROS_B.work1[i].im = 0.0;
    multiModeQuad_ROS_B.dmin = 2.2204460492503131E-16 *
      multiModeQuad_ROS_B.acoeff * multiModeQuad_ROS_B.anorm_l;
    multiModeQuad_ROS_B.temp = (fabs(multiModeQuad_ROS_B.salpha_re) + fabs
      (multiModeQuad_ROS_B.salpha_im)) * 2.2204460492503131E-16;
    if (multiModeQuad_ROS_B.temp > multiModeQuad_ROS_B.dmin) {
      multiModeQuad_ROS_B.dmin = multiModeQuad_ROS_B.temp;
    }

    if (2.2250738585072014E-308 > multiModeQuad_ROS_B.dmin) {
      multiModeQuad_ROS_B.dmin = 2.2250738585072014E-308;
    }

    for (b_x_tmp = 0; b_x_tmp < i; b_x_tmp++) {
      d_re_tmp = b_x_tmp_tmp + b_x_tmp;
      multiModeQuad_ROS_B.work1[b_x_tmp].re = A[d_re_tmp].re *
        multiModeQuad_ROS_B.acoeff;
      multiModeQuad_ROS_B.work1[b_x_tmp].im = A[d_re_tmp].im *
        multiModeQuad_ROS_B.acoeff;
    }

    multiModeQuad_ROS_B.work1[i].re = 1.0;
    multiModeQuad_ROS_B.work1[i].im = 0.0;
    for (b_x_tmp = i - 1; b_x_tmp + 1 > 0; b_x_tmp--) {
      work2_idx_0_re_tmp = b_x_tmp << 2;
      d_re_tmp = work2_idx_0_re_tmp + b_x_tmp;
      multiModeQuad_ROS_B.work2_idx_3_re = A[d_re_tmp].re *
        multiModeQuad_ROS_B.acoeff - multiModeQuad_ROS_B.salpha_re;
      multiModeQuad_ROS_B.scale_d = A[d_re_tmp].im * multiModeQuad_ROS_B.acoeff
        - multiModeQuad_ROS_B.salpha_im;
      if (fabs(multiModeQuad_ROS_B.work2_idx_3_re) + fabs
          (multiModeQuad_ROS_B.scale_d) <= multiModeQuad_ROS_B.dmin) {
        multiModeQuad_ROS_B.work2_idx_3_re = multiModeQuad_ROS_B.dmin;
        multiModeQuad_ROS_B.scale_d = 0.0;
      }

      multiModeQuad_ROS_B.work2_idx_2_im = fabs
        (multiModeQuad_ROS_B.work2_idx_3_re);
      multiModeQuad_ROS_B.e_y = fabs(multiModeQuad_ROS_B.scale_d);
      multiModeQuad_ROS_B.temp = multiModeQuad_ROS_B.work2_idx_2_im +
        multiModeQuad_ROS_B.e_y;
      if (multiModeQuad_ROS_B.temp < 1.0) {
        work2_idx_3_im = fabs(multiModeQuad_ROS_B.work1[b_x_tmp].re) + fabs
          (multiModeQuad_ROS_B.work1[b_x_tmp].im);
        if (work2_idx_3_im >= multiModeQuad_ROS_B.temp * 1.1235582092889474E+307)
        {
          multiModeQuad_ROS_B.temp = 1.0 / work2_idx_3_im;
          for (d_re_tmp = 0; d_re_tmp <= i; d_re_tmp++) {
            multiModeQuad_ROS_B.work1[d_re_tmp].re *= multiModeQuad_ROS_B.temp;
            multiModeQuad_ROS_B.work1[d_re_tmp].im *= multiModeQuad_ROS_B.temp;
          }
        }
      }

      multiModeQuad_ROS_B.temp = multiModeQuad_ROS_B.work1[b_x_tmp].re;
      work2_idx_3_im = multiModeQuad_ROS_B.work1[b_x_tmp].im;
      if (multiModeQuad_ROS_B.scale_d == 0.0) {
        if (-work2_idx_3_im == 0.0) {
          multiModeQuad_ROS_B.work1[b_x_tmp].re = -multiModeQuad_ROS_B.temp /
            multiModeQuad_ROS_B.work2_idx_3_re;
          multiModeQuad_ROS_B.work1[b_x_tmp].im = 0.0;
        } else if (-multiModeQuad_ROS_B.temp == 0.0) {
          multiModeQuad_ROS_B.work1[b_x_tmp].re = 0.0;
          multiModeQuad_ROS_B.work1[b_x_tmp].im = -work2_idx_3_im /
            multiModeQuad_ROS_B.work2_idx_3_re;
        } else {
          multiModeQuad_ROS_B.work1[b_x_tmp].re = -multiModeQuad_ROS_B.temp /
            multiModeQuad_ROS_B.work2_idx_3_re;
          multiModeQuad_ROS_B.work1[b_x_tmp].im = -work2_idx_3_im /
            multiModeQuad_ROS_B.work2_idx_3_re;
        }
      } else if (multiModeQuad_ROS_B.work2_idx_3_re == 0.0) {
        if (-multiModeQuad_ROS_B.temp == 0.0) {
          multiModeQuad_ROS_B.work1[b_x_tmp].re = -work2_idx_3_im /
            multiModeQuad_ROS_B.scale_d;
          multiModeQuad_ROS_B.work1[b_x_tmp].im = 0.0;
        } else if (-work2_idx_3_im == 0.0) {
          multiModeQuad_ROS_B.work1[b_x_tmp].re = 0.0;
          multiModeQuad_ROS_B.work1[b_x_tmp].im = -(-multiModeQuad_ROS_B.temp /
            multiModeQuad_ROS_B.scale_d);
        } else {
          multiModeQuad_ROS_B.work1[b_x_tmp].re = -work2_idx_3_im /
            multiModeQuad_ROS_B.scale_d;
          multiModeQuad_ROS_B.work1[b_x_tmp].im = -(-multiModeQuad_ROS_B.temp /
            multiModeQuad_ROS_B.scale_d);
        }
      } else if (multiModeQuad_ROS_B.work2_idx_2_im > multiModeQuad_ROS_B.e_y) {
        multiModeQuad_ROS_B.work2_idx_2_im = multiModeQuad_ROS_B.scale_d /
          multiModeQuad_ROS_B.work2_idx_3_re;
        multiModeQuad_ROS_B.scale_d = multiModeQuad_ROS_B.work2_idx_2_im *
          multiModeQuad_ROS_B.scale_d + multiModeQuad_ROS_B.work2_idx_3_re;
        multiModeQuad_ROS_B.work1[b_x_tmp].re =
          (multiModeQuad_ROS_B.work2_idx_2_im * -work2_idx_3_im +
           -multiModeQuad_ROS_B.temp) / multiModeQuad_ROS_B.scale_d;
        multiModeQuad_ROS_B.work1[b_x_tmp].im = (-work2_idx_3_im -
          multiModeQuad_ROS_B.work2_idx_2_im * -multiModeQuad_ROS_B.temp) /
          multiModeQuad_ROS_B.scale_d;
      } else if (multiModeQuad_ROS_B.e_y == multiModeQuad_ROS_B.work2_idx_2_im)
      {
        multiModeQuad_ROS_B.work2_idx_3_re = multiModeQuad_ROS_B.work2_idx_3_re >
          0.0 ? 0.5 : -0.5;
        multiModeQuad_ROS_B.scale_d = multiModeQuad_ROS_B.scale_d > 0.0 ? 0.5 :
          -0.5;
        multiModeQuad_ROS_B.work1[b_x_tmp].re = (-multiModeQuad_ROS_B.temp *
          multiModeQuad_ROS_B.work2_idx_3_re + -work2_idx_3_im *
          multiModeQuad_ROS_B.scale_d) / multiModeQuad_ROS_B.work2_idx_2_im;
        multiModeQuad_ROS_B.work1[b_x_tmp].im = (-work2_idx_3_im *
          multiModeQuad_ROS_B.work2_idx_3_re - -multiModeQuad_ROS_B.temp *
          multiModeQuad_ROS_B.scale_d) / multiModeQuad_ROS_B.work2_idx_2_im;
      } else {
        multiModeQuad_ROS_B.work2_idx_2_im = multiModeQuad_ROS_B.work2_idx_3_re /
          multiModeQuad_ROS_B.scale_d;
        multiModeQuad_ROS_B.scale_d += multiModeQuad_ROS_B.work2_idx_2_im *
          multiModeQuad_ROS_B.work2_idx_3_re;
        multiModeQuad_ROS_B.work1[b_x_tmp].re =
          (multiModeQuad_ROS_B.work2_idx_2_im * -multiModeQuad_ROS_B.temp +
           -work2_idx_3_im) / multiModeQuad_ROS_B.scale_d;
        multiModeQuad_ROS_B.work1[b_x_tmp].im =
          (multiModeQuad_ROS_B.work2_idx_2_im * -work2_idx_3_im -
           (-multiModeQuad_ROS_B.temp)) / multiModeQuad_ROS_B.scale_d;
      }

      if (b_x_tmp + 1 > 1) {
        if (fabs(multiModeQuad_ROS_B.work1[b_x_tmp].re) + fabs
            (multiModeQuad_ROS_B.work1[b_x_tmp].im) > 1.0) {
          multiModeQuad_ROS_B.temp = 1.0 / (fabs
            (multiModeQuad_ROS_B.work1[b_x_tmp].re) + fabs
            (multiModeQuad_ROS_B.work1[b_x_tmp].im));
          if (multiModeQuad_ROS_B.acoeff * multiModeQuad_ROS_B.rworka[b_x_tmp] >=
              1.1235582092889474E+307 * multiModeQuad_ROS_B.temp) {
            for (d_re_tmp = 0; d_re_tmp <= i; d_re_tmp++) {
              multiModeQuad_ROS_B.work1[d_re_tmp].re *= multiModeQuad_ROS_B.temp;
              multiModeQuad_ROS_B.work1[d_re_tmp].im *= multiModeQuad_ROS_B.temp;
            }
          }
        }

        multiModeQuad_ROS_B.work2_idx_3_re = multiModeQuad_ROS_B.acoeff *
          multiModeQuad_ROS_B.work1[b_x_tmp].re;
        multiModeQuad_ROS_B.scale_d = multiModeQuad_ROS_B.acoeff *
          multiModeQuad_ROS_B.work1[b_x_tmp].im;
        for (int32_T e_jr = 0; e_jr < b_x_tmp; e_jr++) {
          d_re_tmp = work2_idx_0_re_tmp + e_jr;
          multiModeQuad_ROS_B.temp = A[d_re_tmp].im;
          work2_idx_3_im = A[d_re_tmp].re;
          multiModeQuad_ROS_B.work1[e_jr].re += work2_idx_3_im *
            multiModeQuad_ROS_B.work2_idx_3_re - multiModeQuad_ROS_B.temp *
            multiModeQuad_ROS_B.scale_d;
          multiModeQuad_ROS_B.work1[e_jr].im += multiModeQuad_ROS_B.temp *
            multiModeQuad_ROS_B.work2_idx_3_re + work2_idx_3_im *
            multiModeQuad_ROS_B.scale_d;
        }
      }
    }

    multiModeQuad_ROS_B.salpha_re = 0.0;
    multiModeQuad_ROS_B.salpha_im = 0.0;
    multiModeQuad_ROS_B.acoeff = 0.0;
    multiModeQuad_ROS_B.dmin = 0.0;
    multiModeQuad_ROS_B.scale_d = 0.0;
    multiModeQuad_ROS_B.work2_idx_2_im = 0.0;
    multiModeQuad_ROS_B.work2_idx_3_re = 0.0;
    work2_idx_3_im = 0.0;
    for (b_x_tmp = 0; b_x_tmp <= i; b_x_tmp++) {
      real_T work2_idx_0_re_tmp_0;
      real_T work2_idx_0_re_tmp_1;
      work2_idx_0_re_tmp = b_x_tmp << 2;
      work2_idx_0_re_tmp_0 = V[work2_idx_0_re_tmp].re;
      multiModeQuad_ROS_B.temp = multiModeQuad_ROS_B.work1[b_x_tmp].im;
      work2_idx_0_re_tmp_1 = V[work2_idx_0_re_tmp].im;
      multiModeQuad_ROS_B.e_y = multiModeQuad_ROS_B.work1[b_x_tmp].re;
      multiModeQuad_ROS_B.salpha_re += work2_idx_0_re_tmp_0 *
        multiModeQuad_ROS_B.e_y - work2_idx_0_re_tmp_1 *
        multiModeQuad_ROS_B.temp;
      multiModeQuad_ROS_B.salpha_im += work2_idx_0_re_tmp_0 *
        multiModeQuad_ROS_B.temp + work2_idx_0_re_tmp_1 *
        multiModeQuad_ROS_B.e_y;
      work2_idx_0_re_tmp_0 = V[work2_idx_0_re_tmp + 1].re;
      work2_idx_0_re_tmp_1 = V[work2_idx_0_re_tmp + 1].im;
      multiModeQuad_ROS_B.acoeff += work2_idx_0_re_tmp_0 *
        multiModeQuad_ROS_B.e_y - work2_idx_0_re_tmp_1 *
        multiModeQuad_ROS_B.temp;
      multiModeQuad_ROS_B.dmin += work2_idx_0_re_tmp_0 *
        multiModeQuad_ROS_B.temp + work2_idx_0_re_tmp_1 *
        multiModeQuad_ROS_B.e_y;
      work2_idx_0_re_tmp_0 = V[work2_idx_0_re_tmp + 2].re;
      work2_idx_0_re_tmp_1 = V[work2_idx_0_re_tmp + 2].im;
      multiModeQuad_ROS_B.scale_d += work2_idx_0_re_tmp_0 *
        multiModeQuad_ROS_B.e_y - work2_idx_0_re_tmp_1 *
        multiModeQuad_ROS_B.temp;
      multiModeQuad_ROS_B.work2_idx_2_im += work2_idx_0_re_tmp_0 *
        multiModeQuad_ROS_B.temp + work2_idx_0_re_tmp_1 *
        multiModeQuad_ROS_B.e_y;
      work2_idx_0_re_tmp_0 = V[work2_idx_0_re_tmp + 3].re;
      work2_idx_0_re_tmp_1 = V[work2_idx_0_re_tmp + 3].im;
      multiModeQuad_ROS_B.work2_idx_3_re += work2_idx_0_re_tmp_0 *
        multiModeQuad_ROS_B.e_y - work2_idx_0_re_tmp_1 *
        multiModeQuad_ROS_B.temp;
      work2_idx_3_im += work2_idx_0_re_tmp_0 * multiModeQuad_ROS_B.temp +
        work2_idx_0_re_tmp_1 * multiModeQuad_ROS_B.e_y;
    }

    multiModeQuad_ROS_B.temp = fabs(multiModeQuad_ROS_B.salpha_re) + fabs
      (multiModeQuad_ROS_B.salpha_im);
    multiModeQuad_ROS_B.e_y = fabs(multiModeQuad_ROS_B.acoeff) + fabs
      (multiModeQuad_ROS_B.dmin);
    if (multiModeQuad_ROS_B.e_y > multiModeQuad_ROS_B.temp) {
      multiModeQuad_ROS_B.temp = multiModeQuad_ROS_B.e_y;
    }

    multiModeQuad_ROS_B.e_y = fabs(multiModeQuad_ROS_B.scale_d) + fabs
      (multiModeQuad_ROS_B.work2_idx_2_im);
    if (multiModeQuad_ROS_B.e_y > multiModeQuad_ROS_B.temp) {
      multiModeQuad_ROS_B.temp = multiModeQuad_ROS_B.e_y;
    }

    multiModeQuad_ROS_B.e_y = fabs(multiModeQuad_ROS_B.work2_idx_3_re) + fabs
      (work2_idx_3_im);
    if (multiModeQuad_ROS_B.e_y > multiModeQuad_ROS_B.temp) {
      multiModeQuad_ROS_B.temp = multiModeQuad_ROS_B.e_y;
    }

    if (multiModeQuad_ROS_B.temp > 2.2250738585072014E-308) {
      multiModeQuad_ROS_B.temp = 1.0 / multiModeQuad_ROS_B.temp;
      V[b_x_tmp_tmp].re = multiModeQuad_ROS_B.temp *
        multiModeQuad_ROS_B.salpha_re;
      V[b_x_tmp_tmp].im = multiModeQuad_ROS_B.temp *
        multiModeQuad_ROS_B.salpha_im;
      d_re_tmp = (i << 2) + 1;
      V[d_re_tmp].re = multiModeQuad_ROS_B.temp * multiModeQuad_ROS_B.acoeff;
      V[d_re_tmp].im = multiModeQuad_ROS_B.temp * multiModeQuad_ROS_B.dmin;
      d_re_tmp = (i << 2) + 2;
      V[d_re_tmp].re = multiModeQuad_ROS_B.temp * multiModeQuad_ROS_B.scale_d;
      V[d_re_tmp].im = multiModeQuad_ROS_B.temp *
        multiModeQuad_ROS_B.work2_idx_2_im;
      d_re_tmp = (i << 2) + 3;
      V[d_re_tmp].re = multiModeQuad_ROS_B.temp *
        multiModeQuad_ROS_B.work2_idx_3_re;
      V[d_re_tmp].im = multiModeQuad_ROS_B.temp * work2_idx_3_im;
    } else {
      V[b_x_tmp_tmp].re = 0.0;
      V[b_x_tmp_tmp].im = 0.0;
      V[b_x_tmp_tmp + 1].re = 0.0;
      V[b_x_tmp_tmp + 1].im = 0.0;
      V[b_x_tmp_tmp + 2].re = 0.0;
      V[b_x_tmp_tmp + 2].im = 0.0;
      V[b_x_tmp_tmp + 3].re = 0.0;
      V[b_x_tmp_tmp + 3].im = 0.0;
    }
  }
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_eigStandard(const real_T A[16], creal_T V[16],
  creal_T D[4])
{
  boolean_T exitg1;
  for (multiModeQuad_ROS_B.stemp_re_tmp = 0; multiModeQuad_ROS_B.stemp_re_tmp <
       16; multiModeQuad_ROS_B.stemp_re_tmp++) {
    multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].re =
      A[multiModeQuad_ROS_B.stemp_re_tmp];
    multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].im = 0.0;
  }

  multiModeQuad_ROS_B.anrm = 0.0;
  multiModeQuad_ROS_B.b_k = 0;
  exitg1 = false;
  while ((!exitg1) && (multiModeQuad_ROS_B.b_k < 16)) {
    multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_rt_hypotd_snf
      (multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.b_k].re,
       multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.b_k].im);
    if (rtIsNaN(multiModeQuad_ROS_B.b_absxk)) {
      multiModeQuad_ROS_B.anrm = (rtNaN);
      exitg1 = true;
    } else {
      if (multiModeQuad_ROS_B.b_absxk > multiModeQuad_ROS_B.anrm) {
        multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.b_absxk;
      }

      multiModeQuad_ROS_B.b_k++;
    }
  }

  if (rtIsInf(multiModeQuad_ROS_B.anrm) || rtIsNaN(multiModeQuad_ROS_B.anrm)) {
    D[0].re = (rtNaN);
    D[0].im = 0.0;
    multiModeQuad_ROS_B.y[0].re = (rtNaN);
    multiModeQuad_ROS_B.y[0].im = 0.0;
    D[1].re = (rtNaN);
    D[1].im = 0.0;
    multiModeQuad_ROS_B.y[1].re = (rtNaN);
    multiModeQuad_ROS_B.y[1].im = 0.0;
    D[2].re = (rtNaN);
    D[2].im = 0.0;
    multiModeQuad_ROS_B.y[2].re = (rtNaN);
    multiModeQuad_ROS_B.y[2].im = 0.0;
    D[3].re = (rtNaN);
    D[3].im = 0.0;
    multiModeQuad_ROS_B.y[3].re = (rtNaN);
    multiModeQuad_ROS_B.y[3].im = 0.0;
    for (multiModeQuad_ROS_B.stemp_re_tmp = 0; multiModeQuad_ROS_B.stemp_re_tmp <
         16; multiModeQuad_ROS_B.stemp_re_tmp++) {
      V[multiModeQuad_ROS_B.stemp_re_tmp].re = (rtNaN);
      V[multiModeQuad_ROS_B.stemp_re_tmp].im = 0.0;
    }
  } else {
    boolean_T guard1 = false;
    boolean_T ilascl;
    ilascl = false;
    multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.anrm;
    guard1 = false;
    if ((multiModeQuad_ROS_B.anrm > 0.0) && (multiModeQuad_ROS_B.anrm <
         6.7178761075670888E-139)) {
      multiModeQuad_ROS_B.b_absxk = 6.7178761075670888E-139;
      ilascl = true;
      guard1 = true;
    } else if (multiModeQuad_ROS_B.anrm > 1.4885657073574029E+138) {
      multiModeQuad_ROS_B.b_absxk = 1.4885657073574029E+138;
      ilascl = true;
      guard1 = true;
    }

    if (guard1) {
      boolean_T notdone;
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.anrm;
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.b_absxk;
      notdone = true;
      while (notdone) {
        multiModeQuad_ROS_B.vtemp = multiModeQuad_ROS_B.cfromc *
          2.0041683600089728E-292;
        multiModeQuad_ROS_B.cto1 = multiModeQuad_ROS_B.ctoc /
          4.9896007738368E+291;
        if ((multiModeQuad_ROS_B.vtemp > multiModeQuad_ROS_B.ctoc) &&
            (multiModeQuad_ROS_B.ctoc != 0.0)) {
          multiModeQuad_ROS_B.stemp_im_tmp = 2.0041683600089728E-292;
          multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.vtemp;
        } else if (multiModeQuad_ROS_B.cto1 > multiModeQuad_ROS_B.cfromc) {
          multiModeQuad_ROS_B.stemp_im_tmp = 4.9896007738368E+291;
          multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.cto1;
        } else {
          multiModeQuad_ROS_B.stemp_im_tmp = multiModeQuad_ROS_B.ctoc /
            multiModeQuad_ROS_B.cfromc;
          notdone = false;
        }

        for (multiModeQuad_ROS_B.stemp_re_tmp = 0;
             multiModeQuad_ROS_B.stemp_re_tmp < 16;
             multiModeQuad_ROS_B.stemp_re_tmp++) {
          multiModeQuad_ROS_B.tmp =
            multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp];
          multiModeQuad_ROS_B.tmp.re *= multiModeQuad_ROS_B.stemp_im_tmp;
          multiModeQuad_ROS_B.tmp.im *= multiModeQuad_ROS_B.stemp_im_tmp;
          multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp] =
            multiModeQuad_ROS_B.tmp;
        }
      }
    }

    multiModeQuad_ROS_xzggbal(multiModeQuad_ROS_B.At, &multiModeQuad_ROS_B.c_i_e,
      &multiModeQuad_ROS_B.b_k, multiModeQuad_ROS_B.rscale);
    memset(&V[0], 0, sizeof(creal_T) << 4U);
    V[0].re = 1.0;
    V[0].im = 0.0;
    V[5].re = 1.0;
    V[5].im = 0.0;
    V[10].re = 1.0;
    V[10].im = 0.0;
    V[15].re = 1.0;
    V[15].im = 0.0;
    if (multiModeQuad_ROS_B.b_k >= multiModeQuad_ROS_B.c_i_e + 2) {
      multiModeQuad_ROS_B.jcol = multiModeQuad_ROS_B.c_i_e - 1;
      while (multiModeQuad_ROS_B.jcol + 1 < multiModeQuad_ROS_B.b_k - 1) {
        multiModeQuad_ROS_B.jrow = multiModeQuad_ROS_B.b_k - 1;
        while (multiModeQuad_ROS_B.jrow + 1 > multiModeQuad_ROS_B.jcol + 2) {
          multiModeQuad_ROS_xzlartg(multiModeQuad_ROS_B.At
            [(multiModeQuad_ROS_B.jrow + (multiModeQuad_ROS_B.jcol << 2)) - 1],
            multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.jrow +
            (multiModeQuad_ROS_B.jcol << 2)], &multiModeQuad_ROS_B.cfromc,
            &multiModeQuad_ROS_B.tmp, &multiModeQuad_ROS_B.At
            [(multiModeQuad_ROS_B.jrow + (multiModeQuad_ROS_B.jcol << 2)) - 1]);
          multiModeQuad_ROS_B.stemp_re_tmp = (multiModeQuad_ROS_B.jcol << 2) +
            multiModeQuad_ROS_B.jrow;
          multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].re = 0.0;
          multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].im = 0.0;
          multiModeQuad_ROS_B.b_j = multiModeQuad_ROS_B.jcol + 1;
          while (multiModeQuad_ROS_B.b_j + 1 < 5) {
            multiModeQuad_ROS_B.stemp_re_tmp = (multiModeQuad_ROS_B.b_j << 2) +
              multiModeQuad_ROS_B.jrow;
            multiModeQuad_ROS_B.ctoc =
              multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].im;
            multiModeQuad_ROS_B.vtemp =
              multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].re;
            multiModeQuad_ROS_B.cto1 =
              multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp - 1].re;
            multiModeQuad_ROS_B.stemp_im_tmp =
              multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp - 1].im;
            multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].re =
              multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.cfromc -
              (multiModeQuad_ROS_B.cto1 * multiModeQuad_ROS_B.tmp.re +
               multiModeQuad_ROS_B.stemp_im_tmp * multiModeQuad_ROS_B.tmp.im);
            multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].im =
              multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.cfromc -
              (multiModeQuad_ROS_B.stemp_im_tmp * multiModeQuad_ROS_B.tmp.re -
               multiModeQuad_ROS_B.tmp.im * multiModeQuad_ROS_B.cto1);
            multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp - 1].re =
              (multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.tmp.re -
               multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.tmp.im) +
              multiModeQuad_ROS_B.cto1 * multiModeQuad_ROS_B.cfromc;
            multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp - 1].im =
              (multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.tmp.re +
               multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.tmp.im) +
              multiModeQuad_ROS_B.stemp_im_tmp * multiModeQuad_ROS_B.cfromc;
            multiModeQuad_ROS_B.b_j++;
          }

          multiModeQuad_ROS_B.tmp.re = -multiModeQuad_ROS_B.tmp.re;
          multiModeQuad_ROS_B.tmp.im = -multiModeQuad_ROS_B.tmp.im;
          multiModeQuad_ROS_B.b_j = 0;
          while (multiModeQuad_ROS_B.b_j + 1 <= multiModeQuad_ROS_B.b_k) {
            multiModeQuad_ROS_B.stemp_re_tmp = ((multiModeQuad_ROS_B.jrow - 1) <<
              2) + multiModeQuad_ROS_B.b_j;
            multiModeQuad_ROS_B.ctoc =
              multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].im;
            multiModeQuad_ROS_B.vtemp =
              multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].re;
            multiModeQuad_ROS_B.stemp_re_tmp_b = (multiModeQuad_ROS_B.jrow << 2)
              + multiModeQuad_ROS_B.b_j;
            multiModeQuad_ROS_B.cto1 =
              multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp_b].re;
            multiModeQuad_ROS_B.stemp_im_tmp =
              multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp_b].im;
            multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].re =
              multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.cfromc -
              (multiModeQuad_ROS_B.cto1 * multiModeQuad_ROS_B.tmp.re +
               multiModeQuad_ROS_B.stemp_im_tmp * multiModeQuad_ROS_B.tmp.im);
            multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp].im =
              multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.cfromc -
              (multiModeQuad_ROS_B.stemp_im_tmp * multiModeQuad_ROS_B.tmp.re -
               multiModeQuad_ROS_B.tmp.im * multiModeQuad_ROS_B.cto1);
            multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp_b].re =
              (multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.tmp.re -
               multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.tmp.im) +
              multiModeQuad_ROS_B.cto1 * multiModeQuad_ROS_B.cfromc;
            multiModeQuad_ROS_B.At[multiModeQuad_ROS_B.stemp_re_tmp_b].im =
              (multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.tmp.re +
               multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.tmp.im) +
              multiModeQuad_ROS_B.stemp_im_tmp * multiModeQuad_ROS_B.cfromc;
            multiModeQuad_ROS_B.b_j++;
          }

          multiModeQuad_ROS_B.stemp_re_tmp = (multiModeQuad_ROS_B.jrow - 1) << 2;
          multiModeQuad_ROS_B.ctoc = V[multiModeQuad_ROS_B.stemp_re_tmp].im;
          multiModeQuad_ROS_B.vtemp = V[multiModeQuad_ROS_B.stemp_re_tmp].re;
          multiModeQuad_ROS_B.stemp_re_tmp_b = multiModeQuad_ROS_B.jrow << 2;
          multiModeQuad_ROS_B.cto1 = V[multiModeQuad_ROS_B.stemp_re_tmp_b].re;
          multiModeQuad_ROS_B.stemp_im_tmp =
            V[multiModeQuad_ROS_B.stemp_re_tmp_b].im;
          V[multiModeQuad_ROS_B.stemp_re_tmp].re = multiModeQuad_ROS_B.vtemp *
            multiModeQuad_ROS_B.cfromc - (multiModeQuad_ROS_B.cto1 *
            multiModeQuad_ROS_B.tmp.re + multiModeQuad_ROS_B.stemp_im_tmp *
            multiModeQuad_ROS_B.tmp.im);
          V[multiModeQuad_ROS_B.stemp_re_tmp].im = multiModeQuad_ROS_B.ctoc *
            multiModeQuad_ROS_B.cfromc - (multiModeQuad_ROS_B.stemp_im_tmp *
            multiModeQuad_ROS_B.tmp.re - multiModeQuad_ROS_B.tmp.im *
            multiModeQuad_ROS_B.cto1);
          V[multiModeQuad_ROS_B.stemp_re_tmp_b].re = (multiModeQuad_ROS_B.vtemp *
            multiModeQuad_ROS_B.tmp.re - multiModeQuad_ROS_B.ctoc *
            multiModeQuad_ROS_B.tmp.im) + multiModeQuad_ROS_B.cto1 *
            multiModeQuad_ROS_B.cfromc;
          V[multiModeQuad_ROS_B.stemp_re_tmp_b].im = (multiModeQuad_ROS_B.ctoc *
            multiModeQuad_ROS_B.tmp.re + multiModeQuad_ROS_B.vtemp *
            multiModeQuad_ROS_B.tmp.im) + multiModeQuad_ROS_B.stemp_im_tmp *
            multiModeQuad_ROS_B.cfromc;
          multiModeQuad_ROS_B.ctoc = V[multiModeQuad_ROS_B.stemp_re_tmp + 1].im;
          multiModeQuad_ROS_B.vtemp = V[multiModeQuad_ROS_B.stemp_re_tmp + 1].re;
          multiModeQuad_ROS_B.cto1 = V[multiModeQuad_ROS_B.stemp_re_tmp_b + 1].
            re;
          multiModeQuad_ROS_B.stemp_im_tmp =
            V[multiModeQuad_ROS_B.stemp_re_tmp_b + 1].im;
          V[multiModeQuad_ROS_B.stemp_re_tmp + 1].re = multiModeQuad_ROS_B.vtemp
            * multiModeQuad_ROS_B.cfromc - (multiModeQuad_ROS_B.cto1 *
            multiModeQuad_ROS_B.tmp.re + multiModeQuad_ROS_B.stemp_im_tmp *
            multiModeQuad_ROS_B.tmp.im);
          V[multiModeQuad_ROS_B.stemp_re_tmp + 1].im = multiModeQuad_ROS_B.ctoc *
            multiModeQuad_ROS_B.cfromc - (multiModeQuad_ROS_B.stemp_im_tmp *
            multiModeQuad_ROS_B.tmp.re - multiModeQuad_ROS_B.tmp.im *
            multiModeQuad_ROS_B.cto1);
          V[multiModeQuad_ROS_B.stemp_re_tmp_b + 1].re =
            (multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.tmp.re -
             multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.tmp.im) +
            multiModeQuad_ROS_B.cto1 * multiModeQuad_ROS_B.cfromc;
          V[multiModeQuad_ROS_B.stemp_re_tmp_b + 1].im =
            (multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.tmp.re +
             multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.tmp.im) +
            multiModeQuad_ROS_B.stemp_im_tmp * multiModeQuad_ROS_B.cfromc;
          multiModeQuad_ROS_B.ctoc = V[multiModeQuad_ROS_B.stemp_re_tmp + 2].im;
          multiModeQuad_ROS_B.vtemp = V[multiModeQuad_ROS_B.stemp_re_tmp + 2].re;
          multiModeQuad_ROS_B.cto1 = V[multiModeQuad_ROS_B.stemp_re_tmp_b + 2].
            re;
          multiModeQuad_ROS_B.stemp_im_tmp =
            V[multiModeQuad_ROS_B.stemp_re_tmp_b + 2].im;
          V[multiModeQuad_ROS_B.stemp_re_tmp + 2].re = multiModeQuad_ROS_B.vtemp
            * multiModeQuad_ROS_B.cfromc - (multiModeQuad_ROS_B.cto1 *
            multiModeQuad_ROS_B.tmp.re + multiModeQuad_ROS_B.stemp_im_tmp *
            multiModeQuad_ROS_B.tmp.im);
          V[multiModeQuad_ROS_B.stemp_re_tmp + 2].im = multiModeQuad_ROS_B.ctoc *
            multiModeQuad_ROS_B.cfromc - (multiModeQuad_ROS_B.stemp_im_tmp *
            multiModeQuad_ROS_B.tmp.re - multiModeQuad_ROS_B.tmp.im *
            multiModeQuad_ROS_B.cto1);
          V[multiModeQuad_ROS_B.stemp_re_tmp_b + 2].re =
            (multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.tmp.re -
             multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.tmp.im) +
            multiModeQuad_ROS_B.cto1 * multiModeQuad_ROS_B.cfromc;
          V[multiModeQuad_ROS_B.stemp_re_tmp_b + 2].im =
            (multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.tmp.re +
             multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.tmp.im) +
            multiModeQuad_ROS_B.stemp_im_tmp * multiModeQuad_ROS_B.cfromc;
          multiModeQuad_ROS_B.ctoc = V[multiModeQuad_ROS_B.stemp_re_tmp + 3].im;
          multiModeQuad_ROS_B.vtemp = V[multiModeQuad_ROS_B.stemp_re_tmp + 3].re;
          multiModeQuad_ROS_B.cto1 = V[multiModeQuad_ROS_B.stemp_re_tmp_b + 3].
            re;
          multiModeQuad_ROS_B.stemp_im_tmp =
            V[multiModeQuad_ROS_B.stemp_re_tmp_b + 3].im;
          V[multiModeQuad_ROS_B.stemp_re_tmp + 3].re = multiModeQuad_ROS_B.vtemp
            * multiModeQuad_ROS_B.cfromc - (multiModeQuad_ROS_B.cto1 *
            multiModeQuad_ROS_B.tmp.re + multiModeQuad_ROS_B.stemp_im_tmp *
            multiModeQuad_ROS_B.tmp.im);
          V[multiModeQuad_ROS_B.stemp_re_tmp + 3].im = multiModeQuad_ROS_B.ctoc *
            multiModeQuad_ROS_B.cfromc - (multiModeQuad_ROS_B.stemp_im_tmp *
            multiModeQuad_ROS_B.tmp.re - multiModeQuad_ROS_B.tmp.im *
            multiModeQuad_ROS_B.cto1);
          V[multiModeQuad_ROS_B.stemp_re_tmp_b + 3].re =
            (multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.tmp.re -
             multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.tmp.im) +
            multiModeQuad_ROS_B.cto1 * multiModeQuad_ROS_B.cfromc;
          V[multiModeQuad_ROS_B.stemp_re_tmp_b + 3].im =
            (multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.tmp.re +
             multiModeQuad_ROS_B.vtemp * multiModeQuad_ROS_B.tmp.im) +
            multiModeQuad_ROS_B.stemp_im_tmp * multiModeQuad_ROS_B.cfromc;
          multiModeQuad_ROS_B.jrow--;
        }

        multiModeQuad_ROS_B.jcol++;
      }
    }

    multiModeQuad_ROS_xzhgeqz(multiModeQuad_ROS_B.At, multiModeQuad_ROS_B.c_i_e,
      multiModeQuad_ROS_B.b_k, V, &multiModeQuad_ROS_B.jcol, D,
      multiModeQuad_ROS_B.y);
    if (multiModeQuad_ROS_B.jcol == 0) {
      multiModeQuad_ROS_xztgevc(multiModeQuad_ROS_B.At, V);
      if (multiModeQuad_ROS_B.c_i_e > 1) {
        multiModeQuad_ROS_B.c_i_e -= 2;
        while (multiModeQuad_ROS_B.c_i_e + 1 >= 1) {
          multiModeQuad_ROS_B.jcol =
            multiModeQuad_ROS_B.rscale[multiModeQuad_ROS_B.c_i_e] - 1;
          if (multiModeQuad_ROS_B.c_i_e + 1 !=
              multiModeQuad_ROS_B.rscale[multiModeQuad_ROS_B.c_i_e]) {
            multiModeQuad_ROS_B.tmp = V[multiModeQuad_ROS_B.c_i_e];
            V[multiModeQuad_ROS_B.c_i_e] = V[multiModeQuad_ROS_B.jcol];
            V[multiModeQuad_ROS_B.jcol] = multiModeQuad_ROS_B.tmp;
            multiModeQuad_ROS_B.tmp = V[multiModeQuad_ROS_B.c_i_e + 4];
            V[multiModeQuad_ROS_B.c_i_e + 4] = V[multiModeQuad_ROS_B.jcol + 4];
            V[multiModeQuad_ROS_B.jcol + 4] = multiModeQuad_ROS_B.tmp;
            multiModeQuad_ROS_B.tmp = V[multiModeQuad_ROS_B.c_i_e + 8];
            V[multiModeQuad_ROS_B.c_i_e + 8] = V[multiModeQuad_ROS_B.jcol + 8];
            V[multiModeQuad_ROS_B.jcol + 8] = multiModeQuad_ROS_B.tmp;
            multiModeQuad_ROS_B.tmp = V[multiModeQuad_ROS_B.c_i_e + 12];
            V[multiModeQuad_ROS_B.c_i_e + 12] = V[multiModeQuad_ROS_B.jcol + 12];
            V[multiModeQuad_ROS_B.jcol + 12] = multiModeQuad_ROS_B.tmp;
          }

          multiModeQuad_ROS_B.c_i_e--;
        }
      }

      if (multiModeQuad_ROS_B.b_k < 4) {
        while (multiModeQuad_ROS_B.b_k + 1 < 5) {
          multiModeQuad_ROS_B.jcol =
            multiModeQuad_ROS_B.rscale[multiModeQuad_ROS_B.b_k] - 1;
          if (multiModeQuad_ROS_B.b_k + 1 !=
              multiModeQuad_ROS_B.rscale[multiModeQuad_ROS_B.b_k]) {
            multiModeQuad_ROS_B.tmp = V[multiModeQuad_ROS_B.b_k];
            V[multiModeQuad_ROS_B.b_k] = V[multiModeQuad_ROS_B.jcol];
            V[multiModeQuad_ROS_B.jcol] = multiModeQuad_ROS_B.tmp;
            multiModeQuad_ROS_B.tmp = V[multiModeQuad_ROS_B.b_k + 4];
            V[multiModeQuad_ROS_B.b_k + 4] = V[multiModeQuad_ROS_B.jcol + 4];
            V[multiModeQuad_ROS_B.jcol + 4] = multiModeQuad_ROS_B.tmp;
            multiModeQuad_ROS_B.tmp = V[multiModeQuad_ROS_B.b_k + 8];
            V[multiModeQuad_ROS_B.b_k + 8] = V[multiModeQuad_ROS_B.jcol + 8];
            V[multiModeQuad_ROS_B.jcol + 8] = multiModeQuad_ROS_B.tmp;
            multiModeQuad_ROS_B.tmp = V[multiModeQuad_ROS_B.b_k + 12];
            V[multiModeQuad_ROS_B.b_k + 12] = V[multiModeQuad_ROS_B.jcol + 12];
            V[multiModeQuad_ROS_B.jcol + 12] = multiModeQuad_ROS_B.tmp;
          }

          multiModeQuad_ROS_B.b_k++;
        }
      }

      for (multiModeQuad_ROS_B.b_k = 0; multiModeQuad_ROS_B.b_k < 4;
           multiModeQuad_ROS_B.b_k++) {
        multiModeQuad_ROS_B.c_i_e = multiModeQuad_ROS_B.b_k << 2;
        multiModeQuad_ROS_B.cfromc = V[multiModeQuad_ROS_B.c_i_e].re;
        multiModeQuad_ROS_B.ctoc = V[multiModeQuad_ROS_B.c_i_e].im;
        multiModeQuad_ROS_B.vtemp = fabs(multiModeQuad_ROS_B.cfromc) + fabs
          (multiModeQuad_ROS_B.ctoc);
        multiModeQuad_ROS_B.cto1 = fabs(V[multiModeQuad_ROS_B.c_i_e + 1].re) +
          fabs(V[multiModeQuad_ROS_B.c_i_e + 1].im);
        if (multiModeQuad_ROS_B.cto1 > multiModeQuad_ROS_B.vtemp) {
          multiModeQuad_ROS_B.vtemp = multiModeQuad_ROS_B.cto1;
        }

        multiModeQuad_ROS_B.cto1 = fabs(V[multiModeQuad_ROS_B.c_i_e + 2].re) +
          fabs(V[multiModeQuad_ROS_B.c_i_e + 2].im);
        if (multiModeQuad_ROS_B.cto1 > multiModeQuad_ROS_B.vtemp) {
          multiModeQuad_ROS_B.vtemp = multiModeQuad_ROS_B.cto1;
        }

        multiModeQuad_ROS_B.cto1 = fabs(V[multiModeQuad_ROS_B.c_i_e + 3].re) +
          fabs(V[multiModeQuad_ROS_B.c_i_e + 3].im);
        if (multiModeQuad_ROS_B.cto1 > multiModeQuad_ROS_B.vtemp) {
          multiModeQuad_ROS_B.vtemp = multiModeQuad_ROS_B.cto1;
        }

        if (multiModeQuad_ROS_B.vtemp >= 6.7178761075670888E-139) {
          multiModeQuad_ROS_B.vtemp = 1.0 / multiModeQuad_ROS_B.vtemp;
          V[multiModeQuad_ROS_B.c_i_e].re = multiModeQuad_ROS_B.cfromc *
            multiModeQuad_ROS_B.vtemp;
          V[multiModeQuad_ROS_B.c_i_e].im = multiModeQuad_ROS_B.ctoc *
            multiModeQuad_ROS_B.vtemp;
          V[multiModeQuad_ROS_B.c_i_e + 1].re *= multiModeQuad_ROS_B.vtemp;
          V[multiModeQuad_ROS_B.c_i_e + 1].im *= multiModeQuad_ROS_B.vtemp;
          V[multiModeQuad_ROS_B.c_i_e + 2].re *= multiModeQuad_ROS_B.vtemp;
          V[multiModeQuad_ROS_B.c_i_e + 2].im *= multiModeQuad_ROS_B.vtemp;
          V[multiModeQuad_ROS_B.c_i_e + 3].re *= multiModeQuad_ROS_B.vtemp;
          V[multiModeQuad_ROS_B.c_i_e + 3].im *= multiModeQuad_ROS_B.vtemp;
        }
      }

      if (ilascl) {
        ilascl = true;
        while (ilascl) {
          multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.b_absxk *
            2.0041683600089728E-292;
          multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.anrm /
            4.9896007738368E+291;
          if ((multiModeQuad_ROS_B.cfromc > multiModeQuad_ROS_B.anrm) &&
              (multiModeQuad_ROS_B.anrm != 0.0)) {
            multiModeQuad_ROS_B.vtemp = 2.0041683600089728E-292;
            multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.cfromc;
          } else if (multiModeQuad_ROS_B.ctoc > multiModeQuad_ROS_B.b_absxk) {
            multiModeQuad_ROS_B.vtemp = 4.9896007738368E+291;
            multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.ctoc;
          } else {
            multiModeQuad_ROS_B.vtemp = multiModeQuad_ROS_B.anrm /
              multiModeQuad_ROS_B.b_absxk;
            ilascl = false;
          }

          D[0].re *= multiModeQuad_ROS_B.vtemp;
          D[0].im *= multiModeQuad_ROS_B.vtemp;
          D[1].re *= multiModeQuad_ROS_B.vtemp;
          D[1].im *= multiModeQuad_ROS_B.vtemp;
          D[2].re *= multiModeQuad_ROS_B.vtemp;
          D[2].im *= multiModeQuad_ROS_B.vtemp;
          D[3].re *= multiModeQuad_ROS_B.vtemp;
          D[3].im *= multiModeQuad_ROS_B.vtemp;
        }
      }
    }
  }

  multiModeQuad_ROS_B.anrm = 0.0;
  multiModeQuad_ROS_B.b_absxk = 3.3121686421112381E-170;
  multiModeQuad_ROS_B.b_k = 0;
  while (multiModeQuad_ROS_B.b_k + 1 <= 4) {
    multiModeQuad_ROS_B.cfromc = fabs(V[multiModeQuad_ROS_B.b_k].re);
    if (multiModeQuad_ROS_B.cfromc > multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.anrm *
        multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.ctoc + 1.0;
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.cfromc;
    } else {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.anrm += multiModeQuad_ROS_B.ctoc *
        multiModeQuad_ROS_B.ctoc;
    }

    multiModeQuad_ROS_B.cfromc = fabs(V[multiModeQuad_ROS_B.b_k].im);
    if (multiModeQuad_ROS_B.cfromc > multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.anrm *
        multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.ctoc + 1.0;
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.cfromc;
    } else {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.anrm += multiModeQuad_ROS_B.ctoc *
        multiModeQuad_ROS_B.ctoc;
    }

    multiModeQuad_ROS_B.b_k++;
  }

  multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.b_absxk * sqrt
    (multiModeQuad_ROS_B.anrm);
  multiModeQuad_ROS_B.b_k = 0;
  while (multiModeQuad_ROS_B.b_k + 1 <= 4) {
    multiModeQuad_ROS_B.b_absxk = V[multiModeQuad_ROS_B.b_k].re;
    multiModeQuad_ROS_B.cfromc = V[multiModeQuad_ROS_B.b_k].im;
    if (multiModeQuad_ROS_B.cfromc == 0.0) {
      V[multiModeQuad_ROS_B.b_k].re = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.anrm;
      V[multiModeQuad_ROS_B.b_k].im = 0.0;
    } else if (multiModeQuad_ROS_B.b_absxk == 0.0) {
      V[multiModeQuad_ROS_B.b_k].re = 0.0;
      V[multiModeQuad_ROS_B.b_k].im = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.anrm;
    } else {
      V[multiModeQuad_ROS_B.b_k].re = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.anrm;
      V[multiModeQuad_ROS_B.b_k].im = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.anrm;
    }

    multiModeQuad_ROS_B.b_k++;
  }

  multiModeQuad_ROS_B.anrm = 0.0;
  multiModeQuad_ROS_B.b_absxk = 3.3121686421112381E-170;
  multiModeQuad_ROS_B.b_k = 4;
  while (multiModeQuad_ROS_B.b_k + 1 <= 8) {
    multiModeQuad_ROS_B.cfromc = fabs(V[multiModeQuad_ROS_B.b_k].re);
    if (multiModeQuad_ROS_B.cfromc > multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.anrm *
        multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.ctoc + 1.0;
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.cfromc;
    } else {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.anrm += multiModeQuad_ROS_B.ctoc *
        multiModeQuad_ROS_B.ctoc;
    }

    multiModeQuad_ROS_B.cfromc = fabs(V[multiModeQuad_ROS_B.b_k].im);
    if (multiModeQuad_ROS_B.cfromc > multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.anrm *
        multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.ctoc + 1.0;
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.cfromc;
    } else {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.anrm += multiModeQuad_ROS_B.ctoc *
        multiModeQuad_ROS_B.ctoc;
    }

    multiModeQuad_ROS_B.b_k++;
  }

  multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.b_absxk * sqrt
    (multiModeQuad_ROS_B.anrm);
  multiModeQuad_ROS_B.b_k = 4;
  while (multiModeQuad_ROS_B.b_k + 1 <= 8) {
    multiModeQuad_ROS_B.b_absxk = V[multiModeQuad_ROS_B.b_k].re;
    multiModeQuad_ROS_B.cfromc = V[multiModeQuad_ROS_B.b_k].im;
    if (multiModeQuad_ROS_B.cfromc == 0.0) {
      V[multiModeQuad_ROS_B.b_k].re = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.anrm;
      V[multiModeQuad_ROS_B.b_k].im = 0.0;
    } else if (multiModeQuad_ROS_B.b_absxk == 0.0) {
      V[multiModeQuad_ROS_B.b_k].re = 0.0;
      V[multiModeQuad_ROS_B.b_k].im = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.anrm;
    } else {
      V[multiModeQuad_ROS_B.b_k].re = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.anrm;
      V[multiModeQuad_ROS_B.b_k].im = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.anrm;
    }

    multiModeQuad_ROS_B.b_k++;
  }

  multiModeQuad_ROS_B.anrm = 0.0;
  multiModeQuad_ROS_B.b_absxk = 3.3121686421112381E-170;
  multiModeQuad_ROS_B.b_k = 8;
  while (multiModeQuad_ROS_B.b_k + 1 <= 12) {
    multiModeQuad_ROS_B.cfromc = fabs(V[multiModeQuad_ROS_B.b_k].re);
    if (multiModeQuad_ROS_B.cfromc > multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.anrm *
        multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.ctoc + 1.0;
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.cfromc;
    } else {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.anrm += multiModeQuad_ROS_B.ctoc *
        multiModeQuad_ROS_B.ctoc;
    }

    multiModeQuad_ROS_B.cfromc = fabs(V[multiModeQuad_ROS_B.b_k].im);
    if (multiModeQuad_ROS_B.cfromc > multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.anrm *
        multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.ctoc + 1.0;
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.cfromc;
    } else {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.anrm += multiModeQuad_ROS_B.ctoc *
        multiModeQuad_ROS_B.ctoc;
    }

    multiModeQuad_ROS_B.b_k++;
  }

  multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.b_absxk * sqrt
    (multiModeQuad_ROS_B.anrm);
  multiModeQuad_ROS_B.b_k = 8;
  while (multiModeQuad_ROS_B.b_k + 1 <= 12) {
    multiModeQuad_ROS_B.b_absxk = V[multiModeQuad_ROS_B.b_k].re;
    multiModeQuad_ROS_B.cfromc = V[multiModeQuad_ROS_B.b_k].im;
    if (multiModeQuad_ROS_B.cfromc == 0.0) {
      V[multiModeQuad_ROS_B.b_k].re = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.anrm;
      V[multiModeQuad_ROS_B.b_k].im = 0.0;
    } else if (multiModeQuad_ROS_B.b_absxk == 0.0) {
      V[multiModeQuad_ROS_B.b_k].re = 0.0;
      V[multiModeQuad_ROS_B.b_k].im = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.anrm;
    } else {
      V[multiModeQuad_ROS_B.b_k].re = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.anrm;
      V[multiModeQuad_ROS_B.b_k].im = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.anrm;
    }

    multiModeQuad_ROS_B.b_k++;
  }

  multiModeQuad_ROS_B.anrm = 0.0;
  multiModeQuad_ROS_B.b_absxk = 3.3121686421112381E-170;
  multiModeQuad_ROS_B.b_k = 12;
  while (multiModeQuad_ROS_B.b_k + 1 <= 16) {
    multiModeQuad_ROS_B.cfromc = fabs(V[multiModeQuad_ROS_B.b_k].re);
    if (multiModeQuad_ROS_B.cfromc > multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.anrm *
        multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.ctoc + 1.0;
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.cfromc;
    } else {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.anrm += multiModeQuad_ROS_B.ctoc *
        multiModeQuad_ROS_B.ctoc;
    }

    multiModeQuad_ROS_B.cfromc = fabs(V[multiModeQuad_ROS_B.b_k].im);
    if (multiModeQuad_ROS_B.cfromc > multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.anrm *
        multiModeQuad_ROS_B.ctoc * multiModeQuad_ROS_B.ctoc + 1.0;
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.cfromc;
    } else {
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.anrm += multiModeQuad_ROS_B.ctoc *
        multiModeQuad_ROS_B.ctoc;
    }

    multiModeQuad_ROS_B.b_k++;
  }

  multiModeQuad_ROS_B.anrm = multiModeQuad_ROS_B.b_absxk * sqrt
    (multiModeQuad_ROS_B.anrm);
  multiModeQuad_ROS_B.b_k = 12;
  while (multiModeQuad_ROS_B.b_k + 1 <= 16) {
    multiModeQuad_ROS_B.b_absxk = V[multiModeQuad_ROS_B.b_k].re;
    multiModeQuad_ROS_B.cfromc = V[multiModeQuad_ROS_B.b_k].im;
    if (multiModeQuad_ROS_B.cfromc == 0.0) {
      V[multiModeQuad_ROS_B.b_k].re = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.anrm;
      V[multiModeQuad_ROS_B.b_k].im = 0.0;
    } else if (multiModeQuad_ROS_B.b_absxk == 0.0) {
      V[multiModeQuad_ROS_B.b_k].re = 0.0;
      V[multiModeQuad_ROS_B.b_k].im = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.anrm;
    } else {
      V[multiModeQuad_ROS_B.b_k].re = multiModeQuad_ROS_B.b_absxk /
        multiModeQuad_ROS_B.anrm;
      V[multiModeQuad_ROS_B.b_k].im = multiModeQuad_ROS_B.cfromc /
        multiModeQuad_ROS_B.anrm;
    }

    multiModeQuad_ROS_B.b_k++;
  }

  if (multiModeQuad_ROS_B.y[0].im == 0.0) {
    if (D[0].im == 0.0) {
      multiModeQuad_ROS_B.anrm = D[0].re / multiModeQuad_ROS_B.y[0].re;
      multiModeQuad_ROS_B.b_absxk = 0.0;
    } else if (D[0].re == 0.0) {
      multiModeQuad_ROS_B.anrm = 0.0;
      multiModeQuad_ROS_B.b_absxk = D[0].im / multiModeQuad_ROS_B.y[0].re;
    } else {
      multiModeQuad_ROS_B.anrm = D[0].re / multiModeQuad_ROS_B.y[0].re;
      multiModeQuad_ROS_B.b_absxk = D[0].im / multiModeQuad_ROS_B.y[0].re;
    }
  } else if (multiModeQuad_ROS_B.y[0].re == 0.0) {
    if (D[0].re == 0.0) {
      multiModeQuad_ROS_B.anrm = D[0].im / multiModeQuad_ROS_B.y[0].im;
      multiModeQuad_ROS_B.b_absxk = 0.0;
    } else if (D[0].im == 0.0) {
      multiModeQuad_ROS_B.anrm = 0.0;
      multiModeQuad_ROS_B.b_absxk = -(D[0].re / multiModeQuad_ROS_B.y[0].im);
    } else {
      multiModeQuad_ROS_B.anrm = D[0].im / multiModeQuad_ROS_B.y[0].im;
      multiModeQuad_ROS_B.b_absxk = -(D[0].re / multiModeQuad_ROS_B.y[0].im);
    }
  } else {
    multiModeQuad_ROS_B.b_absxk = fabs(multiModeQuad_ROS_B.y[0].re);
    multiModeQuad_ROS_B.anrm = fabs(multiModeQuad_ROS_B.y[0].im);
    if (multiModeQuad_ROS_B.b_absxk > multiModeQuad_ROS_B.anrm) {
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.y[0].im /
        multiModeQuad_ROS_B.y[0].re;
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.b_absxk *
        multiModeQuad_ROS_B.y[0].im + multiModeQuad_ROS_B.y[0].re;
      multiModeQuad_ROS_B.anrm = (multiModeQuad_ROS_B.b_absxk * D[0].im + D[0].
        re) / multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.b_absxk = (D[0].im - multiModeQuad_ROS_B.b_absxk * D[0]
        .re) / multiModeQuad_ROS_B.cfromc;
    } else if (multiModeQuad_ROS_B.anrm == multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.y[0].re > 0.0 ? 0.5 :
        -0.5;
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.y[0].im > 0.0 ? 0.5 : -0.5;
      multiModeQuad_ROS_B.anrm = (D[0].re * multiModeQuad_ROS_B.cfromc + D[0].im
        * multiModeQuad_ROS_B.ctoc) / multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.b_absxk = (D[0].im * multiModeQuad_ROS_B.cfromc - D[0]
        .re * multiModeQuad_ROS_B.ctoc) / multiModeQuad_ROS_B.b_absxk;
    } else {
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.y[0].re /
        multiModeQuad_ROS_B.y[0].im;
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.b_absxk *
        multiModeQuad_ROS_B.y[0].re + multiModeQuad_ROS_B.y[0].im;
      multiModeQuad_ROS_B.anrm = (multiModeQuad_ROS_B.b_absxk * D[0].re + D[0].
        im) / multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.b_absxk = (multiModeQuad_ROS_B.b_absxk * D[0].im - D[0]
        .re) / multiModeQuad_ROS_B.cfromc;
    }
  }

  D[0].re = multiModeQuad_ROS_B.anrm;
  D[0].im = multiModeQuad_ROS_B.b_absxk;
  if (multiModeQuad_ROS_B.y[1].im == 0.0) {
    if (D[1].im == 0.0) {
      multiModeQuad_ROS_B.anrm = D[1].re / multiModeQuad_ROS_B.y[1].re;
      multiModeQuad_ROS_B.b_absxk = 0.0;
    } else if (D[1].re == 0.0) {
      multiModeQuad_ROS_B.anrm = 0.0;
      multiModeQuad_ROS_B.b_absxk = D[1].im / multiModeQuad_ROS_B.y[1].re;
    } else {
      multiModeQuad_ROS_B.anrm = D[1].re / multiModeQuad_ROS_B.y[1].re;
      multiModeQuad_ROS_B.b_absxk = D[1].im / multiModeQuad_ROS_B.y[1].re;
    }
  } else if (multiModeQuad_ROS_B.y[1].re == 0.0) {
    if (D[1].re == 0.0) {
      multiModeQuad_ROS_B.anrm = D[1].im / multiModeQuad_ROS_B.y[1].im;
      multiModeQuad_ROS_B.b_absxk = 0.0;
    } else if (D[1].im == 0.0) {
      multiModeQuad_ROS_B.anrm = 0.0;
      multiModeQuad_ROS_B.b_absxk = -(D[1].re / multiModeQuad_ROS_B.y[1].im);
    } else {
      multiModeQuad_ROS_B.anrm = D[1].im / multiModeQuad_ROS_B.y[1].im;
      multiModeQuad_ROS_B.b_absxk = -(D[1].re / multiModeQuad_ROS_B.y[1].im);
    }
  } else {
    multiModeQuad_ROS_B.b_absxk = fabs(multiModeQuad_ROS_B.y[1].re);
    multiModeQuad_ROS_B.anrm = fabs(multiModeQuad_ROS_B.y[1].im);
    if (multiModeQuad_ROS_B.b_absxk > multiModeQuad_ROS_B.anrm) {
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.y[1].im /
        multiModeQuad_ROS_B.y[1].re;
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.b_absxk *
        multiModeQuad_ROS_B.y[1].im + multiModeQuad_ROS_B.y[1].re;
      multiModeQuad_ROS_B.anrm = (multiModeQuad_ROS_B.b_absxk * D[1].im + D[1].
        re) / multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.b_absxk = (D[1].im - multiModeQuad_ROS_B.b_absxk * D[1]
        .re) / multiModeQuad_ROS_B.cfromc;
    } else if (multiModeQuad_ROS_B.anrm == multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.y[1].re > 0.0 ? 0.5 :
        -0.5;
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.y[1].im > 0.0 ? 0.5 : -0.5;
      multiModeQuad_ROS_B.anrm = (D[1].re * multiModeQuad_ROS_B.cfromc + D[1].im
        * multiModeQuad_ROS_B.ctoc) / multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.b_absxk = (D[1].im * multiModeQuad_ROS_B.cfromc - D[1]
        .re * multiModeQuad_ROS_B.ctoc) / multiModeQuad_ROS_B.b_absxk;
    } else {
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.y[1].re /
        multiModeQuad_ROS_B.y[1].im;
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.b_absxk *
        multiModeQuad_ROS_B.y[1].re + multiModeQuad_ROS_B.y[1].im;
      multiModeQuad_ROS_B.anrm = (multiModeQuad_ROS_B.b_absxk * D[1].re + D[1].
        im) / multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.b_absxk = (multiModeQuad_ROS_B.b_absxk * D[1].im - D[1]
        .re) / multiModeQuad_ROS_B.cfromc;
    }
  }

  D[1].re = multiModeQuad_ROS_B.anrm;
  D[1].im = multiModeQuad_ROS_B.b_absxk;
  if (multiModeQuad_ROS_B.y[2].im == 0.0) {
    if (D[2].im == 0.0) {
      multiModeQuad_ROS_B.anrm = D[2].re / multiModeQuad_ROS_B.y[2].re;
      multiModeQuad_ROS_B.b_absxk = 0.0;
    } else if (D[2].re == 0.0) {
      multiModeQuad_ROS_B.anrm = 0.0;
      multiModeQuad_ROS_B.b_absxk = D[2].im / multiModeQuad_ROS_B.y[2].re;
    } else {
      multiModeQuad_ROS_B.anrm = D[2].re / multiModeQuad_ROS_B.y[2].re;
      multiModeQuad_ROS_B.b_absxk = D[2].im / multiModeQuad_ROS_B.y[2].re;
    }
  } else if (multiModeQuad_ROS_B.y[2].re == 0.0) {
    if (D[2].re == 0.0) {
      multiModeQuad_ROS_B.anrm = D[2].im / multiModeQuad_ROS_B.y[2].im;
      multiModeQuad_ROS_B.b_absxk = 0.0;
    } else if (D[2].im == 0.0) {
      multiModeQuad_ROS_B.anrm = 0.0;
      multiModeQuad_ROS_B.b_absxk = -(D[2].re / multiModeQuad_ROS_B.y[2].im);
    } else {
      multiModeQuad_ROS_B.anrm = D[2].im / multiModeQuad_ROS_B.y[2].im;
      multiModeQuad_ROS_B.b_absxk = -(D[2].re / multiModeQuad_ROS_B.y[2].im);
    }
  } else {
    multiModeQuad_ROS_B.b_absxk = fabs(multiModeQuad_ROS_B.y[2].re);
    multiModeQuad_ROS_B.anrm = fabs(multiModeQuad_ROS_B.y[2].im);
    if (multiModeQuad_ROS_B.b_absxk > multiModeQuad_ROS_B.anrm) {
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.y[2].im /
        multiModeQuad_ROS_B.y[2].re;
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.b_absxk *
        multiModeQuad_ROS_B.y[2].im + multiModeQuad_ROS_B.y[2].re;
      multiModeQuad_ROS_B.anrm = (multiModeQuad_ROS_B.b_absxk * D[2].im + D[2].
        re) / multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.b_absxk = (D[2].im - multiModeQuad_ROS_B.b_absxk * D[2]
        .re) / multiModeQuad_ROS_B.cfromc;
    } else if (multiModeQuad_ROS_B.anrm == multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.y[2].re > 0.0 ? 0.5 :
        -0.5;
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.y[2].im > 0.0 ? 0.5 : -0.5;
      multiModeQuad_ROS_B.anrm = (D[2].re * multiModeQuad_ROS_B.cfromc + D[2].im
        * multiModeQuad_ROS_B.ctoc) / multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.b_absxk = (D[2].im * multiModeQuad_ROS_B.cfromc - D[2]
        .re * multiModeQuad_ROS_B.ctoc) / multiModeQuad_ROS_B.b_absxk;
    } else {
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.y[2].re /
        multiModeQuad_ROS_B.y[2].im;
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.b_absxk *
        multiModeQuad_ROS_B.y[2].re + multiModeQuad_ROS_B.y[2].im;
      multiModeQuad_ROS_B.anrm = (multiModeQuad_ROS_B.b_absxk * D[2].re + D[2].
        im) / multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.b_absxk = (multiModeQuad_ROS_B.b_absxk * D[2].im - D[2]
        .re) / multiModeQuad_ROS_B.cfromc;
    }
  }

  D[2].re = multiModeQuad_ROS_B.anrm;
  D[2].im = multiModeQuad_ROS_B.b_absxk;
  if (multiModeQuad_ROS_B.y[3].im == 0.0) {
    if (D[3].im == 0.0) {
      multiModeQuad_ROS_B.anrm = D[3].re / multiModeQuad_ROS_B.y[3].re;
      multiModeQuad_ROS_B.b_absxk = 0.0;
    } else if (D[3].re == 0.0) {
      multiModeQuad_ROS_B.anrm = 0.0;
      multiModeQuad_ROS_B.b_absxk = D[3].im / multiModeQuad_ROS_B.y[3].re;
    } else {
      multiModeQuad_ROS_B.anrm = D[3].re / multiModeQuad_ROS_B.y[3].re;
      multiModeQuad_ROS_B.b_absxk = D[3].im / multiModeQuad_ROS_B.y[3].re;
    }
  } else if (multiModeQuad_ROS_B.y[3].re == 0.0) {
    if (D[3].re == 0.0) {
      multiModeQuad_ROS_B.anrm = D[3].im / multiModeQuad_ROS_B.y[3].im;
      multiModeQuad_ROS_B.b_absxk = 0.0;
    } else if (D[3].im == 0.0) {
      multiModeQuad_ROS_B.anrm = 0.0;
      multiModeQuad_ROS_B.b_absxk = -(D[3].re / multiModeQuad_ROS_B.y[3].im);
    } else {
      multiModeQuad_ROS_B.anrm = D[3].im / multiModeQuad_ROS_B.y[3].im;
      multiModeQuad_ROS_B.b_absxk = -(D[3].re / multiModeQuad_ROS_B.y[3].im);
    }
  } else {
    multiModeQuad_ROS_B.b_absxk = fabs(multiModeQuad_ROS_B.y[3].re);
    multiModeQuad_ROS_B.anrm = fabs(multiModeQuad_ROS_B.y[3].im);
    if (multiModeQuad_ROS_B.b_absxk > multiModeQuad_ROS_B.anrm) {
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.y[3].im /
        multiModeQuad_ROS_B.y[3].re;
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.b_absxk *
        multiModeQuad_ROS_B.y[3].im + multiModeQuad_ROS_B.y[3].re;
      multiModeQuad_ROS_B.anrm = (multiModeQuad_ROS_B.b_absxk * D[3].im + D[3].
        re) / multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.b_absxk = (D[3].im - multiModeQuad_ROS_B.b_absxk * D[3]
        .re) / multiModeQuad_ROS_B.cfromc;
    } else if (multiModeQuad_ROS_B.anrm == multiModeQuad_ROS_B.b_absxk) {
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.y[3].re > 0.0 ? 0.5 :
        -0.5;
      multiModeQuad_ROS_B.ctoc = multiModeQuad_ROS_B.y[3].im > 0.0 ? 0.5 : -0.5;
      multiModeQuad_ROS_B.anrm = (D[3].re * multiModeQuad_ROS_B.cfromc + D[3].im
        * multiModeQuad_ROS_B.ctoc) / multiModeQuad_ROS_B.b_absxk;
      multiModeQuad_ROS_B.b_absxk = (D[3].im * multiModeQuad_ROS_B.cfromc - D[3]
        .re * multiModeQuad_ROS_B.ctoc) / multiModeQuad_ROS_B.b_absxk;
    } else {
      multiModeQuad_ROS_B.b_absxk = multiModeQuad_ROS_B.y[3].re /
        multiModeQuad_ROS_B.y[3].im;
      multiModeQuad_ROS_B.cfromc = multiModeQuad_ROS_B.b_absxk *
        multiModeQuad_ROS_B.y[3].re + multiModeQuad_ROS_B.y[3].im;
      multiModeQuad_ROS_B.anrm = (multiModeQuad_ROS_B.b_absxk * D[3].re + D[3].
        im) / multiModeQuad_ROS_B.cfromc;
      multiModeQuad_ROS_B.b_absxk = (multiModeQuad_ROS_B.b_absxk * D[3].im - D[3]
        .re) / multiModeQuad_ROS_B.cfromc;
    }
  }

  D[3].re = multiModeQuad_ROS_B.anrm;
  D[3].im = multiModeQuad_ROS_B.b_absxk;
}

// Function for MATLAB Function: '<S10>/Attitude control'
static real_T multiModeQuad_ROS_xnrm2(int32_T n, const real_T x[16], int32_T ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      real_T scale;
      scale = 3.3121686421112381E-170;
      for (int32_T k = ix0; k <= ix0 + 1; k++) {
        real_T absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_xzlarf(int32_T m, int32_T n, int32_T iv0, real_T
  tau, real_T C[16], int32_T ic0, real_T work[4])
{
  int32_T coltop;
  int32_T jy;
  int32_T lastc;
  int32_T lastv;
  if (tau != 0.0) {
    boolean_T exitg2;
    lastv = m;
    lastc = iv0 + m;
    while ((lastv > 0) && (C[lastc - 2] == 0.0)) {
      lastv--;
      lastc--;
    }

    lastc = n;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      int32_T exitg1;
      coltop = ((lastc - 1) << 2) + ic0;
      jy = coltop;
      do {
        exitg1 = 0;
        if (jy <= (coltop + lastv) - 1) {
          if (C[jy - 1] != 0.0) {
            exitg1 = 1;
          } else {
            jy++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }

    lastc--;
  } else {
    lastv = 0;
    lastc = -1;
  }

  if (lastv > 0) {
    real_T c;
    int32_T d;
    int32_T ia;
    int32_T iac;
    int32_T ix;
    if (lastc + 1 != 0) {
      for (coltop = 0; coltop <= lastc; coltop++) {
        work[coltop] = 0.0;
      }

      coltop = 0;
      jy = (lastc << 2) + ic0;
      for (iac = ic0; iac <= jy; iac += 4) {
        ix = iv0;
        c = 0.0;
        d = (iac + lastv) - 1;
        for (ia = iac; ia <= d; ia++) {
          c += C[ia - 1] * C[ix - 1];
          ix++;
        }

        work[coltop] += c;
        coltop++;
      }
    }

    if (!(-tau == 0.0)) {
      coltop = ic0 - 1;
      jy = 0;
      for (iac = 0; iac <= lastc; iac++) {
        if (work[jy] != 0.0) {
          c = work[jy] * -tau;
          ix = iv0;
          d = coltop;
          ia = lastv + coltop;
          while (d + 1 <= ia) {
            C[d] += C[ix - 1] * c;
            ix++;
            d++;
          }
        }

        jy++;
        coltop += 4;
      }
    }
  }
}

// Function for MATLAB Function: '<S10>/Attitude control'
static real_T multiModeQuad_ROS_xnrm2_l(int32_T n, const real_T x[3])
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[1]);
    } else {
      multiModeQuad_ROS_B.scale_b = 3.3121686421112381E-170;
      multiModeQuad_ROS_B.absxk_l = fabs(x[1]);
      if (multiModeQuad_ROS_B.absxk_l > 3.3121686421112381E-170) {
        y = 1.0;
        multiModeQuad_ROS_B.scale_b = multiModeQuad_ROS_B.absxk_l;
      } else {
        multiModeQuad_ROS_B.t_h = multiModeQuad_ROS_B.absxk_l /
          3.3121686421112381E-170;
        y = multiModeQuad_ROS_B.t_h * multiModeQuad_ROS_B.t_h;
      }

      multiModeQuad_ROS_B.absxk_l = fabs(x[2]);
      if (multiModeQuad_ROS_B.absxk_l > multiModeQuad_ROS_B.scale_b) {
        multiModeQuad_ROS_B.t_h = multiModeQuad_ROS_B.scale_b /
          multiModeQuad_ROS_B.absxk_l;
        y = y * multiModeQuad_ROS_B.t_h * multiModeQuad_ROS_B.t_h + 1.0;
        multiModeQuad_ROS_B.scale_b = multiModeQuad_ROS_B.absxk_l;
      } else {
        multiModeQuad_ROS_B.t_h = multiModeQuad_ROS_B.absxk_l /
          multiModeQuad_ROS_B.scale_b;
        y += multiModeQuad_ROS_B.t_h * multiModeQuad_ROS_B.t_h;
      }

      y = multiModeQuad_ROS_B.scale_b * sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S10>/Attitude control'
static real_T multiModeQuad_ROS_xzlarfg(int32_T n, real_T *alpha1, real_T x[3])
{
  real_T tau;
  tau = 0.0;
  if (n > 0) {
    multiModeQuad_ROS_B.xnorm = multiModeQuad_ROS_xnrm2_l(n - 1, x);
    if (multiModeQuad_ROS_B.xnorm != 0.0) {
      multiModeQuad_ROS_B.xnorm = multiModeQuad_ROS_rt_hypotd_snf(*alpha1,
        multiModeQuad_ROS_B.xnorm);
      if (*alpha1 >= 0.0) {
        multiModeQuad_ROS_B.xnorm = -multiModeQuad_ROS_B.xnorm;
      }

      if (fabs(multiModeQuad_ROS_B.xnorm) < 1.0020841800044864E-292) {
        int32_T c_k;
        int32_T knt;
        knt = -1;
        do {
          knt++;
          for (c_k = 1; c_k < n; c_k++) {
            x[c_k] *= 9.9792015476736E+291;
          }

          multiModeQuad_ROS_B.xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(multiModeQuad_ROS_B.xnorm) >= 1.0020841800044864E-292));

        multiModeQuad_ROS_B.xnorm = multiModeQuad_ROS_rt_hypotd_snf(*alpha1,
          multiModeQuad_ROS_xnrm2_l(n - 1, x));
        if (*alpha1 >= 0.0) {
          multiModeQuad_ROS_B.xnorm = -multiModeQuad_ROS_B.xnorm;
        }

        tau = (multiModeQuad_ROS_B.xnorm - *alpha1) / multiModeQuad_ROS_B.xnorm;
        multiModeQuad_ROS_B.a_n = 1.0 / (*alpha1 - multiModeQuad_ROS_B.xnorm);
        for (c_k = 1; c_k < n; c_k++) {
          x[c_k] *= multiModeQuad_ROS_B.a_n;
        }

        for (c_k = 0; c_k <= knt; c_k++) {
          multiModeQuad_ROS_B.xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = multiModeQuad_ROS_B.xnorm;
      } else {
        tau = (multiModeQuad_ROS_B.xnorm - *alpha1) / multiModeQuad_ROS_B.xnorm;
        multiModeQuad_ROS_B.a_n = 1.0 / (*alpha1 - multiModeQuad_ROS_B.xnorm);
        for (int32_T knt = 1; knt < n; knt++) {
          x[knt] *= multiModeQuad_ROS_B.a_n;
        }

        *alpha1 = multiModeQuad_ROS_B.xnorm;
      }
    }
  }

  return tau;
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_xdlanv2(real_T *a, real_T *b, real_T *c, real_T *d,
  real_T *rt1r, real_T *rt1i, real_T *rt2r, real_T *rt2i, real_T *cs, real_T *sn)
{
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    multiModeQuad_ROS_B.bcmax = *d;
    *d = *a;
    *a = multiModeQuad_ROS_B.bcmax;
    *b = -*c;
    *c = 0.0;
  } else {
    multiModeQuad_ROS_B.tau = *a - *d;
    if ((multiModeQuad_ROS_B.tau == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
      *cs = 1.0;
      *sn = 0.0;
    } else {
      int32_T b_0;
      int32_T c_0;
      boolean_T tmp;
      multiModeQuad_ROS_B.p = 0.5 * multiModeQuad_ROS_B.tau;
      multiModeQuad_ROS_B.bcmis = fabs(*b);
      multiModeQuad_ROS_B.z = fabs(*c);
      tmp = rtIsNaN(multiModeQuad_ROS_B.z);
      if ((multiModeQuad_ROS_B.bcmis >= multiModeQuad_ROS_B.z) || tmp) {
        multiModeQuad_ROS_B.bcmax = multiModeQuad_ROS_B.bcmis;
      } else {
        multiModeQuad_ROS_B.bcmax = multiModeQuad_ROS_B.z;
      }

      if ((multiModeQuad_ROS_B.bcmis <= multiModeQuad_ROS_B.z) || tmp) {
        multiModeQuad_ROS_B.z = multiModeQuad_ROS_B.bcmis;
      }

      if (!(*b < 0.0)) {
        b_0 = 1;
      } else {
        b_0 = -1;
      }

      if (!(*c < 0.0)) {
        c_0 = 1;
      } else {
        c_0 = -1;
      }

      multiModeQuad_ROS_B.bcmis = multiModeQuad_ROS_B.z * static_cast<real_T>
        (b_0) * static_cast<real_T>(c_0);
      multiModeQuad_ROS_B.scale_p = fabs(multiModeQuad_ROS_B.p);
      if ((!(multiModeQuad_ROS_B.scale_p >= multiModeQuad_ROS_B.bcmax)) &&
          (!rtIsNaN(multiModeQuad_ROS_B.bcmax))) {
        multiModeQuad_ROS_B.scale_p = multiModeQuad_ROS_B.bcmax;
      }

      multiModeQuad_ROS_B.z = multiModeQuad_ROS_B.p /
        multiModeQuad_ROS_B.scale_p * multiModeQuad_ROS_B.p +
        multiModeQuad_ROS_B.bcmax / multiModeQuad_ROS_B.scale_p *
        multiModeQuad_ROS_B.bcmis;
      if (multiModeQuad_ROS_B.z >= 8.8817841970012523E-16) {
        if (!(multiModeQuad_ROS_B.p < 0.0)) {
          multiModeQuad_ROS_B.tau = sqrt(multiModeQuad_ROS_B.scale_p) * sqrt
            (multiModeQuad_ROS_B.z);
        } else {
          multiModeQuad_ROS_B.tau = -(sqrt(multiModeQuad_ROS_B.scale_p) * sqrt
            (multiModeQuad_ROS_B.z));
        }

        multiModeQuad_ROS_B.z = multiModeQuad_ROS_B.p + multiModeQuad_ROS_B.tau;
        *a = *d + multiModeQuad_ROS_B.z;
        *d -= multiModeQuad_ROS_B.bcmax / multiModeQuad_ROS_B.z *
          multiModeQuad_ROS_B.bcmis;
        multiModeQuad_ROS_B.tau = multiModeQuad_ROS_rt_hypotd_snf(*c,
          multiModeQuad_ROS_B.z);
        *cs = multiModeQuad_ROS_B.z / multiModeQuad_ROS_B.tau;
        *sn = *c / multiModeQuad_ROS_B.tau;
        *b -= *c;
        *c = 0.0;
      } else {
        multiModeQuad_ROS_B.bcmax = *b + *c;
        multiModeQuad_ROS_B.tau = multiModeQuad_ROS_rt_hypotd_snf
          (multiModeQuad_ROS_B.bcmax, multiModeQuad_ROS_B.tau);
        *cs = sqrt((fabs(multiModeQuad_ROS_B.bcmax) / multiModeQuad_ROS_B.tau +
                    1.0) * 0.5);
        if (!(multiModeQuad_ROS_B.bcmax < 0.0)) {
          b_0 = 1;
        } else {
          b_0 = -1;
        }

        *sn = -(multiModeQuad_ROS_B.p / (multiModeQuad_ROS_B.tau * *cs)) *
          static_cast<real_T>(b_0);
        multiModeQuad_ROS_B.p = *a * *cs + *b * *sn;
        multiModeQuad_ROS_B.tau = -*a * *sn + *b * *cs;
        multiModeQuad_ROS_B.bcmax = *c * *cs + *d * *sn;
        multiModeQuad_ROS_B.bcmis = -*c * *sn + *d * *cs;
        *b = multiModeQuad_ROS_B.tau * *cs + multiModeQuad_ROS_B.bcmis * *sn;
        *c = -multiModeQuad_ROS_B.p * *sn + multiModeQuad_ROS_B.bcmax * *cs;
        multiModeQuad_ROS_B.bcmax = ((multiModeQuad_ROS_B.p * *cs +
          multiModeQuad_ROS_B.bcmax * *sn) + (-multiModeQuad_ROS_B.tau * *sn +
          multiModeQuad_ROS_B.bcmis * *cs)) * 0.5;
        *a = multiModeQuad_ROS_B.bcmax;
        *d = multiModeQuad_ROS_B.bcmax;
        if (*c != 0.0) {
          if (*b != 0.0) {
            if ((*b < 0.0) == (*c < 0.0)) {
              multiModeQuad_ROS_B.z = sqrt(fabs(*b));
              multiModeQuad_ROS_B.bcmis = sqrt(fabs(*c));
              if (!(*c < 0.0)) {
                multiModeQuad_ROS_B.p = multiModeQuad_ROS_B.z *
                  multiModeQuad_ROS_B.bcmis;
              } else {
                multiModeQuad_ROS_B.p = -(multiModeQuad_ROS_B.z *
                  multiModeQuad_ROS_B.bcmis);
              }

              multiModeQuad_ROS_B.tau = 1.0 / sqrt(fabs(*b + *c));
              *a = multiModeQuad_ROS_B.bcmax + multiModeQuad_ROS_B.p;
              *d = multiModeQuad_ROS_B.bcmax - multiModeQuad_ROS_B.p;
              *b -= *c;
              *c = 0.0;
              multiModeQuad_ROS_B.p = multiModeQuad_ROS_B.z *
                multiModeQuad_ROS_B.tau;
              multiModeQuad_ROS_B.tau *= multiModeQuad_ROS_B.bcmis;
              multiModeQuad_ROS_B.bcmax = *cs * multiModeQuad_ROS_B.p - *sn *
                multiModeQuad_ROS_B.tau;
              *sn = *cs * multiModeQuad_ROS_B.tau + *sn * multiModeQuad_ROS_B.p;
              *cs = multiModeQuad_ROS_B.bcmax;
            }
          } else {
            *b = -*c;
            *c = 0.0;
            multiModeQuad_ROS_B.bcmax = *cs;
            *cs = -*sn;
            *sn = multiModeQuad_ROS_B.bcmax;
          }
        }
      }
    }
  }

  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = sqrt(fabs(*b)) * sqrt(fabs(*c));
    *rt2i = -*rt1i;
  }
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_xrot(int32_T n, real_T x[16], int32_T ix0, int32_T
  iy0, real_T c, real_T s)
{
  if (n >= 1) {
    int32_T ix;
    int32_T iy;
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (int32_T k = 0; k < n; k++) {
      multiModeQuad_ROS_B.temp_b = c * x[ix] + s * x[iy];
      x[iy] = c * x[iy] - s * x[ix];
      x[ix] = multiModeQuad_ROS_B.temp_b;
      iy++;
      ix++;
    }
  }
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_xrot_i(real_T x[16], int32_T ix0, int32_T iy0,
  real_T c, real_T s)
{
  real_T temp_tmp;
  multiModeQuad_ROS_B.temp_d = x[iy0 - 1];
  temp_tmp = x[ix0 - 1];
  x[iy0 - 1] = multiModeQuad_ROS_B.temp_d * c - temp_tmp * s;
  x[ix0 - 1] = temp_tmp * c + multiModeQuad_ROS_B.temp_d * s;
  multiModeQuad_ROS_B.temp_d = x[ix0] * c + x[iy0] * s;
  x[iy0] = x[iy0] * c - x[ix0] * s;
  x[ix0] = multiModeQuad_ROS_B.temp_d;
  multiModeQuad_ROS_B.temp_d = x[iy0 + 1];
  temp_tmp = x[ix0 + 1];
  x[iy0 + 1] = multiModeQuad_ROS_B.temp_d * c - temp_tmp * s;
  x[ix0 + 1] = temp_tmp * c + multiModeQuad_ROS_B.temp_d * s;
  multiModeQuad_ROS_B.temp_d = x[iy0 + 2];
  temp_tmp = x[ix0 + 2];
  x[iy0 + 2] = multiModeQuad_ROS_B.temp_d * c - temp_tmp * s;
  x[ix0 + 2] = temp_tmp * c + multiModeQuad_ROS_B.temp_d * s;
}

// Function for MATLAB Function: '<S10>/Attitude control'
static int32_T multiModeQuad_ROS_xhseqr(real_T h[16], real_T z[16])
{
  int32_T info;
  boolean_T exitg1;
  info = 0;
  multiModeQuad_ROS_B.v[0] = 0.0;
  multiModeQuad_ROS_B.v[1] = 0.0;
  multiModeQuad_ROS_B.v[2] = 0.0;
  h[2] = 0.0;
  h[3] = 0.0;
  h[7] = 0.0;
  multiModeQuad_ROS_B.i = 3;
  exitg1 = false;
  while ((!exitg1) && (multiModeQuad_ROS_B.i + 1 >= 1)) {
    int32_T ix;
    int32_T k;
    int32_T m;
    boolean_T exitg2;
    boolean_T goto150;
    multiModeQuad_ROS_B.L = 1;
    goto150 = false;
    ix = 0;
    exitg2 = false;
    while ((!exitg2) && (ix < 301)) {
      int32_T m_tmp;
      int32_T s_tmp;
      boolean_T exitg3;
      k = multiModeQuad_ROS_B.i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > multiModeQuad_ROS_B.L)) {
        s_tmp = ((k - 1) << 2) + k;
        if (fabs(h[s_tmp]) <= 4.0083367200179456E-292) {
          exitg3 = true;
        } else {
          m_tmp = (k << 2) + k;
          multiModeQuad_ROS_B.tst = fabs(h[s_tmp - 1]) + fabs(h[m_tmp]);
          if (multiModeQuad_ROS_B.tst == 0.0) {
            if (k - 1 >= 1) {
              multiModeQuad_ROS_B.tst = fabs(h[(((k - 2) << 2) + k) - 1]);
            }

            if (k + 2 <= 4) {
              multiModeQuad_ROS_B.tst += fabs(h[m_tmp + 1]);
            }
          }

          if (fabs(h[s_tmp]) <= 2.2204460492503131E-16 * multiModeQuad_ROS_B.tst)
          {
            multiModeQuad_ROS_B.htmp1 = fabs(h[s_tmp]);
            multiModeQuad_ROS_B.tst = fabs(h[m_tmp - 1]);
            if (multiModeQuad_ROS_B.htmp1 > multiModeQuad_ROS_B.tst) {
              multiModeQuad_ROS_B.ab = multiModeQuad_ROS_B.htmp1;
              multiModeQuad_ROS_B.ba = multiModeQuad_ROS_B.tst;
            } else {
              multiModeQuad_ROS_B.ab = multiModeQuad_ROS_B.tst;
              multiModeQuad_ROS_B.ba = multiModeQuad_ROS_B.htmp1;
            }

            multiModeQuad_ROS_B.tst = h[m_tmp];
            multiModeQuad_ROS_B.htmp1 = fabs(multiModeQuad_ROS_B.tst);
            multiModeQuad_ROS_B.tst = fabs(h[s_tmp - 1] -
              multiModeQuad_ROS_B.tst);
            if (multiModeQuad_ROS_B.htmp1 > multiModeQuad_ROS_B.tst) {
              multiModeQuad_ROS_B.aa = multiModeQuad_ROS_B.htmp1;
              multiModeQuad_ROS_B.htmp1 = multiModeQuad_ROS_B.tst;
            } else {
              multiModeQuad_ROS_B.aa = multiModeQuad_ROS_B.tst;
            }

            multiModeQuad_ROS_B.tst = multiModeQuad_ROS_B.aa +
              multiModeQuad_ROS_B.ab;
            multiModeQuad_ROS_B.htmp1 = multiModeQuad_ROS_B.aa /
              multiModeQuad_ROS_B.tst * multiModeQuad_ROS_B.htmp1 *
              2.2204460492503131E-16;
            if ((4.0083367200179456E-292 >= multiModeQuad_ROS_B.htmp1) ||
                rtIsNaN(multiModeQuad_ROS_B.htmp1)) {
              multiModeQuad_ROS_B.htmp1 = 4.0083367200179456E-292;
            }

            if (multiModeQuad_ROS_B.ab / multiModeQuad_ROS_B.tst *
                multiModeQuad_ROS_B.ba <= multiModeQuad_ROS_B.htmp1) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }

      multiModeQuad_ROS_B.L = k + 1;
      if (k + 1 > 1) {
        h[k + ((k - 1) << 2)] = 0.0;
      }

      if (k + 1 >= multiModeQuad_ROS_B.i) {
        goto150 = true;
        exitg2 = true;
      } else {
        switch (ix) {
         case 10:
          s_tmp = (k << 2) + k;
          multiModeQuad_ROS_B.tst = fabs(h[(((k + 1) << 2) + k) + 2]) + fabs
            (h[s_tmp + 1]);
          multiModeQuad_ROS_B.ab = 0.75 * multiModeQuad_ROS_B.tst + h[s_tmp];
          multiModeQuad_ROS_B.h12 = -0.4375 * multiModeQuad_ROS_B.tst;
          multiModeQuad_ROS_B.aa = multiModeQuad_ROS_B.tst;
          multiModeQuad_ROS_B.htmp1 = multiModeQuad_ROS_B.ab;
          break;

         case 20:
          multiModeQuad_ROS_B.tst = fabs(h[(((multiModeQuad_ROS_B.i - 2) << 2) +
            multiModeQuad_ROS_B.i) - 1]) + fabs(h[((multiModeQuad_ROS_B.i - 1) <<
            2) + multiModeQuad_ROS_B.i]);
          multiModeQuad_ROS_B.ab = h[(multiModeQuad_ROS_B.i << 2) +
            multiModeQuad_ROS_B.i] + 0.75 * multiModeQuad_ROS_B.tst;
          multiModeQuad_ROS_B.h12 = -0.4375 * multiModeQuad_ROS_B.tst;
          multiModeQuad_ROS_B.aa = multiModeQuad_ROS_B.tst;
          multiModeQuad_ROS_B.htmp1 = multiModeQuad_ROS_B.ab;
          break;

         default:
          m = ((multiModeQuad_ROS_B.i - 1) << 2) + multiModeQuad_ROS_B.i;
          multiModeQuad_ROS_B.ab = h[m - 1];
          multiModeQuad_ROS_B.aa = h[m];
          multiModeQuad_ROS_B.h12 = h[((multiModeQuad_ROS_B.i << 2) +
            multiModeQuad_ROS_B.i) - 1];
          multiModeQuad_ROS_B.htmp1 = h[(multiModeQuad_ROS_B.i << 2) +
            multiModeQuad_ROS_B.i];
          break;
        }

        multiModeQuad_ROS_B.tst = ((fabs(multiModeQuad_ROS_B.ab) + fabs
          (multiModeQuad_ROS_B.h12)) + fabs(multiModeQuad_ROS_B.aa)) + fabs
          (multiModeQuad_ROS_B.htmp1);
        if (multiModeQuad_ROS_B.tst == 0.0) {
          multiModeQuad_ROS_B.ab = 0.0;
          multiModeQuad_ROS_B.htmp1 = 0.0;
          multiModeQuad_ROS_B.ba = 0.0;
          multiModeQuad_ROS_B.aa = 0.0;
        } else {
          multiModeQuad_ROS_B.ab /= multiModeQuad_ROS_B.tst;
          multiModeQuad_ROS_B.htmp1 /= multiModeQuad_ROS_B.tst;
          multiModeQuad_ROS_B.ba = (multiModeQuad_ROS_B.ab +
            multiModeQuad_ROS_B.htmp1) / 2.0;
          multiModeQuad_ROS_B.ab = (multiModeQuad_ROS_B.ab -
            multiModeQuad_ROS_B.ba) * (multiModeQuad_ROS_B.htmp1 -
            multiModeQuad_ROS_B.ba) - multiModeQuad_ROS_B.h12 /
            multiModeQuad_ROS_B.tst * (multiModeQuad_ROS_B.aa /
            multiModeQuad_ROS_B.tst);
          multiModeQuad_ROS_B.aa = sqrt(fabs(multiModeQuad_ROS_B.ab));
          if (multiModeQuad_ROS_B.ab >= 0.0) {
            multiModeQuad_ROS_B.ab = multiModeQuad_ROS_B.ba *
              multiModeQuad_ROS_B.tst;
            multiModeQuad_ROS_B.ba = multiModeQuad_ROS_B.ab;
            multiModeQuad_ROS_B.htmp1 = multiModeQuad_ROS_B.aa *
              multiModeQuad_ROS_B.tst;
            multiModeQuad_ROS_B.aa = -multiModeQuad_ROS_B.htmp1;
          } else {
            multiModeQuad_ROS_B.ab = multiModeQuad_ROS_B.ba +
              multiModeQuad_ROS_B.aa;
            multiModeQuad_ROS_B.ba -= multiModeQuad_ROS_B.aa;
            if (fabs(multiModeQuad_ROS_B.ab - multiModeQuad_ROS_B.htmp1) <= fabs
                (multiModeQuad_ROS_B.ba - multiModeQuad_ROS_B.htmp1)) {
              multiModeQuad_ROS_B.ab *= multiModeQuad_ROS_B.tst;
              multiModeQuad_ROS_B.ba = multiModeQuad_ROS_B.ab;
            } else {
              multiModeQuad_ROS_B.ba *= multiModeQuad_ROS_B.tst;
              multiModeQuad_ROS_B.ab = multiModeQuad_ROS_B.ba;
            }

            multiModeQuad_ROS_B.htmp1 = 0.0;
            multiModeQuad_ROS_B.aa = 0.0;
          }
        }

        m = multiModeQuad_ROS_B.i - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          s_tmp = ((m - 1) << 2) + m;
          multiModeQuad_ROS_B.h21s = h[s_tmp];
          multiModeQuad_ROS_B.h12 = h[s_tmp - 1];
          multiModeQuad_ROS_B.s_tmp = multiModeQuad_ROS_B.h12 -
            multiModeQuad_ROS_B.ba;
          multiModeQuad_ROS_B.tst = (fabs(multiModeQuad_ROS_B.s_tmp) + fabs
            (multiModeQuad_ROS_B.aa)) + fabs(multiModeQuad_ROS_B.h21s);
          multiModeQuad_ROS_B.h21s /= multiModeQuad_ROS_B.tst;
          s_tmp = (m << 2) + m;
          multiModeQuad_ROS_B.v[0] = (multiModeQuad_ROS_B.s_tmp /
            multiModeQuad_ROS_B.tst * (multiModeQuad_ROS_B.h12 -
            multiModeQuad_ROS_B.ab) + h[s_tmp - 1] * multiModeQuad_ROS_B.h21s) -
            multiModeQuad_ROS_B.aa / multiModeQuad_ROS_B.tst *
            multiModeQuad_ROS_B.htmp1;
          multiModeQuad_ROS_B.s_tmp = h[s_tmp];
          multiModeQuad_ROS_B.v[1] = (((multiModeQuad_ROS_B.h12 +
            multiModeQuad_ROS_B.s_tmp) - multiModeQuad_ROS_B.ab) -
            multiModeQuad_ROS_B.ba) * multiModeQuad_ROS_B.h21s;
          multiModeQuad_ROS_B.v[2] = h[s_tmp + 1] * multiModeQuad_ROS_B.h21s;
          multiModeQuad_ROS_B.tst = (fabs(multiModeQuad_ROS_B.v[0]) + fabs
            (multiModeQuad_ROS_B.v[1])) + fabs(multiModeQuad_ROS_B.v[2]);
          multiModeQuad_ROS_B.v[0] /= multiModeQuad_ROS_B.tst;
          multiModeQuad_ROS_B.v[1] /= multiModeQuad_ROS_B.tst;
          multiModeQuad_ROS_B.v[2] /= multiModeQuad_ROS_B.tst;
          if (k + 1 == m) {
            exitg3 = true;
          } else {
            s_tmp = ((m - 2) << 2) + m;
            if (fabs(h[s_tmp - 1]) * (fabs(multiModeQuad_ROS_B.v[1]) + fabs
                 (multiModeQuad_ROS_B.v[2])) <= ((fabs(h[s_tmp - 2]) + fabs
                  (multiModeQuad_ROS_B.h12)) + fabs(multiModeQuad_ROS_B.s_tmp)) *
                (2.2204460492503131E-16 * fabs(multiModeQuad_ROS_B.v[0]))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
        }

        for (s_tmp = m; s_tmp <= multiModeQuad_ROS_B.i; s_tmp++) {
          int32_T hoffset;
          int32_T nr;
          nr = (multiModeQuad_ROS_B.i - s_tmp) + 2;
          if (3 <= nr) {
            nr = 3;
          }

          if (s_tmp > m) {
            hoffset = ((s_tmp - 2) << 2) + s_tmp;
            for (m_tmp = 0; m_tmp < nr; m_tmp++) {
              multiModeQuad_ROS_B.v[m_tmp] = h[(m_tmp + hoffset) - 1];
            }
          }

          multiModeQuad_ROS_B.ab = multiModeQuad_ROS_B.v[0];
          multiModeQuad_ROS_B.tst = multiModeQuad_ROS_xzlarfg(nr,
            &multiModeQuad_ROS_B.ab, multiModeQuad_ROS_B.v);
          multiModeQuad_ROS_B.v[0] = multiModeQuad_ROS_B.ab;
          if (s_tmp > m) {
            h[(s_tmp + ((s_tmp - 2) << 2)) - 1] = multiModeQuad_ROS_B.ab;
            h[s_tmp + ((s_tmp - 2) << 2)] = 0.0;
            if (s_tmp < multiModeQuad_ROS_B.i) {
              h[s_tmp + 1] = 0.0;
            }
          } else if (m > k + 1) {
            h[s_tmp - 1] *= 1.0 - multiModeQuad_ROS_B.tst;
          }

          multiModeQuad_ROS_B.ab = multiModeQuad_ROS_B.v[1];
          multiModeQuad_ROS_B.ba = multiModeQuad_ROS_B.tst *
            multiModeQuad_ROS_B.v[1];
          switch (nr) {
           case 3:
            {
              int32_T sum1_tmp;
              multiModeQuad_ROS_B.aa = multiModeQuad_ROS_B.v[2];
              multiModeQuad_ROS_B.h12 = multiModeQuad_ROS_B.tst *
                multiModeQuad_ROS_B.v[2];
              for (hoffset = s_tmp - 1; hoffset + 1 < 5; hoffset++) {
                nr = (hoffset << 2) + s_tmp;
                multiModeQuad_ROS_B.htmp1 = (h[nr - 1] + h[nr] *
                  multiModeQuad_ROS_B.ab) + h[nr + 1] * multiModeQuad_ROS_B.aa;
                h[nr - 1] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.tst;
                h[nr] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.ba;
                h[nr + 1] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.h12;
              }

              if (s_tmp + 3 <= multiModeQuad_ROS_B.i + 1) {
                m_tmp = s_tmp + 3;
              } else {
                m_tmp = multiModeQuad_ROS_B.i + 1;
              }

              for (int32_T c_j = 0; c_j < m_tmp; c_j++) {
                nr = ((s_tmp - 1) << 2) + c_j;
                hoffset = (s_tmp << 2) + c_j;
                sum1_tmp = ((s_tmp + 1) << 2) + c_j;
                multiModeQuad_ROS_B.htmp1 = (h[hoffset] * multiModeQuad_ROS_B.ab
                  + h[nr]) + h[sum1_tmp] * multiModeQuad_ROS_B.aa;
                h[nr] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.tst;
                h[hoffset] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.ba;
                h[sum1_tmp] -= multiModeQuad_ROS_B.htmp1 *
                  multiModeQuad_ROS_B.h12;
              }

              for (m_tmp = 0; m_tmp < 4; m_tmp++) {
                nr = ((s_tmp - 1) << 2) + m_tmp;
                hoffset = (s_tmp << 2) + m_tmp;
                sum1_tmp = ((s_tmp + 1) << 2) + m_tmp;
                multiModeQuad_ROS_B.htmp1 = (z[hoffset] * multiModeQuad_ROS_B.ab
                  + z[nr]) + z[sum1_tmp] * multiModeQuad_ROS_B.aa;
                z[nr] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.tst;
                z[hoffset] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.ba;
                z[sum1_tmp] -= multiModeQuad_ROS_B.htmp1 *
                  multiModeQuad_ROS_B.h12;
              }
            }
            break;

           case 2:
            for (hoffset = s_tmp - 1; hoffset + 1 < 5; hoffset++) {
              nr = (hoffset << 2) + s_tmp;
              multiModeQuad_ROS_B.aa = h[nr - 1];
              multiModeQuad_ROS_B.htmp1 = h[nr] * multiModeQuad_ROS_B.ab +
                multiModeQuad_ROS_B.aa;
              h[nr - 1] = multiModeQuad_ROS_B.aa - multiModeQuad_ROS_B.htmp1 *
                multiModeQuad_ROS_B.tst;
              h[nr] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.ba;
            }

            for (m_tmp = 0; m_tmp <= multiModeQuad_ROS_B.i; m_tmp++) {
              nr = ((s_tmp - 1) << 2) + m_tmp;
              hoffset = (s_tmp << 2) + m_tmp;
              multiModeQuad_ROS_B.htmp1 = h[hoffset] * multiModeQuad_ROS_B.ab +
                h[nr];
              h[nr] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.tst;
              h[hoffset] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.ba;
            }

            for (m_tmp = 0; m_tmp < 4; m_tmp++) {
              nr = ((s_tmp - 1) << 2) + m_tmp;
              multiModeQuad_ROS_B.aa = z[nr];
              hoffset = (s_tmp << 2) + m_tmp;
              multiModeQuad_ROS_B.htmp1 = z[hoffset] * multiModeQuad_ROS_B.ab +
                multiModeQuad_ROS_B.aa;
              z[nr] = multiModeQuad_ROS_B.aa - multiModeQuad_ROS_B.htmp1 *
                multiModeQuad_ROS_B.tst;
              z[hoffset] -= multiModeQuad_ROS_B.htmp1 * multiModeQuad_ROS_B.ba;
            }
            break;
          }
        }

        ix++;
      }
    }

    if (!goto150) {
      info = multiModeQuad_ROS_B.i + 1;
      exitg1 = true;
    } else {
      if ((multiModeQuad_ROS_B.i + 1 != multiModeQuad_ROS_B.L) &&
          (multiModeQuad_ROS_B.L == multiModeQuad_ROS_B.i)) {
        ix = (multiModeQuad_ROS_B.i << 2) + multiModeQuad_ROS_B.i;
        multiModeQuad_ROS_B.ba = h[ix - 1];
        k = ((multiModeQuad_ROS_B.i - 1) << 2) + multiModeQuad_ROS_B.i;
        multiModeQuad_ROS_B.htmp1 = h[k];
        multiModeQuad_ROS_B.aa = h[ix];
        multiModeQuad_ROS_xdlanv2(&h[(multiModeQuad_ROS_B.i +
          ((multiModeQuad_ROS_B.i - 1) << 2)) - 1], &multiModeQuad_ROS_B.ba,
          &multiModeQuad_ROS_B.htmp1, &multiModeQuad_ROS_B.aa,
          &multiModeQuad_ROS_B.h12, &multiModeQuad_ROS_B.s_tmp,
          &multiModeQuad_ROS_B.h21s, &multiModeQuad_ROS_B.a__4,
          &multiModeQuad_ROS_B.tst, &multiModeQuad_ROS_B.ab);
        h[ix - 1] = multiModeQuad_ROS_B.ba;
        h[k] = multiModeQuad_ROS_B.htmp1;
        h[ix] = multiModeQuad_ROS_B.aa;
        if (4 > multiModeQuad_ROS_B.i + 1) {
          k = ((multiModeQuad_ROS_B.i + 1) << 2) + multiModeQuad_ROS_B.i;
          ix = k - 1;
          for (m = 0; m <= 2 - multiModeQuad_ROS_B.i; m++) {
            multiModeQuad_ROS_B.ba = multiModeQuad_ROS_B.tst * h[ix] +
              multiModeQuad_ROS_B.ab * h[k];
            h[k] = multiModeQuad_ROS_B.tst * h[k] - multiModeQuad_ROS_B.ab *
              h[ix];
            h[ix] = multiModeQuad_ROS_B.ba;
            k += 4;
            ix += 4;
          }
        }

        multiModeQuad_ROS_xrot(multiModeQuad_ROS_B.i - 1, h,
          ((multiModeQuad_ROS_B.i - 1) << 2) + 1, (multiModeQuad_ROS_B.i << 2) +
          1, multiModeQuad_ROS_B.tst, multiModeQuad_ROS_B.ab);
        multiModeQuad_ROS_xrot_i(z, ((multiModeQuad_ROS_B.i - 1) << 2) + 1,
          (multiModeQuad_ROS_B.i << 2) + 1, multiModeQuad_ROS_B.tst,
          multiModeQuad_ROS_B.ab);
      }

      multiModeQuad_ROS_B.i = multiModeQuad_ROS_B.L - 2;
    }
  }

  h[3] = 0.0;
  return info;
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_schur(const real_T A[16], real_T V[16], real_T T
  [16])
{
  if (multiModeQuad_ROS_anyNonFinite(A)) {
    for (multiModeQuad_ROS_B.knt = 0; multiModeQuad_ROS_B.knt < 16;
         multiModeQuad_ROS_B.knt++) {
      V[multiModeQuad_ROS_B.knt] = (rtNaN);
    }

    multiModeQuad_ROS_B.knt = 2;
    while (multiModeQuad_ROS_B.knt < 5) {
      V[multiModeQuad_ROS_B.knt - 1] = 0.0;
      multiModeQuad_ROS_B.knt++;
    }

    multiModeQuad_ROS_B.knt = 3;
    while (multiModeQuad_ROS_B.knt < 5) {
      V[multiModeQuad_ROS_B.knt + 3] = 0.0;
      multiModeQuad_ROS_B.knt++;
    }

    V[11] = 0.0;
    for (multiModeQuad_ROS_B.knt = 0; multiModeQuad_ROS_B.knt < 16;
         multiModeQuad_ROS_B.knt++) {
      T[multiModeQuad_ROS_B.knt] = (rtNaN);
    }
  } else {
    int32_T exitg1;
    boolean_T exitg2;
    memcpy(&T[0], &A[0], sizeof(real_T) << 4U);
    multiModeQuad_ROS_B.work[0] = 0.0;
    multiModeQuad_ROS_B.work[1] = 0.0;
    multiModeQuad_ROS_B.work[2] = 0.0;
    multiModeQuad_ROS_B.work[3] = 0.0;
    multiModeQuad_ROS_B.alpha1 = T[1];
    multiModeQuad_ROS_B.tau_idx_0 = 0.0;
    multiModeQuad_ROS_B.beta1 = multiModeQuad_ROS_xnrm2(2, T, 3);
    if (multiModeQuad_ROS_B.beta1 != 0.0) {
      multiModeQuad_ROS_B.beta1 = multiModeQuad_ROS_rt_hypotd_snf(T[1],
        multiModeQuad_ROS_B.beta1);
      if (T[1] >= 0.0) {
        multiModeQuad_ROS_B.beta1 = -multiModeQuad_ROS_B.beta1;
      }

      if (fabs(multiModeQuad_ROS_B.beta1) < 1.0020841800044864E-292) {
        multiModeQuad_ROS_B.knt = -1;
        do {
          multiModeQuad_ROS_B.knt++;
          multiModeQuad_ROS_B.lastc = 3;
          while (multiModeQuad_ROS_B.lastc <= 4) {
            T[multiModeQuad_ROS_B.lastc - 1] *= 9.9792015476736E+291;
            multiModeQuad_ROS_B.lastc++;
          }

          multiModeQuad_ROS_B.beta1 *= 9.9792015476736E+291;
          multiModeQuad_ROS_B.alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(multiModeQuad_ROS_B.beta1) >= 1.0020841800044864E-292));

        multiModeQuad_ROS_B.beta1 = multiModeQuad_ROS_rt_hypotd_snf
          (multiModeQuad_ROS_B.alpha1, multiModeQuad_ROS_xnrm2(2, T, 3));
        if (multiModeQuad_ROS_B.alpha1 >= 0.0) {
          multiModeQuad_ROS_B.beta1 = -multiModeQuad_ROS_B.beta1;
        }

        multiModeQuad_ROS_B.tau_idx_0 = (multiModeQuad_ROS_B.beta1 -
          multiModeQuad_ROS_B.alpha1) / multiModeQuad_ROS_B.beta1;
        multiModeQuad_ROS_B.alpha1 = 1.0 / (multiModeQuad_ROS_B.alpha1 -
          multiModeQuad_ROS_B.beta1);
        multiModeQuad_ROS_B.lastc = 3;
        while (multiModeQuad_ROS_B.lastc <= 4) {
          T[multiModeQuad_ROS_B.lastc - 1] *= multiModeQuad_ROS_B.alpha1;
          multiModeQuad_ROS_B.lastc++;
        }

        multiModeQuad_ROS_B.lastc = 0;
        while (multiModeQuad_ROS_B.lastc <= multiModeQuad_ROS_B.knt) {
          multiModeQuad_ROS_B.beta1 *= 1.0020841800044864E-292;
          multiModeQuad_ROS_B.lastc++;
        }

        multiModeQuad_ROS_B.alpha1 = multiModeQuad_ROS_B.beta1;
      } else {
        multiModeQuad_ROS_B.tau_idx_0 = (multiModeQuad_ROS_B.beta1 - T[1]) /
          multiModeQuad_ROS_B.beta1;
        multiModeQuad_ROS_B.alpha1 = 1.0 / (T[1] - multiModeQuad_ROS_B.beta1);
        multiModeQuad_ROS_B.knt = 3;
        while (multiModeQuad_ROS_B.knt <= 4) {
          T[multiModeQuad_ROS_B.knt - 1] *= multiModeQuad_ROS_B.alpha1;
          multiModeQuad_ROS_B.knt++;
        }

        multiModeQuad_ROS_B.alpha1 = multiModeQuad_ROS_B.beta1;
      }
    }

    T[1] = 1.0;
    if (multiModeQuad_ROS_B.tau_idx_0 != 0.0) {
      multiModeQuad_ROS_B.knt = 2;
      multiModeQuad_ROS_B.lastc = 3;
      while ((multiModeQuad_ROS_B.knt + 1 > 0) && (T[multiModeQuad_ROS_B.lastc] ==
              0.0)) {
        multiModeQuad_ROS_B.knt--;
        multiModeQuad_ROS_B.lastc--;
      }

      multiModeQuad_ROS_B.lastc = 4;
      exitg2 = false;
      while ((!exitg2) && (multiModeQuad_ROS_B.lastc > 0)) {
        multiModeQuad_ROS_B.ix = multiModeQuad_ROS_B.lastc + 4;
        do {
          exitg1 = 0;
          if (multiModeQuad_ROS_B.ix <= ((multiModeQuad_ROS_B.knt << 2) +
               multiModeQuad_ROS_B.lastc) + 4) {
            if (T[multiModeQuad_ROS_B.ix - 1] != 0.0) {
              exitg1 = 1;
            } else {
              multiModeQuad_ROS_B.ix += 4;
            }
          } else {
            multiModeQuad_ROS_B.lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      multiModeQuad_ROS_B.knt = -1;
      multiModeQuad_ROS_B.lastc = 0;
    }

    if (multiModeQuad_ROS_B.knt + 1 > 0) {
      if (multiModeQuad_ROS_B.lastc != 0) {
        multiModeQuad_ROS_B.ix = 0;
        while (multiModeQuad_ROS_B.ix <= multiModeQuad_ROS_B.lastc - 1) {
          multiModeQuad_ROS_B.work[multiModeQuad_ROS_B.ix] = 0.0;
          multiModeQuad_ROS_B.ix++;
        }

        multiModeQuad_ROS_B.ix = 1;
        multiModeQuad_ROS_B.jy = (multiModeQuad_ROS_B.knt << 2) + 5;
        multiModeQuad_ROS_B.iac = 5;
        while (multiModeQuad_ROS_B.iac <= multiModeQuad_ROS_B.jy) {
          multiModeQuad_ROS_B.b_ix = 0;
          multiModeQuad_ROS_B.g = (multiModeQuad_ROS_B.iac +
            multiModeQuad_ROS_B.lastc) - 1;
          multiModeQuad_ROS_B.b_ia = multiModeQuad_ROS_B.iac;
          while (multiModeQuad_ROS_B.b_ia <= multiModeQuad_ROS_B.g) {
            multiModeQuad_ROS_B.work[multiModeQuad_ROS_B.b_ix] +=
              T[multiModeQuad_ROS_B.b_ia - 1] * T[multiModeQuad_ROS_B.ix];
            multiModeQuad_ROS_B.b_ix++;
            multiModeQuad_ROS_B.b_ia++;
          }

          multiModeQuad_ROS_B.ix++;
          multiModeQuad_ROS_B.iac += 4;
        }
      }

      if (!(-multiModeQuad_ROS_B.tau_idx_0 == 0.0)) {
        multiModeQuad_ROS_B.ix = 4;
        multiModeQuad_ROS_B.jy = 1;
        multiModeQuad_ROS_B.iac = 0;
        while (multiModeQuad_ROS_B.iac <= multiModeQuad_ROS_B.knt) {
          if (T[multiModeQuad_ROS_B.jy] != 0.0) {
            multiModeQuad_ROS_B.beta1 = T[multiModeQuad_ROS_B.jy] *
              -multiModeQuad_ROS_B.tau_idx_0;
            multiModeQuad_ROS_B.b_ix = 0;
            multiModeQuad_ROS_B.g = multiModeQuad_ROS_B.ix;
            multiModeQuad_ROS_B.b_ia = multiModeQuad_ROS_B.lastc +
              multiModeQuad_ROS_B.ix;
            while (multiModeQuad_ROS_B.g + 1 <= multiModeQuad_ROS_B.b_ia) {
              T[multiModeQuad_ROS_B.g] +=
                multiModeQuad_ROS_B.work[multiModeQuad_ROS_B.b_ix] *
                multiModeQuad_ROS_B.beta1;
              multiModeQuad_ROS_B.b_ix++;
              multiModeQuad_ROS_B.g++;
            }
          }

          multiModeQuad_ROS_B.jy++;
          multiModeQuad_ROS_B.ix += 4;
          multiModeQuad_ROS_B.iac++;
        }
      }
    }

    multiModeQuad_ROS_xzlarf(3, 3, 2, multiModeQuad_ROS_B.tau_idx_0, T, 6,
      multiModeQuad_ROS_B.work);
    T[1] = multiModeQuad_ROS_B.alpha1;
    multiModeQuad_ROS_B.alpha1 = T[6];
    multiModeQuad_ROS_B.tau_idx_1 = 0.0;
    multiModeQuad_ROS_B.beta1 = multiModeQuad_ROS_xnrm2(1, T, 8);
    if (multiModeQuad_ROS_B.beta1 != 0.0) {
      multiModeQuad_ROS_B.beta1 = multiModeQuad_ROS_rt_hypotd_snf(T[6],
        multiModeQuad_ROS_B.beta1);
      if (T[6] >= 0.0) {
        multiModeQuad_ROS_B.beta1 = -multiModeQuad_ROS_B.beta1;
      }

      if (fabs(multiModeQuad_ROS_B.beta1) < 1.0020841800044864E-292) {
        multiModeQuad_ROS_B.knt = -1;
        do {
          multiModeQuad_ROS_B.knt++;
          T[7] *= 9.9792015476736E+291;
          multiModeQuad_ROS_B.beta1 *= 9.9792015476736E+291;
          multiModeQuad_ROS_B.alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(multiModeQuad_ROS_B.beta1) >= 1.0020841800044864E-292));

        multiModeQuad_ROS_B.beta1 = multiModeQuad_ROS_rt_hypotd_snf
          (multiModeQuad_ROS_B.alpha1, multiModeQuad_ROS_xnrm2(1, T, 8));
        if (multiModeQuad_ROS_B.alpha1 >= 0.0) {
          multiModeQuad_ROS_B.beta1 = -multiModeQuad_ROS_B.beta1;
        }

        multiModeQuad_ROS_B.tau_idx_1 = (multiModeQuad_ROS_B.beta1 -
          multiModeQuad_ROS_B.alpha1) / multiModeQuad_ROS_B.beta1;
        T[7] *= 1.0 / (multiModeQuad_ROS_B.alpha1 - multiModeQuad_ROS_B.beta1);
        multiModeQuad_ROS_B.lastc = 0;
        while (multiModeQuad_ROS_B.lastc <= multiModeQuad_ROS_B.knt) {
          multiModeQuad_ROS_B.beta1 *= 1.0020841800044864E-292;
          multiModeQuad_ROS_B.lastc++;
        }

        multiModeQuad_ROS_B.alpha1 = multiModeQuad_ROS_B.beta1;
      } else {
        multiModeQuad_ROS_B.tau_idx_1 = (multiModeQuad_ROS_B.beta1 - T[6]) /
          multiModeQuad_ROS_B.beta1;
        T[7] *= 1.0 / (T[6] - multiModeQuad_ROS_B.beta1);
        multiModeQuad_ROS_B.alpha1 = multiModeQuad_ROS_B.beta1;
      }
    }

    T[6] = 1.0;
    if (multiModeQuad_ROS_B.tau_idx_1 != 0.0) {
      multiModeQuad_ROS_B.knt = 1;
      multiModeQuad_ROS_B.lastc = 7;
      while ((multiModeQuad_ROS_B.knt + 1 > 0) && (T[multiModeQuad_ROS_B.lastc] ==
              0.0)) {
        multiModeQuad_ROS_B.knt--;
        multiModeQuad_ROS_B.lastc--;
      }

      multiModeQuad_ROS_B.lastc = 4;
      exitg2 = false;
      while ((!exitg2) && (multiModeQuad_ROS_B.lastc > 0)) {
        multiModeQuad_ROS_B.ix = multiModeQuad_ROS_B.lastc + 8;
        do {
          exitg1 = 0;
          if (multiModeQuad_ROS_B.ix <= ((multiModeQuad_ROS_B.knt << 2) +
               multiModeQuad_ROS_B.lastc) + 8) {
            if (T[multiModeQuad_ROS_B.ix - 1] != 0.0) {
              exitg1 = 1;
            } else {
              multiModeQuad_ROS_B.ix += 4;
            }
          } else {
            multiModeQuad_ROS_B.lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      multiModeQuad_ROS_B.knt = -1;
      multiModeQuad_ROS_B.lastc = 0;
    }

    if (multiModeQuad_ROS_B.knt + 1 > 0) {
      if (multiModeQuad_ROS_B.lastc != 0) {
        multiModeQuad_ROS_B.ix = 0;
        while (multiModeQuad_ROS_B.ix <= multiModeQuad_ROS_B.lastc - 1) {
          multiModeQuad_ROS_B.work[multiModeQuad_ROS_B.ix] = 0.0;
          multiModeQuad_ROS_B.ix++;
        }

        multiModeQuad_ROS_B.ix = 6;
        multiModeQuad_ROS_B.jy = (multiModeQuad_ROS_B.knt << 2) + 9;
        multiModeQuad_ROS_B.iac = 9;
        while (multiModeQuad_ROS_B.iac <= multiModeQuad_ROS_B.jy) {
          multiModeQuad_ROS_B.b_ix = 0;
          multiModeQuad_ROS_B.g = (multiModeQuad_ROS_B.iac +
            multiModeQuad_ROS_B.lastc) - 1;
          multiModeQuad_ROS_B.b_ia = multiModeQuad_ROS_B.iac;
          while (multiModeQuad_ROS_B.b_ia <= multiModeQuad_ROS_B.g) {
            multiModeQuad_ROS_B.work[multiModeQuad_ROS_B.b_ix] +=
              T[multiModeQuad_ROS_B.b_ia - 1] * T[multiModeQuad_ROS_B.ix];
            multiModeQuad_ROS_B.b_ix++;
            multiModeQuad_ROS_B.b_ia++;
          }

          multiModeQuad_ROS_B.ix++;
          multiModeQuad_ROS_B.iac += 4;
        }
      }

      if (!(-multiModeQuad_ROS_B.tau_idx_1 == 0.0)) {
        multiModeQuad_ROS_B.ix = 8;
        multiModeQuad_ROS_B.jy = 6;
        multiModeQuad_ROS_B.iac = 0;
        while (multiModeQuad_ROS_B.iac <= multiModeQuad_ROS_B.knt) {
          if (T[multiModeQuad_ROS_B.jy] != 0.0) {
            multiModeQuad_ROS_B.beta1 = T[multiModeQuad_ROS_B.jy] *
              -multiModeQuad_ROS_B.tau_idx_1;
            multiModeQuad_ROS_B.b_ix = 0;
            multiModeQuad_ROS_B.g = multiModeQuad_ROS_B.ix;
            multiModeQuad_ROS_B.b_ia = multiModeQuad_ROS_B.lastc +
              multiModeQuad_ROS_B.ix;
            while (multiModeQuad_ROS_B.g + 1 <= multiModeQuad_ROS_B.b_ia) {
              T[multiModeQuad_ROS_B.g] +=
                multiModeQuad_ROS_B.work[multiModeQuad_ROS_B.b_ix] *
                multiModeQuad_ROS_B.beta1;
              multiModeQuad_ROS_B.b_ix++;
              multiModeQuad_ROS_B.g++;
            }
          }

          multiModeQuad_ROS_B.jy++;
          multiModeQuad_ROS_B.ix += 4;
          multiModeQuad_ROS_B.iac++;
        }
      }
    }

    multiModeQuad_ROS_xzlarf(2, 2, 7, multiModeQuad_ROS_B.tau_idx_1, T, 11,
      multiModeQuad_ROS_B.work);
    T[6] = multiModeQuad_ROS_B.alpha1;
    multiModeQuad_ROS_B.alpha1 = T[11];
    multiModeQuad_ROS_B.tau_idx_2 = 0.0;
    multiModeQuad_ROS_B.xnorm_tmp_tmp = multiModeQuad_ROS_xnrm2(0, T, 12);
    if (multiModeQuad_ROS_B.xnorm_tmp_tmp != 0.0) {
      multiModeQuad_ROS_B.beta1 = multiModeQuad_ROS_rt_hypotd_snf(T[11],
        multiModeQuad_ROS_B.xnorm_tmp_tmp);
      if (T[11] >= 0.0) {
        multiModeQuad_ROS_B.beta1 = -multiModeQuad_ROS_B.beta1;
      }

      if (fabs(multiModeQuad_ROS_B.beta1) < 1.0020841800044864E-292) {
        multiModeQuad_ROS_B.knt = -1;
        do {
          multiModeQuad_ROS_B.knt++;
          multiModeQuad_ROS_B.beta1 *= 9.9792015476736E+291;
          multiModeQuad_ROS_B.alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(multiModeQuad_ROS_B.beta1) >= 1.0020841800044864E-292));

        multiModeQuad_ROS_B.beta1 = multiModeQuad_ROS_rt_hypotd_snf
          (multiModeQuad_ROS_B.alpha1, multiModeQuad_ROS_B.xnorm_tmp_tmp);
        if (multiModeQuad_ROS_B.alpha1 >= 0.0) {
          multiModeQuad_ROS_B.beta1 = -multiModeQuad_ROS_B.beta1;
        }

        multiModeQuad_ROS_B.tau_idx_2 = (multiModeQuad_ROS_B.beta1 -
          multiModeQuad_ROS_B.alpha1) / multiModeQuad_ROS_B.beta1;
        multiModeQuad_ROS_B.lastc = 0;
        while (multiModeQuad_ROS_B.lastc <= multiModeQuad_ROS_B.knt) {
          multiModeQuad_ROS_B.beta1 *= 1.0020841800044864E-292;
          multiModeQuad_ROS_B.lastc++;
        }

        multiModeQuad_ROS_B.alpha1 = multiModeQuad_ROS_B.beta1;
      } else {
        multiModeQuad_ROS_B.tau_idx_2 = (multiModeQuad_ROS_B.beta1 - T[11]) /
          multiModeQuad_ROS_B.beta1;
        multiModeQuad_ROS_B.alpha1 = multiModeQuad_ROS_B.beta1;
      }
    }

    T[11] = 1.0;
    if (multiModeQuad_ROS_B.tau_idx_2 != 0.0) {
      multiModeQuad_ROS_B.knt = 0;
      multiModeQuad_ROS_B.lastc = 11;
      while ((multiModeQuad_ROS_B.knt + 1 > 0) && (T[multiModeQuad_ROS_B.lastc] ==
              0.0)) {
        multiModeQuad_ROS_B.knt--;
        multiModeQuad_ROS_B.lastc--;
      }

      multiModeQuad_ROS_B.lastc = 4;
      exitg2 = false;
      while ((!exitg2) && (multiModeQuad_ROS_B.lastc > 0)) {
        multiModeQuad_ROS_B.ix = multiModeQuad_ROS_B.lastc + 12;
        do {
          exitg1 = 0;
          if (multiModeQuad_ROS_B.ix <= ((multiModeQuad_ROS_B.knt << 2) +
               multiModeQuad_ROS_B.lastc) + 12) {
            if (T[multiModeQuad_ROS_B.ix - 1] != 0.0) {
              exitg1 = 1;
            } else {
              multiModeQuad_ROS_B.ix += 4;
            }
          } else {
            multiModeQuad_ROS_B.lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      multiModeQuad_ROS_B.knt = -1;
      multiModeQuad_ROS_B.lastc = 0;
    }

    if (multiModeQuad_ROS_B.knt + 1 > 0) {
      if (multiModeQuad_ROS_B.lastc != 0) {
        multiModeQuad_ROS_B.ix = 0;
        while (multiModeQuad_ROS_B.ix <= multiModeQuad_ROS_B.lastc - 1) {
          multiModeQuad_ROS_B.work[multiModeQuad_ROS_B.ix] = 0.0;
          multiModeQuad_ROS_B.ix++;
        }

        multiModeQuad_ROS_B.ix = 11;
        multiModeQuad_ROS_B.jy = (multiModeQuad_ROS_B.knt << 2) + 13;
        multiModeQuad_ROS_B.iac = 13;
        while (multiModeQuad_ROS_B.iac <= multiModeQuad_ROS_B.jy) {
          multiModeQuad_ROS_B.b_ix = 0;
          multiModeQuad_ROS_B.g = (multiModeQuad_ROS_B.iac +
            multiModeQuad_ROS_B.lastc) - 1;
          multiModeQuad_ROS_B.b_ia = multiModeQuad_ROS_B.iac;
          while (multiModeQuad_ROS_B.b_ia <= multiModeQuad_ROS_B.g) {
            multiModeQuad_ROS_B.work[multiModeQuad_ROS_B.b_ix] +=
              T[multiModeQuad_ROS_B.b_ia - 1] * T[multiModeQuad_ROS_B.ix];
            multiModeQuad_ROS_B.b_ix++;
            multiModeQuad_ROS_B.b_ia++;
          }

          multiModeQuad_ROS_B.ix++;
          multiModeQuad_ROS_B.iac += 4;
        }
      }

      if (!(-multiModeQuad_ROS_B.tau_idx_2 == 0.0)) {
        multiModeQuad_ROS_B.ix = 12;
        multiModeQuad_ROS_B.jy = 11;
        multiModeQuad_ROS_B.iac = 0;
        while (multiModeQuad_ROS_B.iac <= multiModeQuad_ROS_B.knt) {
          if (T[multiModeQuad_ROS_B.jy] != 0.0) {
            multiModeQuad_ROS_B.beta1 = T[multiModeQuad_ROS_B.jy] *
              -multiModeQuad_ROS_B.tau_idx_2;
            multiModeQuad_ROS_B.b_ix = 0;
            multiModeQuad_ROS_B.g = multiModeQuad_ROS_B.ix;
            multiModeQuad_ROS_B.b_ia = multiModeQuad_ROS_B.lastc +
              multiModeQuad_ROS_B.ix;
            while (multiModeQuad_ROS_B.g + 1 <= multiModeQuad_ROS_B.b_ia) {
              T[multiModeQuad_ROS_B.g] +=
                multiModeQuad_ROS_B.work[multiModeQuad_ROS_B.b_ix] *
                multiModeQuad_ROS_B.beta1;
              multiModeQuad_ROS_B.b_ix++;
              multiModeQuad_ROS_B.g++;
            }
          }

          multiModeQuad_ROS_B.jy++;
          multiModeQuad_ROS_B.ix += 4;
          multiModeQuad_ROS_B.iac++;
        }
      }
    }

    multiModeQuad_ROS_xzlarf(1, 1, 12, multiModeQuad_ROS_B.tau_idx_2, T, 16,
      multiModeQuad_ROS_B.work);
    T[11] = multiModeQuad_ROS_B.alpha1;
    memcpy(&V[0], &T[0], sizeof(real_T) << 4U);
    multiModeQuad_ROS_B.knt = 0;
    while (multiModeQuad_ROS_B.knt <= 2) {
      V[multiModeQuad_ROS_B.knt + 12] = 0.0;
      multiModeQuad_ROS_B.knt++;
    }

    multiModeQuad_ROS_B.knt = 0;
    while (multiModeQuad_ROS_B.knt <= 1) {
      V[multiModeQuad_ROS_B.knt + 8] = 0.0;
      multiModeQuad_ROS_B.knt++;
    }

    multiModeQuad_ROS_B.knt = 1;
    while (multiModeQuad_ROS_B.knt + 3 < 5) {
      V[multiModeQuad_ROS_B.knt + 10] = V[multiModeQuad_ROS_B.knt + 6];
      multiModeQuad_ROS_B.knt++;
    }

    V[4] = 0.0;
    multiModeQuad_ROS_B.knt = 0;
    while (multiModeQuad_ROS_B.knt + 3 < 5) {
      V[multiModeQuad_ROS_B.knt + 6] = V[multiModeQuad_ROS_B.knt + 2];
      multiModeQuad_ROS_B.knt++;
    }

    multiModeQuad_ROS_B.work[0] = 0.0;
    V[1] = 0.0;
    multiModeQuad_ROS_B.work[1] = 0.0;
    V[2] = 0.0;
    multiModeQuad_ROS_B.work[2] = 0.0;
    V[3] = 0.0;
    multiModeQuad_ROS_B.work[3] = 0.0;
    V[0] = 1.0;
    V[15] = 1.0 - multiModeQuad_ROS_B.tau_idx_2;
    multiModeQuad_ROS_B.knt = 0;
    while (multiModeQuad_ROS_B.knt <= 1) {
      V[14 - multiModeQuad_ROS_B.knt] = 0.0;
      multiModeQuad_ROS_B.knt++;
    }

    V[10] = 1.0;
    multiModeQuad_ROS_xzlarf(2, 1, 11, multiModeQuad_ROS_B.tau_idx_1, V, 15,
      multiModeQuad_ROS_B.work);
    multiModeQuad_ROS_B.knt = 11;
    while (multiModeQuad_ROS_B.knt + 1 <= 12) {
      V[multiModeQuad_ROS_B.knt] *= -multiModeQuad_ROS_B.tau_idx_1;
      multiModeQuad_ROS_B.knt++;
    }

    V[10] = 1.0 - multiModeQuad_ROS_B.tau_idx_1;
    V[9] = 0.0;
    V[5] = 1.0;
    multiModeQuad_ROS_xzlarf(3, 2, 6, multiModeQuad_ROS_B.tau_idx_0, V, 10,
      multiModeQuad_ROS_B.work);
    multiModeQuad_ROS_B.knt = 6;
    while (multiModeQuad_ROS_B.knt + 1 <= 8) {
      V[multiModeQuad_ROS_B.knt] *= -multiModeQuad_ROS_B.tau_idx_0;
      multiModeQuad_ROS_B.knt++;
    }

    V[5] = 1.0 - multiModeQuad_ROS_B.tau_idx_0;
    multiModeQuad_ROS_xhseqr(T, V);
  }
}

// Function for MATLAB Function: '<S10>/Attitude control'
static void multiModeQuad_ROS_quat2axang(real_T q[4], real_T axang[4])
{
  real_T c_b;
  real_T theta;
  int8_T c_data[3];
  theta = 1.0 / sqrt(((q[0] * q[0] + q[1] * q[1]) + q[2] * q[2]) + q[3] * q[3]);
  q[0] *= theta;
  q[1] *= theta;
  q[2] *= theta;
  theta *= q[3];
  c_b = 1.0 / sqrt((q[1] * q[1] + q[2] * q[2]) + theta * theta);
  multiModeQuad_ROS_B.v_c[0] = q[1] * c_b;
  multiModeQuad_ROS_B.v_c[1] = q[2] * c_b;
  multiModeQuad_ROS_B.v_c[2] = theta * c_b;
  theta = 2.0 * acos(q[0]);
  if (fabs(theta) > 3.1415926535897931) {
    if (rtIsNaN(theta + 3.1415926535897931)) {
      c_b = (rtNaN);
    } else if (rtIsInf(theta + 3.1415926535897931)) {
      c_b = (rtNaN);
    } else if (theta + 3.1415926535897931 == 0.0) {
      c_b = 0.0;
    } else {
      boolean_T rEQ0;
      c_b = fmod(theta + 3.1415926535897931, 6.2831853071795862);
      rEQ0 = (c_b == 0.0);
      if (!rEQ0) {
        real_T b_q;
        b_q = fabs((theta + 3.1415926535897931) / 6.2831853071795862);
        rEQ0 = !(fabs(b_q - floor(b_q + 0.5)) > 2.2204460492503131E-16 * b_q);
      }

      if (rEQ0) {
        c_b = 0.0;
      } else if (theta + 3.1415926535897931 < 0.0) {
        c_b += 6.2831853071795862;
      }
    }

    if ((c_b == 0.0) && (theta + 3.1415926535897931 > 0.0)) {
      c_b = 6.2831853071795862;
    }

    theta = c_b - 3.1415926535897931;
  }

  if (fabs(theta) < 2.2204460492503131E-15) {
    c_data[0] = 0;
    c_data[1] = 0;
    c_data[2] = 1;
    for (int32_T i = 0; i < 3; i++) {
      multiModeQuad_ROS_B.v_c[i] = c_data[i];
    }

    theta = 0.0;
  }

  axang[0] = multiModeQuad_ROS_B.v_c[0];
  axang[1] = multiModeQuad_ROS_B.v_c[1];
  axang[2] = multiModeQuad_ROS_B.v_c[2];
  axang[3] = theta;
}

static void rt_mrdivide_U1d1x3_U2d_9vOrDY_i(const real_T u0[3], const real_T u1
  [9], real_T y[3])
{
  real_T a21;
  real_T maxval;
  int32_T r1;
  int32_T r2;
  int32_T r3;
  memcpy(&multiModeQuad_ROS_B.A_k[0], &u1[0], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = fabs(u1[0]);
  a21 = fabs(u1[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (fabs(u1[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  multiModeQuad_ROS_B.A_k[r2] = u1[r2] / u1[r1];
  multiModeQuad_ROS_B.A_k[r3] /= multiModeQuad_ROS_B.A_k[r1];
  multiModeQuad_ROS_B.A_k[r2 + 3] -= multiModeQuad_ROS_B.A_k[r1 + 3] *
    multiModeQuad_ROS_B.A_k[r2];
  multiModeQuad_ROS_B.A_k[r3 + 3] -= multiModeQuad_ROS_B.A_k[r1 + 3] *
    multiModeQuad_ROS_B.A_k[r3];
  multiModeQuad_ROS_B.A_k[r2 + 6] -= multiModeQuad_ROS_B.A_k[r1 + 6] *
    multiModeQuad_ROS_B.A_k[r2];
  multiModeQuad_ROS_B.A_k[r3 + 6] -= multiModeQuad_ROS_B.A_k[r1 + 6] *
    multiModeQuad_ROS_B.A_k[r3];
  if (fabs(multiModeQuad_ROS_B.A_k[r3 + 3]) > fabs(multiModeQuad_ROS_B.A_k[r2 +
       3])) {
    int32_T rtemp;
    rtemp = r2 + 1;
    r2 = r3;
    r3 = rtemp - 1;
  }

  multiModeQuad_ROS_B.A_k[r3 + 3] /= multiModeQuad_ROS_B.A_k[r2 + 3];
  multiModeQuad_ROS_B.A_k[r3 + 6] -= multiModeQuad_ROS_B.A_k[r3 + 3] *
    multiModeQuad_ROS_B.A_k[r2 + 6];
  y[r1] = u0[0] / multiModeQuad_ROS_B.A_k[r1];
  y[r2] = u0[1] - multiModeQuad_ROS_B.A_k[r1 + 3] * y[r1];
  y[r3] = u0[2] - multiModeQuad_ROS_B.A_k[r1 + 6] * y[r1];
  y[r2] /= multiModeQuad_ROS_B.A_k[r2 + 3];
  y[r3] -= multiModeQuad_ROS_B.A_k[r2 + 6] * y[r2];
  y[r3] /= multiModeQuad_ROS_B.A_k[r3 + 6];
  y[r2] -= multiModeQuad_ROS_B.A_k[r3 + 3] * y[r3];
  y[r1] -= y[r3] * multiModeQuad_ROS_B.A_k[r3];
  y[r1] -= y[r2] * multiModeQuad_ROS_B.A_k[r2];
}

// Model step function
void multiModeQuad_ROS_step(void)
{
  SL_Bus_multiModeQuad_ROS_std_msgs_Int16 b_varargout_2;
  boolean_T exitg2;
  boolean_T tmp;
  if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
    // set solver stop time
    rtsiSetSolverStopTime(&multiModeQuad_ROS_M->solverInfo,
                          ((multiModeQuad_ROS_M->Timing.clockTick0+1)*
      multiModeQuad_ROS_M->Timing.stepSize0));
  }                                    // end MajorTimeStep

  // Update absolute time of base rate at minor time step
  if (rtmIsMinorTimeStep(multiModeQuad_ROS_M)) {
    multiModeQuad_ROS_M->Timing.t[0] = rtsiGetT(&multiModeQuad_ROS_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
    for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 3;
         multiModeQuad_ROS_B.c_i++) {
      // Concatenate: '<S37>/Vector Concatenate' incorporates:
      //   Constant: '<S37>/Constant1'
      //   Constant: '<S37>/Constant2'
      //   Selector: '<S36>/Selector1'

      multiModeQuad_ROS_B.VectorConcatenate[6 * multiModeQuad_ROS_B.c_i] =
        multiModeQuad_ROS_P.uDOFEulerAngles2_inertia[3 * multiModeQuad_ROS_B.c_i];
      multiModeQuad_ROS_B.c_j = 6 * multiModeQuad_ROS_B.c_i + 3;
      multiModeQuad_ROS_B.VectorConcatenate[multiModeQuad_ROS_B.c_j] =
        multiModeQuad_ROS_P.Constant2_Value_n[3 * multiModeQuad_ROS_B.c_i];

      // Selector: '<S36>/Selector' incorporates:
      //   Concatenate: '<S37>/Vector Concatenate'
      //   Selector: '<S36>/Selector2'

      multiModeQuad_ROS_B.Selector_tmp = multiModeQuad_ROS_B.VectorConcatenate[6
        * multiModeQuad_ROS_B.c_i];

      // Selector: '<S36>/Selector'
      multiModeQuad_ROS_B.Selector[3 * multiModeQuad_ROS_B.c_i] =
        multiModeQuad_ROS_B.Selector_tmp;

      // Selector: '<S36>/Selector1' incorporates:
      //   Concatenate: '<S37>/Vector Concatenate'

      multiModeQuad_ROS_B.Selector1[3 * multiModeQuad_ROS_B.c_i] =
        multiModeQuad_ROS_B.VectorConcatenate[multiModeQuad_ROS_B.c_j];

      // Selector: '<S36>/Selector2'
      multiModeQuad_ROS_B.Selector2[3 * multiModeQuad_ROS_B.c_i] =
        multiModeQuad_ROS_B.Selector_tmp;

      // Concatenate: '<S37>/Vector Concatenate' incorporates:
      //   Constant: '<S37>/Constant1'
      //   Constant: '<S37>/Constant2'
      //   Selector: '<S36>/Selector'
      //   Selector: '<S36>/Selector1'
      //   Selector: '<S36>/Selector2'

      multiModeQuad_ROS_B.c_j = 3 * multiModeQuad_ROS_B.c_i + 1;
      multiModeQuad_ROS_B.sgn = 6 * multiModeQuad_ROS_B.c_i + 1;
      multiModeQuad_ROS_B.VectorConcatenate[multiModeQuad_ROS_B.sgn] =
        multiModeQuad_ROS_P.uDOFEulerAngles2_inertia[multiModeQuad_ROS_B.c_j];
      multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp = 6 *
        multiModeQuad_ROS_B.c_i + 4;
      multiModeQuad_ROS_B.VectorConcatenate[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp]
        = multiModeQuad_ROS_P.Constant2_Value_n[multiModeQuad_ROS_B.c_j];

      // Selector: '<S36>/Selector' incorporates:
      //   Concatenate: '<S37>/Vector Concatenate'
      //   Selector: '<S36>/Selector2'

      multiModeQuad_ROS_B.Selector_tmp =
        multiModeQuad_ROS_B.VectorConcatenate[multiModeQuad_ROS_B.sgn];

      // Selector: '<S36>/Selector'
      multiModeQuad_ROS_B.Selector[multiModeQuad_ROS_B.c_j] =
        multiModeQuad_ROS_B.Selector_tmp;

      // Selector: '<S36>/Selector1' incorporates:
      //   Concatenate: '<S37>/Vector Concatenate'

      multiModeQuad_ROS_B.Selector1[multiModeQuad_ROS_B.c_j] =
        multiModeQuad_ROS_B.VectorConcatenate[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp];

      // Selector: '<S36>/Selector2'
      multiModeQuad_ROS_B.Selector2[multiModeQuad_ROS_B.c_j] =
        multiModeQuad_ROS_B.Selector_tmp;

      // Concatenate: '<S37>/Vector Concatenate' incorporates:
      //   Constant: '<S37>/Constant1'
      //   Constant: '<S37>/Constant2'
      //   Selector: '<S36>/Selector'
      //   Selector: '<S36>/Selector1'
      //   Selector: '<S36>/Selector2'

      multiModeQuad_ROS_B.c_j = 3 * multiModeQuad_ROS_B.c_i + 2;
      multiModeQuad_ROS_B.sgn = 6 * multiModeQuad_ROS_B.c_i + 2;
      multiModeQuad_ROS_B.VectorConcatenate[multiModeQuad_ROS_B.sgn] =
        multiModeQuad_ROS_P.uDOFEulerAngles2_inertia[multiModeQuad_ROS_B.c_j];
      multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp = 6 *
        multiModeQuad_ROS_B.c_i + 5;
      multiModeQuad_ROS_B.VectorConcatenate[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp]
        = multiModeQuad_ROS_P.Constant2_Value_n[multiModeQuad_ROS_B.c_j];

      // Selector: '<S36>/Selector' incorporates:
      //   Concatenate: '<S37>/Vector Concatenate'
      //   Selector: '<S36>/Selector2'

      multiModeQuad_ROS_B.Selector_tmp =
        multiModeQuad_ROS_B.VectorConcatenate[multiModeQuad_ROS_B.sgn];

      // Selector: '<S36>/Selector'
      multiModeQuad_ROS_B.Selector[multiModeQuad_ROS_B.c_j] =
        multiModeQuad_ROS_B.Selector_tmp;

      // Selector: '<S36>/Selector1' incorporates:
      //   Concatenate: '<S37>/Vector Concatenate'

      multiModeQuad_ROS_B.Selector1[multiModeQuad_ROS_B.c_j] =
        multiModeQuad_ROS_B.VectorConcatenate[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp];

      // Selector: '<S36>/Selector2'
      multiModeQuad_ROS_B.Selector2[multiModeQuad_ROS_B.c_j] =
        multiModeQuad_ROS_B.Selector_tmp;
    }
  }

  // Integrator: '<S35>/phi theta psi'
  multiModeQuad_ROS_B.sincos_o1_l[0] = multiModeQuad_ROS_X.phithetapsi_CSTATE[0];

  // Trigonometry: '<S43>/sincos' incorporates:
  //   Integrator: '<S35>/phi theta psi'
  //   SignalConversion generated from: '<S43>/sincos'

  multiModeQuad_ROS_B.rtb_Sum2_idx_0 = cos
    (multiModeQuad_ROS_X.phithetapsi_CSTATE[2]);
  multiModeQuad_ROS_B.rtb_Switch_idx_0 = sin
    (multiModeQuad_ROS_X.phithetapsi_CSTATE[2]);

  // Integrator: '<S35>/phi theta psi'
  multiModeQuad_ROS_B.sincos_o1_l[1] = multiModeQuad_ROS_X.phithetapsi_CSTATE[1];

  // Trigonometry: '<S43>/sincos' incorporates:
  //   Integrator: '<S35>/phi theta psi'
  //   SignalConversion generated from: '<S43>/sincos'

  multiModeQuad_ROS_B.rtb_Sum2_idx_1 = cos
    (multiModeQuad_ROS_X.phithetapsi_CSTATE[1]);
  multiModeQuad_ROS_B.Selector_tmp = sin(multiModeQuad_ROS_X.phithetapsi_CSTATE
    [1]);

  // Integrator: '<S35>/phi theta psi'
  multiModeQuad_ROS_B.sincos_o1_l[2] = multiModeQuad_ROS_X.phithetapsi_CSTATE[2];

  // Trigonometry: '<S43>/sincos' incorporates:
  //   Integrator: '<S35>/phi theta psi'
  //   SignalConversion generated from: '<S43>/sincos'

  multiModeQuad_ROS_B.rtb_Sum2_idx_2 = cos
    (multiModeQuad_ROS_X.phithetapsi_CSTATE[0]);
  multiModeQuad_ROS_B.rtb_Switch_idx_1 = sin
    (multiModeQuad_ROS_X.phithetapsi_CSTATE[0]);

  // Fcn: '<S43>/Fcn11' incorporates:
  //   Concatenate: '<S45>/Vector Concatenate'

  multiModeQuad_ROS_B.VectorConcatenate_m[0] =
    multiModeQuad_ROS_B.rtb_Sum2_idx_0 * multiModeQuad_ROS_B.rtb_Sum2_idx_1;

  // Fcn: '<S43>/Fcn21' incorporates:
  //   Concatenate: '<S45>/Vector Concatenate'
  //   Fcn: '<S43>/Fcn22'

  multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 = multiModeQuad_ROS_B.Selector_tmp *
    multiModeQuad_ROS_B.rtb_Switch_idx_1;
  multiModeQuad_ROS_B.VectorConcatenate_m[1] =
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 * multiModeQuad_ROS_B.rtb_Sum2_idx_0
    - multiModeQuad_ROS_B.rtb_Switch_idx_0 * multiModeQuad_ROS_B.rtb_Sum2_idx_2;

  // Fcn: '<S43>/Fcn31' incorporates:
  //   Concatenate: '<S45>/Vector Concatenate'
  //   Fcn: '<S43>/Fcn32'

  multiModeQuad_ROS_B.rtb_Switch_idx_3 = multiModeQuad_ROS_B.Selector_tmp *
    multiModeQuad_ROS_B.rtb_Sum2_idx_2;
  multiModeQuad_ROS_B.VectorConcatenate_m[2] =
    multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_B.rtb_Sum2_idx_0 +
    multiModeQuad_ROS_B.rtb_Switch_idx_0 * multiModeQuad_ROS_B.rtb_Switch_idx_1;

  // Fcn: '<S43>/Fcn12' incorporates:
  //   Concatenate: '<S45>/Vector Concatenate'

  multiModeQuad_ROS_B.VectorConcatenate_m[3] =
    multiModeQuad_ROS_B.rtb_Switch_idx_0 * multiModeQuad_ROS_B.rtb_Sum2_idx_1;

  // Fcn: '<S43>/Fcn22' incorporates:
  //   Concatenate: '<S45>/Vector Concatenate'

  multiModeQuad_ROS_B.VectorConcatenate_m[4] =
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 *
    multiModeQuad_ROS_B.rtb_Switch_idx_0 + multiModeQuad_ROS_B.rtb_Sum2_idx_0 *
    multiModeQuad_ROS_B.rtb_Sum2_idx_2;

  // Fcn: '<S43>/Fcn32' incorporates:
  //   Concatenate: '<S45>/Vector Concatenate'

  multiModeQuad_ROS_B.VectorConcatenate_m[5] =
    multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_B.rtb_Switch_idx_0
    - multiModeQuad_ROS_B.rtb_Sum2_idx_0 * multiModeQuad_ROS_B.rtb_Switch_idx_1;

  // Fcn: '<S43>/Fcn13' incorporates:
  //   Concatenate: '<S45>/Vector Concatenate'

  multiModeQuad_ROS_B.VectorConcatenate_m[6] = -multiModeQuad_ROS_B.Selector_tmp;

  // Fcn: '<S43>/Fcn23' incorporates:
  //   Concatenate: '<S45>/Vector Concatenate'

  multiModeQuad_ROS_B.VectorConcatenate_m[7] =
    multiModeQuad_ROS_B.rtb_Sum2_idx_1 * multiModeQuad_ROS_B.rtb_Switch_idx_1;

  // Fcn: '<S43>/Fcn33' incorporates:
  //   Concatenate: '<S45>/Vector Concatenate'

  multiModeQuad_ROS_B.VectorConcatenate_m[8] =
    multiModeQuad_ROS_B.rtb_Sum2_idx_1 * multiModeQuad_ROS_B.rtb_Sum2_idx_2;
  if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
    // Outputs for Atomic SubSystem: '<Root>/Flight mode'
    // MATLABSystem: '<S3>/SourceBlock' incorporates:
    //   Inport: '<S13>/In1'

    tmp = Sub_multiModeQuad_ROS_472.getLatestMessage(&b_varargout_2);

    // Outputs for Enabled SubSystem: '<S3>/Enabled Subsystem' incorporates:
    //   EnablePort: '<S13>/Enable'

    if (tmp) {
      multiModeQuad_ROS_B.In1_e = b_varargout_2;
    }

    // End of MATLABSystem: '<S3>/SourceBlock'
    // End of Outputs for SubSystem: '<S3>/Enabled Subsystem'
    // End of Outputs for SubSystem: '<Root>/Flight mode'

    // Outputs for Atomic SubSystem: '<Root>/Sub setpoint attitude'
    // MATLABSystem: '<S6>/SourceBlock' incorporates:
    //   Inport: '<S14>/In1'

    tmp = Sub_multiModeQuad_ROS_496.getLatestMessage
      (&multiModeQuad_ROS_B.BusAssignment);

    // Outputs for Enabled SubSystem: '<S6>/Enabled Subsystem' incorporates:
    //   EnablePort: '<S14>/Enable'

    if (tmp) {
      multiModeQuad_ROS_B.In1 = multiModeQuad_ROS_B.BusAssignment;
    }

    // End of MATLABSystem: '<S6>/SourceBlock'
    // End of Outputs for SubSystem: '<S6>/Enabled Subsystem'
    // End of Outputs for SubSystem: '<Root>/Sub setpoint attitude'

    // Outputs for Atomic SubSystem: '<Root>/Sub setpoint velocity'
    // MATLABSystem: '<S9>/SourceBlock' incorporates:
    //   Inport: '<S17>/In1'

    tmp = Sub_multiModeQuad_ROS_426.getLatestMessage
      (&multiModeQuad_ROS_B.b_varargout_2_c);

    // Outputs for Enabled SubSystem: '<S9>/Enabled Subsystem' incorporates:
    //   EnablePort: '<S17>/Enable'

    if (tmp) {
      multiModeQuad_ROS_B.In1_h = multiModeQuad_ROS_B.b_varargout_2_c;
    }

    // End of MATLABSystem: '<S9>/SourceBlock'
    // End of Outputs for SubSystem: '<S9>/Enabled Subsystem'
    // End of Outputs for SubSystem: '<Root>/Sub setpoint velocity'
  }

  // Outputs for IfAction SubSystem: '<S55>/If Warning//Error' incorporates:
  //   ActionPort: '<S79>/if'

  for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 3;
       multiModeQuad_ROS_B.c_i++) {
    // If: '<S55>/If1' incorporates:
    //   Concatenate: '<S45>/Vector Concatenate'
    //   MATLAB Function: '<S10>/Attitude control'
    //   Math: '<S18>/Transpose'
    //   Math: '<S82>/Math Function'

    multiModeQuad_ROS_B.Product_tmp[3 * multiModeQuad_ROS_B.c_i] =
      multiModeQuad_ROS_B.VectorConcatenate_m[multiModeQuad_ROS_B.c_i];
    multiModeQuad_ROS_B.Product_tmp[3 * multiModeQuad_ROS_B.c_i + 1] =
      multiModeQuad_ROS_B.VectorConcatenate_m[multiModeQuad_ROS_B.c_i + 3];
    multiModeQuad_ROS_B.Product_tmp[3 * multiModeQuad_ROS_B.c_i + 2] =
      multiModeQuad_ROS_B.VectorConcatenate_m[multiModeQuad_ROS_B.c_i + 6];
  }

  // End of Outputs for SubSystem: '<S55>/If Warning//Error'
  for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 3;
       multiModeQuad_ROS_B.c_i++) {
    // Product: '<S42>/Product' incorporates:
    //   Integrator: '<S18>/ub,vb,wb'
    //   Math: '<S18>/Transpose'

    multiModeQuad_ROS_B.Product[multiModeQuad_ROS_B.c_i] = 0.0;
    multiModeQuad_ROS_B.Product[multiModeQuad_ROS_B.c_i] +=
      multiModeQuad_ROS_B.Product_tmp[multiModeQuad_ROS_B.c_i] *
      multiModeQuad_ROS_X.ubvbwb_CSTATE[0];
    multiModeQuad_ROS_B.Product[multiModeQuad_ROS_B.c_i] +=
      multiModeQuad_ROS_B.Product_tmp[multiModeQuad_ROS_B.c_i + 3] *
      multiModeQuad_ROS_X.ubvbwb_CSTATE[1];
    multiModeQuad_ROS_B.Product[multiModeQuad_ROS_B.c_i] +=
      multiModeQuad_ROS_B.Product_tmp[multiModeQuad_ROS_B.c_i + 6] *
      multiModeQuad_ROS_X.ubvbwb_CSTATE[2];
  }

  // Sum: '<S10>/Sum2' incorporates:
  //   Gain: '<S33>/Gain'

  multiModeQuad_ROS_B.rtb_Sum2_idx_0 = multiModeQuad_ROS_B.In1_h.Linear.X -
    multiModeQuad_ROS_B.Product[0];
  multiModeQuad_ROS_B.rtb_Sum2_idx_1 = multiModeQuad_ROS_B.In1_h.Linear.Y -
    multiModeQuad_ROS_B.Product[1];
  multiModeQuad_ROS_B.rtb_Sum2_idx_2 = multiModeQuad_ROS_B.In1_h.Linear.Z -
    multiModeQuad_ROS_P.Gain_Gain_ip * multiModeQuad_ROS_B.Product[2];

  // Gain: '<S125>/Filter Coefficient' incorporates:
  //   Gain: '<S116>/Derivative Gain'
  //   Integrator: '<S117>/Filter'
  //   Sum: '<S117>/SumD'

  multiModeQuad_ROS_B.FilterCoefficient = (multiModeQuad_ROS_P.PIDVelocityx_D *
    multiModeQuad_ROS_B.rtb_Sum2_idx_0 - multiModeQuad_ROS_X.Filter_CSTATE) *
    multiModeQuad_ROS_P.PIDVelocityx_N;

  // Gain: '<S173>/Filter Coefficient' incorporates:
  //   Gain: '<S164>/Derivative Gain'
  //   Integrator: '<S165>/Filter'
  //   Sum: '<S165>/SumD'

  multiModeQuad_ROS_B.FilterCoefficient_m = (multiModeQuad_ROS_P.PIDVelocityy_D *
    multiModeQuad_ROS_B.rtb_Sum2_idx_1 - multiModeQuad_ROS_X.Filter_CSTATE_k) *
    multiModeQuad_ROS_P.PIDVelocityy_N;

  // SignalConversion generated from: '<S31>/ SFunction ' incorporates:
  //   Gain: '<S127>/Proportional Gain'
  //   Gain: '<S175>/Proportional Gain'
  //   Integrator: '<S122>/Integrator'
  //   Integrator: '<S170>/Integrator'
  //   MATLAB Function: '<S10>/R_EW1'
  //   Sum: '<S131>/Sum'
  //   Sum: '<S179>/Sum'

  multiModeQuad_ROS_B.ubvbwb[0] = (multiModeQuad_ROS_P.PIDVelocityx_P *
    multiModeQuad_ROS_B.rtb_Sum2_idx_0 + multiModeQuad_ROS_X.Integrator_CSTATE)
    + multiModeQuad_ROS_B.FilterCoefficient;
  multiModeQuad_ROS_B.ubvbwb[1] = (multiModeQuad_ROS_P.PIDVelocityy_P *
    multiModeQuad_ROS_B.rtb_Sum2_idx_1 + multiModeQuad_ROS_X.Integrator_CSTATE_o)
    + multiModeQuad_ROS_B.FilterCoefficient_m;

  // MATLAB Function: '<S10>/R_EW1' incorporates:
  //   Concatenate: '<S45>/Vector Concatenate'
  //   SignalConversion generated from: '<S31>/ SFunction '

  for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 3;
       multiModeQuad_ROS_B.c_i++) {
    multiModeQuad_ROS_B.att_sp_Body_temp[multiModeQuad_ROS_B.c_i] =
      ((multiModeQuad_ROS_B.VectorConcatenate_m[multiModeQuad_ROS_B.c_i + 3] *
        multiModeQuad_ROS_B.ubvbwb[1] +
        multiModeQuad_ROS_B.VectorConcatenate_m[multiModeQuad_ROS_B.c_i] *
        multiModeQuad_ROS_B.ubvbwb[0]) +
       multiModeQuad_ROS_B.VectorConcatenate_m[multiModeQuad_ROS_B.c_i + 6] *
       multiModeQuad_ROS_B.In1_h.Angular.Z) * 0.017453292519943295;
  }

  // Switch: '<S10>/Switch' incorporates:
  //   Fcn: '<S32>/q0'
  //   Fcn: '<S32>/q1'
  //   Fcn: '<S32>/q2'
  //   Fcn: '<S32>/q3'
  //   Trigonometry: '<S32>/sincos'

  if (multiModeQuad_ROS_B.In1_e.Data >= multiModeQuad_ROS_P.Switch_Threshold) {
    multiModeQuad_ROS_B.rtb_Switch_idx_0 =
      multiModeQuad_ROS_B.In1.Pose.Orientation.W;
    multiModeQuad_ROS_B.rtb_Switch_idx_1 =
      multiModeQuad_ROS_B.In1.Pose.Orientation.X;
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 =
      multiModeQuad_ROS_B.In1.Pose.Orientation.Y;
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_B.In1.Pose.Orientation.Z;
  } else {
    // Gain: '<S32>/1//2' incorporates:
    //   MATLAB Function: '<S10>/R_EW1'

    multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_P.u2_Gain *
      multiModeQuad_ROS_B.att_sp_Body_temp[1];

    // Trigonometry: '<S32>/sincos'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 = sin(multiModeQuad_ROS_B.K12);
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 = cos(multiModeQuad_ROS_B.K12);

    // Fcn: '<S32>/q0'
    multiModeQuad_ROS_B.rtb_Switch_idx_0 = multiModeQuad_ROS_B.rtb_Switch_idx_3;
    multiModeQuad_ROS_B.sincos_o1[0] = multiModeQuad_ROS_B.rtb_Switch_idx_3;

    // Gain: '<S32>/1//2' incorporates:
    //   MATLAB Function: '<S10>/R_EW1'
    //   Trigonometry: '<S32>/sincos'

    multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_P.u2_Gain *
      -multiModeQuad_ROS_B.att_sp_Body_temp[0];

    // Trigonometry: '<S32>/sincos'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 = sin(multiModeQuad_ROS_B.K12);
    multiModeQuad_ROS_B.K12_f = cos(multiModeQuad_ROS_B.K12);

    // Fcn: '<S32>/q0'
    multiModeQuad_ROS_B.rtb_Switch_idx_1 = multiModeQuad_ROS_B.rtb_Switch_idx_3;
    multiModeQuad_ROS_B.sincos_o1[1] = multiModeQuad_ROS_B.rtb_Switch_idx_3;

    // Gain: '<S32>/1//2' incorporates:
    //   MATLAB Function: '<S10>/R_EW1'
    //   SignalConversion generated from: '<S31>/ SFunction '
    //   Trigonometry: '<S32>/sincos'

    multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_P.u2_Gain *
      multiModeQuad_ROS_B.In1_h.Angular.Z;

    // Trigonometry: '<S32>/sincos'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 = sin(multiModeQuad_ROS_B.K12);
    multiModeQuad_ROS_B.K12 = cos(multiModeQuad_ROS_B.K12);

    // Fcn: '<S32>/q0' incorporates:
    //   Fcn: '<S32>/q1'

    multiModeQuad_ROS_B.K14 = multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 *
      multiModeQuad_ROS_B.K12_f;
    multiModeQuad_ROS_B.rtb_Switch_idx_0 = multiModeQuad_ROS_B.rtb_Switch_idx_0 *
      multiModeQuad_ROS_B.rtb_Switch_idx_1 *
      multiModeQuad_ROS_B.rtb_Switch_idx_3 + multiModeQuad_ROS_B.K14 *
      multiModeQuad_ROS_B.K12;
    multiModeQuad_ROS_B.rtb_Switch_idx_1 = multiModeQuad_ROS_B.K14 *
      multiModeQuad_ROS_B.rtb_Switch_idx_3 - multiModeQuad_ROS_B.sincos_o1[0] *
      multiModeQuad_ROS_B.sincos_o1[1] * multiModeQuad_ROS_B.K12;

    // Fcn: '<S32>/q2' incorporates:
    //   Fcn: '<S32>/q0'
    //   Fcn: '<S32>/q1'
    //   Fcn: '<S32>/q3'

    multiModeQuad_ROS_B.K12_f *= multiModeQuad_ROS_B.sincos_o1[0];
    multiModeQuad_ROS_B.K14 = multiModeQuad_ROS_B.sincos_o1[1] *
      multiModeQuad_ROS_B.rtb_sincos_o2_idx_0;
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 = multiModeQuad_ROS_B.K14 *
      multiModeQuad_ROS_B.K12 + multiModeQuad_ROS_B.K12_f *
      multiModeQuad_ROS_B.rtb_Switch_idx_3;
    multiModeQuad_ROS_B.rtb_Switch_idx_3 = multiModeQuad_ROS_B.K12_f *
      multiModeQuad_ROS_B.K12 - multiModeQuad_ROS_B.K14 *
      multiModeQuad_ROS_B.rtb_Switch_idx_3;
  }

  // End of Switch: '<S10>/Switch'

  // MATLAB Function: '<S10>/Attitude control'
  multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_B.Product_tmp[1] +
    multiModeQuad_ROS_B.Product_tmp[3];
  multiModeQuad_ROS_B.K12_f = multiModeQuad_ROS_B.Product_tmp[2] +
    multiModeQuad_ROS_B.Product_tmp[6];
  multiModeQuad_ROS_B.K14 = multiModeQuad_ROS_B.Product_tmp[5] -
    multiModeQuad_ROS_B.Product_tmp[7];
  multiModeQuad_ROS_B.K23 = multiModeQuad_ROS_B.Product_tmp[5] +
    multiModeQuad_ROS_B.Product_tmp[7];
  multiModeQuad_ROS_B.K24 = multiModeQuad_ROS_B.Product_tmp[6] -
    multiModeQuad_ROS_B.Product_tmp[2];
  multiModeQuad_ROS_B.K34 = multiModeQuad_ROS_B.Product_tmp[1] -
    multiModeQuad_ROS_B.Product_tmp[3];
  multiModeQuad_ROS_B.A[0] = ((multiModeQuad_ROS_B.Product_tmp[0] -
    multiModeQuad_ROS_B.Product_tmp[4]) - multiModeQuad_ROS_B.Product_tmp[8]) /
    3.0;
  multiModeQuad_ROS_B.A[4] = multiModeQuad_ROS_B.K12 / 3.0;
  multiModeQuad_ROS_B.A[8] = multiModeQuad_ROS_B.K12_f / 3.0;
  multiModeQuad_ROS_B.A[12] = multiModeQuad_ROS_B.K14 / 3.0;
  multiModeQuad_ROS_B.A[1] = multiModeQuad_ROS_B.K12 / 3.0;
  multiModeQuad_ROS_B.A[5] = ((multiModeQuad_ROS_B.Product_tmp[4] -
    multiModeQuad_ROS_B.Product_tmp[0]) - multiModeQuad_ROS_B.Product_tmp[8]) /
    3.0;
  multiModeQuad_ROS_B.A[9] = multiModeQuad_ROS_B.K23 / 3.0;
  multiModeQuad_ROS_B.A[13] = multiModeQuad_ROS_B.K24 / 3.0;
  multiModeQuad_ROS_B.A[2] = multiModeQuad_ROS_B.K12_f / 3.0;
  multiModeQuad_ROS_B.A[6] = multiModeQuad_ROS_B.K23 / 3.0;
  multiModeQuad_ROS_B.A[10] = ((multiModeQuad_ROS_B.Product_tmp[8] -
    multiModeQuad_ROS_B.Product_tmp[0]) - multiModeQuad_ROS_B.Product_tmp[4]) /
    3.0;
  multiModeQuad_ROS_B.A[14] = multiModeQuad_ROS_B.K34 / 3.0;
  multiModeQuad_ROS_B.A[3] = multiModeQuad_ROS_B.K14 / 3.0;
  multiModeQuad_ROS_B.A[7] = multiModeQuad_ROS_B.K24 / 3.0;
  multiModeQuad_ROS_B.A[11] = multiModeQuad_ROS_B.K34 / 3.0;
  multiModeQuad_ROS_B.A[15] = ((multiModeQuad_ROS_B.Product_tmp[0] +
    multiModeQuad_ROS_B.Product_tmp[4]) + multiModeQuad_ROS_B.Product_tmp[8]) /
    3.0;
  if (multiModeQuad_ROS_anyNonFinite(multiModeQuad_ROS_B.A)) {
    for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 16;
         multiModeQuad_ROS_B.c_i++) {
      multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].re = (rtNaN);
      multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].im = 0.0;
    }

    multiModeQuad_ROS_B.D[0].re = (rtNaN);
    multiModeQuad_ROS_B.D[1].re = (rtNaN);
    multiModeQuad_ROS_B.D[2].re = (rtNaN);
    multiModeQuad_ROS_B.D[3].re = (rtNaN);
  } else {
    int32_T exitg1;
    tmp = true;
    multiModeQuad_ROS_B.c_j = 0;
    exitg2 = false;
    while ((!exitg2) && (multiModeQuad_ROS_B.c_j < 4)) {
      multiModeQuad_ROS_B.c_i = 0;
      do {
        exitg1 = 0;
        if (multiModeQuad_ROS_B.c_i <= multiModeQuad_ROS_B.c_j) {
          if (!(multiModeQuad_ROS_B.A[(multiModeQuad_ROS_B.c_j << 2) +
                multiModeQuad_ROS_B.c_i] == multiModeQuad_ROS_B.A
                [(multiModeQuad_ROS_B.c_i << 2) + multiModeQuad_ROS_B.c_j])) {
            tmp = false;
            exitg1 = 1;
          } else {
            multiModeQuad_ROS_B.c_i++;
          }
        } else {
          multiModeQuad_ROS_B.c_j++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }

    if (tmp) {
      multiModeQuad_ROS_schur(multiModeQuad_ROS_B.A, multiModeQuad_ROS_B.U,
        multiModeQuad_ROS_B.T);
      for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 16;
           multiModeQuad_ROS_B.c_i++) {
        multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].re =
          multiModeQuad_ROS_B.U[multiModeQuad_ROS_B.c_i];
        multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].im = 0.0;
      }

      multiModeQuad_ROS_B.D[0].re = multiModeQuad_ROS_B.T[0];
      multiModeQuad_ROS_B.D[1].re = multiModeQuad_ROS_B.T[5];
      multiModeQuad_ROS_B.D[2].re = multiModeQuad_ROS_B.T[10];
      multiModeQuad_ROS_B.D[3].re = multiModeQuad_ROS_B.T[15];
    } else {
      tmp = true;
      multiModeQuad_ROS_B.c_i = 0;
      exitg2 = false;
      while ((!exitg2) && (multiModeQuad_ROS_B.c_i < 4)) {
        multiModeQuad_ROS_B.c_j = 0;
        do {
          exitg1 = 0;
          if (multiModeQuad_ROS_B.c_j <= multiModeQuad_ROS_B.c_i) {
            if (!(multiModeQuad_ROS_B.A[(multiModeQuad_ROS_B.c_i << 2) +
                  multiModeQuad_ROS_B.c_j] == -multiModeQuad_ROS_B.A
                  [(multiModeQuad_ROS_B.c_j << 2) + multiModeQuad_ROS_B.c_i])) {
              tmp = false;
              exitg1 = 1;
            } else {
              multiModeQuad_ROS_B.c_j++;
            }
          } else {
            multiModeQuad_ROS_B.c_i++;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }

      if (tmp) {
        multiModeQuad_ROS_schur(multiModeQuad_ROS_B.A, multiModeQuad_ROS_B.U,
          multiModeQuad_ROS_B.T);
        multiModeQuad_ROS_B.c_i = 1;
        do {
          exitg1 = 0;
          if (multiModeQuad_ROS_B.c_i <= 4) {
            boolean_T guard1 = false;
            guard1 = false;
            if (multiModeQuad_ROS_B.c_i != 4) {
              multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_B.T
                [((multiModeQuad_ROS_B.c_i - 1) << 2) + multiModeQuad_ROS_B.c_i];
              if (multiModeQuad_ROS_B.K12 != 0.0) {
                multiModeQuad_ROS_B.K12 = fabs(multiModeQuad_ROS_B.K12);
                multiModeQuad_ROS_B.D[multiModeQuad_ROS_B.c_i - 1].re = 0.0;
                multiModeQuad_ROS_B.D[multiModeQuad_ROS_B.c_i - 1].im =
                  multiModeQuad_ROS_B.K12;
                multiModeQuad_ROS_B.D[multiModeQuad_ROS_B.c_i].re = 0.0;
                multiModeQuad_ROS_B.D[multiModeQuad_ROS_B.c_i].im =
                  -multiModeQuad_ROS_B.K12;
                multiModeQuad_ROS_B.c_i += 2;
              } else {
                guard1 = true;
              }
            } else {
              guard1 = true;
            }

            if (guard1) {
              multiModeQuad_ROS_B.D[multiModeQuad_ROS_B.c_i - 1].re = 0.0;
              multiModeQuad_ROS_B.D[multiModeQuad_ROS_B.c_i - 1].im = 0.0;
              multiModeQuad_ROS_B.c_i++;
            }
          } else {
            exitg1 = 1;
          }
        } while (exitg1 == 0);

        for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 16;
             multiModeQuad_ROS_B.c_i++) {
          multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].re =
            multiModeQuad_ROS_B.U[multiModeQuad_ROS_B.c_i];
          multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].im = 0.0;
        }

        multiModeQuad_ROS_B.c_j = 1;
        do {
          exitg1 = 0;
          if (multiModeQuad_ROS_B.c_j <= 4) {
            if (multiModeQuad_ROS_B.c_j != 4) {
              multiModeQuad_ROS_B.c_i = (multiModeQuad_ROS_B.c_j - 1) << 2;
              multiModeQuad_ROS_B.K12 =
                multiModeQuad_ROS_B.T[multiModeQuad_ROS_B.c_i +
                multiModeQuad_ROS_B.c_j];
              if (multiModeQuad_ROS_B.K12 != 0.0) {
                if (multiModeQuad_ROS_B.K12 < 0.0) {
                  multiModeQuad_ROS_B.sgn = 1;
                } else {
                  multiModeQuad_ROS_B.sgn = -1;
                }

                multiModeQuad_ROS_B.K12 =
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].re;
                multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp =
                  multiModeQuad_ROS_B.c_j << 2;
                multiModeQuad_ROS_B.K12_f =
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp]
                  .re * static_cast<real_T>(multiModeQuad_ROS_B.sgn);
                if (multiModeQuad_ROS_B.K12_f == 0.0) {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].re =
                    multiModeQuad_ROS_B.K12 / 1.4142135623730951;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].im = 0.0;
                } else if (multiModeQuad_ROS_B.K12 == 0.0) {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].re = 0.0;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].im =
                    multiModeQuad_ROS_B.K12_f / 1.4142135623730951;
                } else {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].re =
                    multiModeQuad_ROS_B.K12 / 1.4142135623730951;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].im =
                    multiModeQuad_ROS_B.K12_f / 1.4142135623730951;
                }

                multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp]
                  .re = multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].re;
                multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp]
                  .im = -multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].im;
                multiModeQuad_ROS_B.K12 =
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 1].re;
                multiModeQuad_ROS_B.K12_f =
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp
                  + 1].re * static_cast<real_T>(multiModeQuad_ROS_B.sgn);
                if (multiModeQuad_ROS_B.K12_f == 0.0) {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 1].re =
                    multiModeQuad_ROS_B.K12 / 1.4142135623730951;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 1].im = 0.0;
                } else if (multiModeQuad_ROS_B.K12 == 0.0) {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 1].re = 0.0;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 1].im =
                    multiModeQuad_ROS_B.K12_f / 1.4142135623730951;
                } else {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 1].re =
                    multiModeQuad_ROS_B.K12 / 1.4142135623730951;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 1].im =
                    multiModeQuad_ROS_B.K12_f / 1.4142135623730951;
                }

                multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp
                  + 1].re = multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 1].
                  re;
                multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp
                  + 1].im = -multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 1].
                  im;
                multiModeQuad_ROS_B.K12 =
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 2].re;
                multiModeQuad_ROS_B.K12_f =
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp
                  + 2].re * static_cast<real_T>(multiModeQuad_ROS_B.sgn);
                if (multiModeQuad_ROS_B.K12_f == 0.0) {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 2].re =
                    multiModeQuad_ROS_B.K12 / 1.4142135623730951;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 2].im = 0.0;
                } else if (multiModeQuad_ROS_B.K12 == 0.0) {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 2].re = 0.0;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 2].im =
                    multiModeQuad_ROS_B.K12_f / 1.4142135623730951;
                } else {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 2].re =
                    multiModeQuad_ROS_B.K12 / 1.4142135623730951;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 2].im =
                    multiModeQuad_ROS_B.K12_f / 1.4142135623730951;
                }

                multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp
                  + 2].re = multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 2].
                  re;
                multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp
                  + 2].im = -multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 2].
                  im;
                multiModeQuad_ROS_B.K12 =
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 3].re;
                multiModeQuad_ROS_B.K12_f =
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp
                  + 3].re * static_cast<real_T>(multiModeQuad_ROS_B.sgn);
                if (multiModeQuad_ROS_B.K12_f == 0.0) {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 3].re =
                    multiModeQuad_ROS_B.K12 / 1.4142135623730951;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 3].im = 0.0;
                } else if (multiModeQuad_ROS_B.K12 == 0.0) {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 3].re = 0.0;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 3].im =
                    multiModeQuad_ROS_B.K12_f / 1.4142135623730951;
                } else {
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 3].re =
                    multiModeQuad_ROS_B.K12 / 1.4142135623730951;
                  multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 3].im =
                    multiModeQuad_ROS_B.K12_f / 1.4142135623730951;
                }

                multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp
                  + 3].re = multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 3].
                  re;
                multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.rtb_VectorConcatenate_tmp
                  + 3].im = -multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 3].
                  im;
                multiModeQuad_ROS_B.c_j += 2;
              } else {
                multiModeQuad_ROS_B.c_j++;
              }
            } else {
              multiModeQuad_ROS_B.c_j++;
            }
          } else {
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      } else {
        multiModeQuad_ROS_eigStandard(multiModeQuad_ROS_B.A,
          multiModeQuad_ROS_B.V, multiModeQuad_ROS_B.D);
      }
    }
  }

  multiModeQuad_ROS_B.feedback_quat[0] = multiModeQuad_ROS_B.D[0].re;
  multiModeQuad_ROS_B.feedback_quat[1] = multiModeQuad_ROS_B.D[1].re;
  multiModeQuad_ROS_B.feedback_quat[2] = multiModeQuad_ROS_B.D[2].re;
  multiModeQuad_ROS_B.feedback_quat[3] = multiModeQuad_ROS_B.D[3].re;
  if (!rtIsNaN(multiModeQuad_ROS_B.D[0].re)) {
    multiModeQuad_ROS_B.c_i = 1;
  } else {
    multiModeQuad_ROS_B.c_i = 0;
    multiModeQuad_ROS_B.c_j = 2;
    exitg2 = false;
    while ((!exitg2) && (multiModeQuad_ROS_B.c_j < 5)) {
      if (!rtIsNaN(multiModeQuad_ROS_B.feedback_quat[multiModeQuad_ROS_B.c_j - 1]))
      {
        multiModeQuad_ROS_B.c_i = multiModeQuad_ROS_B.c_j;
        exitg2 = true;
      } else {
        multiModeQuad_ROS_B.c_j++;
      }
    }
  }

  if (multiModeQuad_ROS_B.c_i == 0) {
    multiModeQuad_ROS_B.c_j = 0;
  } else {
    multiModeQuad_ROS_B.K12 =
      multiModeQuad_ROS_B.feedback_quat[multiModeQuad_ROS_B.c_i - 1];
    multiModeQuad_ROS_B.c_j = multiModeQuad_ROS_B.c_i - 1;
    while (multiModeQuad_ROS_B.c_i + 1 < 5) {
      if (multiModeQuad_ROS_B.K12 <
          multiModeQuad_ROS_B.feedback_quat[multiModeQuad_ROS_B.c_i]) {
        multiModeQuad_ROS_B.K12 =
          multiModeQuad_ROS_B.feedback_quat[multiModeQuad_ROS_B.c_i];
        multiModeQuad_ROS_B.c_j = multiModeQuad_ROS_B.c_i;
      }

      multiModeQuad_ROS_B.c_i++;
    }
  }

  multiModeQuad_ROS_B.c_i = multiModeQuad_ROS_B.c_j << 2;
  multiModeQuad_ROS_B.feedback_quat[0] =
    multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 3].re;
  multiModeQuad_ROS_B.feedback_quat[1] =
    multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i].re;
  multiModeQuad_ROS_B.feedback_quat[2] =
    multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 1].re;
  multiModeQuad_ROS_B.feedback_quat[3] =
    multiModeQuad_ROS_B.V[multiModeQuad_ROS_B.c_i + 2].re;
  if (multiModeQuad_ROS_B.feedback_quat[0] < 0.0) {
    multiModeQuad_ROS_B.feedback_quat[0] = -multiModeQuad_ROS_B.feedback_quat[0];
    multiModeQuad_ROS_B.feedback_quat[1] = -multiModeQuad_ROS_B.feedback_quat[1];
    multiModeQuad_ROS_B.feedback_quat[2] = -multiModeQuad_ROS_B.feedback_quat[2];
    multiModeQuad_ROS_B.feedback_quat[3] = -multiModeQuad_ROS_B.feedback_quat[3];
  }

  multiModeQuad_ROS_B.K23 = ((multiModeQuad_ROS_B.feedback_quat[0] *
    multiModeQuad_ROS_B.feedback_quat[0] + multiModeQuad_ROS_B.feedback_quat[1] *
    multiModeQuad_ROS_B.feedback_quat[1]) + multiModeQuad_ROS_B.feedback_quat[2]
    * multiModeQuad_ROS_B.feedback_quat[2]) + multiModeQuad_ROS_B.feedback_quat
    [3] * multiModeQuad_ROS_B.feedback_quat[3];
  multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_B.feedback_quat[0] /
    multiModeQuad_ROS_B.K23;
  multiModeQuad_ROS_B.K12_f = -multiModeQuad_ROS_B.feedback_quat[1] /
    multiModeQuad_ROS_B.K23;
  multiModeQuad_ROS_B.K14 = -multiModeQuad_ROS_B.feedback_quat[2] /
    multiModeQuad_ROS_B.K23;
  multiModeQuad_ROS_B.K23 = -multiModeQuad_ROS_B.feedback_quat[3] /
    multiModeQuad_ROS_B.K23;
  multiModeQuad_ROS_B.quat_diff[0] = ((multiModeQuad_ROS_B.K12 *
    multiModeQuad_ROS_B.rtb_Switch_idx_0 - multiModeQuad_ROS_B.K12_f *
    multiModeQuad_ROS_B.rtb_Switch_idx_1) - multiModeQuad_ROS_B.K14 *
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0) - multiModeQuad_ROS_B.K23 *
    multiModeQuad_ROS_B.rtb_Switch_idx_3;
  multiModeQuad_ROS_B.quat_diff[1] = (multiModeQuad_ROS_B.K12 *
    multiModeQuad_ROS_B.rtb_Switch_idx_1 + multiModeQuad_ROS_B.rtb_Switch_idx_0 *
    multiModeQuad_ROS_B.K12_f) + (multiModeQuad_ROS_B.K14 *
    multiModeQuad_ROS_B.rtb_Switch_idx_3 -
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 * multiModeQuad_ROS_B.K23);
  multiModeQuad_ROS_B.quat_diff[2] = (multiModeQuad_ROS_B.K12 *
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 +
    multiModeQuad_ROS_B.rtb_Switch_idx_0 * multiModeQuad_ROS_B.K14) +
    (multiModeQuad_ROS_B.rtb_Switch_idx_1 * multiModeQuad_ROS_B.K23 -
     multiModeQuad_ROS_B.K12_f * multiModeQuad_ROS_B.rtb_Switch_idx_3);
  multiModeQuad_ROS_B.quat_diff[3] = (multiModeQuad_ROS_B.K12 *
    multiModeQuad_ROS_B.rtb_Switch_idx_3 + multiModeQuad_ROS_B.rtb_Switch_idx_0 *
    multiModeQuad_ROS_B.K23) + (multiModeQuad_ROS_B.K12_f *
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 -
    multiModeQuad_ROS_B.rtb_Switch_idx_1 * multiModeQuad_ROS_B.K14);
  multiModeQuad_ROS_quat2axang(multiModeQuad_ROS_B.quat_diff,
    multiModeQuad_ROS_B.feedback_quat);
  if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
    // Outputs for Atomic SubSystem: '<Root>/Sub setpoint rate'
    // MATLABSystem: '<S7>/SourceBlock' incorporates:
    //   Inport: '<S15>/In1'

    tmp = Sub_multiModeQuad_ROS_497.getLatestMessage
      (&multiModeQuad_ROS_B.b_varargout_2);

    // Outputs for Enabled SubSystem: '<S7>/Enabled Subsystem' incorporates:
    //   EnablePort: '<S15>/Enable'

    if (tmp) {
      multiModeQuad_ROS_B.In1_o = multiModeQuad_ROS_B.b_varargout_2;
    }

    // End of MATLABSystem: '<S7>/SourceBlock'
    // End of Outputs for SubSystem: '<S7>/Enabled Subsystem'
    // End of Outputs for SubSystem: '<Root>/Sub setpoint rate'
  }

  // Switch: '<S10>/Switch1' incorporates:
  //   MATLAB Function: '<S10>/Attitude control'

  if (multiModeQuad_ROS_B.In1_e.Data >= multiModeQuad_ROS_P.Switch1_Threshold) {
    multiModeQuad_ROS_B.att_sp_Body[0] =
      multiModeQuad_ROS_B.In1_o.Twist.Angular.X;
    multiModeQuad_ROS_B.att_sp_Body[1] =
      multiModeQuad_ROS_B.In1_o.Twist.Angular.Y;
    multiModeQuad_ROS_B.att_sp_Body[2] =
      multiModeQuad_ROS_B.In1_o.Twist.Angular.Z;
  } else {
    multiModeQuad_ROS_B.att_sp_Body[0] = multiModeQuad_ROS_B.feedback_quat[0] *
      multiModeQuad_ROS_B.feedback_quat[3];
    multiModeQuad_ROS_B.att_sp_Body[1] = multiModeQuad_ROS_B.feedback_quat[1] *
      multiModeQuad_ROS_B.feedback_quat[3];

    // Switch: '<S10>/Switch3' incorporates:
    //   MATLAB Function: '<S10>/Attitude control'

    if (multiModeQuad_ROS_B.In1_e.Data > multiModeQuad_ROS_P.Switch3_Threshold)
    {
      multiModeQuad_ROS_B.att_sp_Body[2] = multiModeQuad_ROS_B.feedback_quat[2] *
        multiModeQuad_ROS_B.feedback_quat[3];
    } else {
      multiModeQuad_ROS_B.att_sp_Body[2] = multiModeQuad_ROS_B.In1_h.Angular.Z;
    }

    // End of Switch: '<S10>/Switch3'
  }

  // End of Switch: '<S10>/Switch1'

  // TransferFcn: '<S10>/Transfer Fcn2'
  multiModeQuad_ROS_B.rtb_Switch_idx_0 = multiModeQuad_ROS_P.TransferFcn2_C *
    multiModeQuad_ROS_X.TransferFcn2_CSTATE;

  // Saturate: '<S10>/Saturation'
  if (multiModeQuad_ROS_B.rtb_Switch_idx_0 >
      multiModeQuad_ROS_P.Saturation_UpperSat) {
    multiModeQuad_ROS_B.rtb_Switch_idx_1 =
      multiModeQuad_ROS_P.Saturation_UpperSat;
  } else if (multiModeQuad_ROS_B.rtb_Switch_idx_0 <
             multiModeQuad_ROS_P.Saturation_LowerSat) {
    multiModeQuad_ROS_B.rtb_Switch_idx_1 =
      multiModeQuad_ROS_P.Saturation_LowerSat;
  } else {
    multiModeQuad_ROS_B.rtb_Switch_idx_1 = multiModeQuad_ROS_B.rtb_Switch_idx_0;
  }

  // End of Saturate: '<S10>/Saturation'

  // TransferFcn: '<S10>/Transfer Fcn1'
  multiModeQuad_ROS_B.rtb_Switch_idx_0 = multiModeQuad_ROS_P.TransferFcn1_C *
    multiModeQuad_ROS_X.TransferFcn1_CSTATE;

  // Saturate: '<S10>/Saturation1'
  if (multiModeQuad_ROS_B.rtb_Switch_idx_0 >
      multiModeQuad_ROS_P.Saturation1_UpperSat) {
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 =
      multiModeQuad_ROS_P.Saturation1_UpperSat;
  } else if (multiModeQuad_ROS_B.rtb_Switch_idx_0 <
             multiModeQuad_ROS_P.Saturation1_LowerSat) {
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 =
      multiModeQuad_ROS_P.Saturation1_LowerSat;
  } else {
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 =
      multiModeQuad_ROS_B.rtb_Switch_idx_0;
  }

  // End of Saturate: '<S10>/Saturation1'

  // TransferFcn: '<S10>/Transfer Fcn3'
  multiModeQuad_ROS_B.rtb_Switch_idx_0 = multiModeQuad_ROS_P.TransferFcn3_C *
    multiModeQuad_ROS_X.TransferFcn3_CSTATE;

  // Saturate: '<S10>/Saturation2'
  if (multiModeQuad_ROS_B.rtb_Switch_idx_0 >
      multiModeQuad_ROS_P.Saturation2_UpperSat) {
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_P.Saturation2_UpperSat;
  } else if (multiModeQuad_ROS_B.rtb_Switch_idx_0 <
             multiModeQuad_ROS_P.Saturation2_LowerSat) {
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_P.Saturation2_LowerSat;
  } else {
    multiModeQuad_ROS_B.rtb_Switch_idx_3 = multiModeQuad_ROS_B.rtb_Switch_idx_0;
  }

  // End of Saturate: '<S10>/Saturation2'

  // TransferFcn: '<S10>/Transfer Fcn4'
  multiModeQuad_ROS_B.rtb_Switch_idx_0 = multiModeQuad_ROS_P.TransferFcn4_C *
    multiModeQuad_ROS_X.TransferFcn4_CSTATE;

  // Saturate: '<S10>/Saturation3'
  if (multiModeQuad_ROS_B.rtb_Switch_idx_0 >
      multiModeQuad_ROS_P.Saturation3_UpperSat) {
    multiModeQuad_ROS_B.rtb_Switch_idx_0 =
      multiModeQuad_ROS_P.Saturation3_UpperSat;
  } else if (multiModeQuad_ROS_B.rtb_Switch_idx_0 <
             multiModeQuad_ROS_P.Saturation3_LowerSat) {
    multiModeQuad_ROS_B.rtb_Switch_idx_0 =
      multiModeQuad_ROS_P.Saturation3_LowerSat;
  }

  // End of Saturate: '<S10>/Saturation3'
  for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 3;
       multiModeQuad_ROS_B.c_i++) {
    // Sum: '<S10>/Sum7' incorporates:
    //   Integrator: '<S18>/p,q,r '

    multiModeQuad_ROS_B.sincos_o1[multiModeQuad_ROS_B.c_i] =
      multiModeQuad_ROS_B.att_sp_Body[multiModeQuad_ROS_B.c_i] -
      multiModeQuad_ROS_X.pqr_CSTATE[multiModeQuad_ROS_B.c_i];

    // Product: '<S47>/Product' incorporates:
    //   Integrator: '<S18>/p,q,r '
    //   Selector: '<S36>/Selector'
    //   Sum: '<Root>/Sum'

    multiModeQuad_ROS_B.att_sp_Body[multiModeQuad_ROS_B.c_i] =
      (multiModeQuad_ROS_B.Selector[multiModeQuad_ROS_B.c_i + 3] *
       multiModeQuad_ROS_X.pqr_CSTATE[1] +
       multiModeQuad_ROS_B.Selector[multiModeQuad_ROS_B.c_i] *
       multiModeQuad_ROS_X.pqr_CSTATE[0]) +
      multiModeQuad_ROS_B.Selector[multiModeQuad_ROS_B.c_i + 6] *
      multiModeQuad_ROS_X.pqr_CSTATE[2];

    // Product: '<S48>/Product' incorporates:
    //   Integrator: '<S18>/p,q,r '
    //   Selector: '<S36>/Selector1'

    multiModeQuad_ROS_B.att_sp_Body_temp[multiModeQuad_ROS_B.c_i] =
      (multiModeQuad_ROS_B.Selector1[multiModeQuad_ROS_B.c_i + 3] *
       multiModeQuad_ROS_X.pqr_CSTATE[1] +
       multiModeQuad_ROS_B.Selector1[multiModeQuad_ROS_B.c_i] *
       multiModeQuad_ROS_X.pqr_CSTATE[0]) +
      multiModeQuad_ROS_B.Selector1[multiModeQuad_ROS_B.c_i + 6] *
      multiModeQuad_ROS_X.pqr_CSTATE[2];
  }

  // Sum: '<S36>/Sum2' incorporates:
  //   Integrator: '<S18>/p,q,r '
  //   MATLAB Function: '<S10>/Mapping'
  //   Product: '<S48>/Product'
  //   Product: '<S49>/i x j'
  //   Product: '<S49>/j x k'
  //   Product: '<S49>/k x i'
  //   Product: '<S50>/i x k'
  //   Product: '<S50>/j x i'
  //   Product: '<S50>/k x j'
  //   Sum: '<S46>/Sum'

  multiModeQuad_ROS_B.rtb_Saturation2_p[0] =
    ((((multiModeQuad_ROS_B.rtb_Switch_idx_3 +
        multiModeQuad_ROS_B.rtb_Switch_idx_0) -
       multiModeQuad_ROS_B.rtb_Switch_idx_1) -
      multiModeQuad_ROS_B.rtb_sincos_o2_idx_0) * 0.15 -
     multiModeQuad_ROS_B.att_sp_Body_temp[0]) - (multiModeQuad_ROS_X.pqr_CSTATE
    [1] * multiModeQuad_ROS_B.att_sp_Body[2] - multiModeQuad_ROS_B.att_sp_Body[1]
    * multiModeQuad_ROS_X.pqr_CSTATE[2]);
  multiModeQuad_ROS_B.rtb_Saturation2_p[1] =
    ((((multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 +
        multiModeQuad_ROS_B.rtb_Switch_idx_0) -
       multiModeQuad_ROS_B.rtb_Switch_idx_1) -
      multiModeQuad_ROS_B.rtb_Switch_idx_3) * 0.15 -
     multiModeQuad_ROS_B.att_sp_Body_temp[1]) -
    (multiModeQuad_ROS_B.att_sp_Body[0] * multiModeQuad_ROS_X.pqr_CSTATE[2] -
     multiModeQuad_ROS_X.pqr_CSTATE[0] * multiModeQuad_ROS_B.att_sp_Body[2]);
  multiModeQuad_ROS_B.rtb_Saturation2_p[2] =
    ((((multiModeQuad_ROS_B.rtb_sincos_o2_idx_0 +
        multiModeQuad_ROS_B.rtb_Switch_idx_3) -
       multiModeQuad_ROS_B.rtb_Switch_idx_1) -
      multiModeQuad_ROS_B.rtb_Switch_idx_0) * 0.2 * 0.2 -
     multiModeQuad_ROS_B.att_sp_Body_temp[2]) - (multiModeQuad_ROS_X.pqr_CSTATE
    [0] * multiModeQuad_ROS_B.att_sp_Body[1] - multiModeQuad_ROS_B.att_sp_Body[0]
    * multiModeQuad_ROS_X.pqr_CSTATE[1]);

  // Product: '<S36>/Product2' incorporates:
  //   Selector: '<S36>/Selector2'

  rt_mrdivide_U1d1x3_U2d_9vOrDY_i(multiModeQuad_ROS_B.rtb_Saturation2_p,
    multiModeQuad_ROS_B.Selector2, multiModeQuad_ROS_B.Product2);
  if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
    // Reshape: '<S21>/Reshape' incorporates:
    //   Constant: '<S21>/Constant4'
    //   Constant: '<S21>/Constant5'
    //   Constant: '<S21>/Mass'
    //   Gain: '<S21>/g'

    multiModeQuad_ROS_B.Reshape[0] = multiModeQuad_ROS_P.Constant4_Value;
    multiModeQuad_ROS_B.Reshape[1] = multiModeQuad_ROS_P.Constant5_Value;
    multiModeQuad_ROS_B.Reshape[2] = multiModeQuad_ROS_P.g_Gain *
      multiModeQuad_ROS_P.Mass_Value;
  }

  // Sum: '<S38>/Sum' incorporates:
  //   Integrator: '<S18>/p,q,r '
  //   Integrator: '<S18>/ub,vb,wb'
  //   Product: '<S51>/i x j'
  //   Product: '<S51>/j x k'
  //   Product: '<S51>/k x i'
  //   Product: '<S52>/i x k'
  //   Product: '<S52>/j x i'
  //   Product: '<S52>/k x j'

  multiModeQuad_ROS_B.ubvbwb[0] = multiModeQuad_ROS_X.ubvbwb_CSTATE[1] *
    multiModeQuad_ROS_X.pqr_CSTATE[2] - multiModeQuad_ROS_X.pqr_CSTATE[1] *
    multiModeQuad_ROS_X.ubvbwb_CSTATE[2];
  multiModeQuad_ROS_B.ubvbwb[1] = multiModeQuad_ROS_X.pqr_CSTATE[0] *
    multiModeQuad_ROS_X.ubvbwb_CSTATE[2] - multiModeQuad_ROS_X.ubvbwb_CSTATE[0] *
    multiModeQuad_ROS_X.pqr_CSTATE[2];
  multiModeQuad_ROS_B.ubvbwb[2] = multiModeQuad_ROS_X.ubvbwb_CSTATE[0] *
    multiModeQuad_ROS_X.pqr_CSTATE[1] - multiModeQuad_ROS_X.pqr_CSTATE[0] *
    multiModeQuad_ROS_X.ubvbwb_CSTATE[1];

  // Sum: '<S10>/Sum15' incorporates:
  //   MATLAB Function: '<S10>/Mapping'

  multiModeQuad_ROS_B.att_sp_Body_temp[0] = 0.0;
  multiModeQuad_ROS_B.att_sp_Body_temp[1] = 0.0;
  multiModeQuad_ROS_B.att_sp_Body_temp[2] =
    -(((multiModeQuad_ROS_B.rtb_Switch_idx_1 +
        multiModeQuad_ROS_B.rtb_sincos_o2_idx_0) +
       multiModeQuad_ROS_B.rtb_Switch_idx_3) +
      multiModeQuad_ROS_B.rtb_Switch_idx_0);
  for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 3;
       multiModeQuad_ROS_B.c_i++) {
    // Trigonometry: '<S44>/sincos'
    multiModeQuad_ROS_B.rtb_Switch_idx_1 =
      multiModeQuad_ROS_B.sincos_o1_l[multiModeQuad_ROS_B.c_i];

    // Sum: '<S18>/Sum' incorporates:
    //   Concatenate: '<S45>/Vector Concatenate'
    //   Constant: '<S37>/Constant'
    //   Product: '<S18>/Product'
    //   Product: '<S21>/Product'
    //   Reshape: '<S21>/Reshape'
    //   Sum: '<S10>/Sum15'

    multiModeQuad_ROS_B.Sum[multiModeQuad_ROS_B.c_i] =
      (((multiModeQuad_ROS_B.VectorConcatenate_m[multiModeQuad_ROS_B.c_i + 3] *
         multiModeQuad_ROS_B.Reshape[1] +
         multiModeQuad_ROS_B.VectorConcatenate_m[multiModeQuad_ROS_B.c_i] *
         multiModeQuad_ROS_B.Reshape[0]) +
        multiModeQuad_ROS_B.VectorConcatenate_m[multiModeQuad_ROS_B.c_i + 6] *
        multiModeQuad_ROS_B.Reshape[2]) +
       multiModeQuad_ROS_B.att_sp_Body_temp[multiModeQuad_ROS_B.c_i]) /
      multiModeQuad_ROS_P.uDOFEulerAngles2_mass_0 +
      multiModeQuad_ROS_B.ubvbwb[multiModeQuad_ROS_B.c_i];

    // Trigonometry: '<S44>/sincos'
    multiModeQuad_ROS_B.sincos_o1_l[multiModeQuad_ROS_B.c_i] = sin
      (multiModeQuad_ROS_B.rtb_Switch_idx_1);
    multiModeQuad_ROS_B.ubvbwb[multiModeQuad_ROS_B.c_i] = cos
      (multiModeQuad_ROS_B.rtb_Switch_idx_1);
  }

  // MATLABSystem: '<Root>/Get Parameter'
  ParamGet_multiModeQuad_ROS_459.get_parameter
    (&multiModeQuad_ROS_B.rtb_Switch_idx_0);

  // MATLABSystem: '<Root>/Get Parameter1'
  ParamGet_multiModeQuad_ROS_460.get_parameter
    (&multiModeQuad_ROS_B.rtb_Switch_idx_1);

  // MATLABSystem: '<Root>/Get Parameter2'
  ParamGet_multiModeQuad_ROS_461.get_parameter
    (&multiModeQuad_ROS_B.rtb_sincos_o2_idx_0);

  // Fcn: '<S44>/phidot' incorporates:
  //   Fcn: '<S44>/psidot'
  //   Integrator: '<S18>/p,q,r '

  multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_B.sincos_o1_l[0] *
    multiModeQuad_ROS_X.pqr_CSTATE[1] + multiModeQuad_ROS_B.ubvbwb[0] *
    multiModeQuad_ROS_X.pqr_CSTATE[2];

  // Saturate: '<S10>/Saturation7' incorporates:
  //   Fcn: '<S44>/phidot'
  //   Integrator: '<S18>/p,q,r '

  multiModeQuad_ROS_B.rtb_Switch_idx_3 = multiModeQuad_ROS_B.sincos_o1_l[1] /
    multiModeQuad_ROS_B.ubvbwb[1] * multiModeQuad_ROS_B.K12 +
    multiModeQuad_ROS_X.pqr_CSTATE[0];

  // Integrator: '<S213>/Filter' incorporates:
  //   Fcn: '<S44>/thetadot'
  //   Integrator: '<S18>/p,q,r '

  multiModeQuad_ROS_B.K12_f = multiModeQuad_ROS_B.ubvbwb[0] *
    multiModeQuad_ROS_X.pqr_CSTATE[1] - multiModeQuad_ROS_B.sincos_o1_l[0] *
    multiModeQuad_ROS_X.pqr_CSTATE[2];

  // SignalConversion generated from: '<S35>/phi theta psi' incorporates:
  //   Fcn: '<S44>/psidot'

  multiModeQuad_ROS_B.TmpSignalConversionAtphithetaps[0] =
    multiModeQuad_ROS_B.rtb_Switch_idx_3;
  multiModeQuad_ROS_B.TmpSignalConversionAtphithetaps[1] =
    multiModeQuad_ROS_B.K12_f;
  multiModeQuad_ROS_B.TmpSignalConversionAtphithetaps[2] =
    multiModeQuad_ROS_B.K12 / multiModeQuad_ROS_B.ubvbwb[1];

  // Sum: '<S56>/Add'
  multiModeQuad_ROS_B.K12 = (multiModeQuad_ROS_B.VectorConcatenate_m[0] +
    multiModeQuad_ROS_B.VectorConcatenate_m[4]) +
    multiModeQuad_ROS_B.VectorConcatenate_m[8];

  // If: '<S20>/If' incorporates:
  //   Sum: '<S56>/Add'

  if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
    multiModeQuad_ROS_DW.If_ActiveSubsystem = static_cast<int8_T>
      (!(multiModeQuad_ROS_B.K12 > 0.0));
  }

  switch (multiModeQuad_ROS_DW.If_ActiveSubsystem) {
   case 0:
    // Outputs for IfAction SubSystem: '<S20>/Positive Trace' incorporates:
    //   ActionPort: '<S54>/Action Port'

    // Sqrt: '<S54>/sqrt' incorporates:
    //   Constant: '<S54>/Constant'
    //   Sum: '<S54>/Sum'
    //   Sum: '<S56>/Add'

    multiModeQuad_ROS_B.rtb_Switch_idx_3 = sqrt(multiModeQuad_ROS_B.K12 +
      multiModeQuad_ROS_P.Constant_Value_o);

    // Gain: '<S54>/Gain' incorporates:
    //   Merge: '<S20>/Merge'

    multiModeQuad_ROS_B.Merge[0] = multiModeQuad_ROS_P.Gain_Gain *
      multiModeQuad_ROS_B.rtb_Switch_idx_3;

    // Gain: '<S54>/Gain1'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 *= multiModeQuad_ROS_P.Gain1_Gain;

    // Product: '<S54>/Product' incorporates:
    //   Fcn: '<S43>/Fcn13'
    //   Merge: '<S20>/Merge'
    //   Sum: '<S76>/Add'
    //   Sum: '<S77>/Add'
    //   Sum: '<S78>/Add'

    multiModeQuad_ROS_B.Merge[1] = (multiModeQuad_ROS_B.VectorConcatenate_m[7] -
      multiModeQuad_ROS_B.VectorConcatenate_m[5]) /
      multiModeQuad_ROS_B.rtb_Switch_idx_3;
    multiModeQuad_ROS_B.Merge[2] = (multiModeQuad_ROS_B.VectorConcatenate_m[2] -
      (-multiModeQuad_ROS_B.Selector_tmp)) /
      multiModeQuad_ROS_B.rtb_Switch_idx_3;
    multiModeQuad_ROS_B.Merge[3] = (multiModeQuad_ROS_B.VectorConcatenate_m[3] -
      multiModeQuad_ROS_B.VectorConcatenate_m[1]) /
      multiModeQuad_ROS_B.rtb_Switch_idx_3;

    // End of Outputs for SubSystem: '<S20>/Positive Trace'
    break;

   case 1:
    // Outputs for IfAction SubSystem: '<S20>/Negative Trace' incorporates:
    //   ActionPort: '<S53>/Action Port'

    // If: '<S53>/Find Maximum Diagonal Value'
    if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
      if ((multiModeQuad_ROS_B.VectorConcatenate_m[4] >
           multiModeQuad_ROS_B.VectorConcatenate_m[0]) &&
          (multiModeQuad_ROS_B.VectorConcatenate_m[4] >
           multiModeQuad_ROS_B.VectorConcatenate_m[8])) {
        multiModeQuad_ROS_DW.FindMaximumDiagonalValue_Active = 0;
      } else if (multiModeQuad_ROS_B.VectorConcatenate_m[8] >
                 multiModeQuad_ROS_B.VectorConcatenate_m[0]) {
        multiModeQuad_ROS_DW.FindMaximumDiagonalValue_Active = 1;
      } else {
        multiModeQuad_ROS_DW.FindMaximumDiagonalValue_Active = 2;
      }
    }

    switch (multiModeQuad_ROS_DW.FindMaximumDiagonalValue_Active) {
     case 0:
      // Outputs for IfAction SubSystem: '<S53>/Maximum Value at DCM(2,2)' incorporates:
      //   ActionPort: '<S58>/Action Port'

      // Sqrt: '<S58>/sqrt' incorporates:
      //   Constant: '<S70>/Constant'
      //   Sum: '<S70>/Add'

      multiModeQuad_ROS_B.K12_f = sqrt
        (((multiModeQuad_ROS_B.VectorConcatenate_m[4] -
           multiModeQuad_ROS_B.VectorConcatenate_m[0]) -
          multiModeQuad_ROS_B.VectorConcatenate_m[8]) +
         multiModeQuad_ROS_P.Constant_Value_n);

      // Switch: '<S69>/Switch' incorporates:
      //   Constant: '<S69>/Constant1'
      //   Constant: '<S69>/Constant2'

      if (multiModeQuad_ROS_B.K12_f != 0.0) {
        multiModeQuad_ROS_B.rtb_Switch_idx_3 =
          multiModeQuad_ROS_P.Constant1_Value;
        multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_B.K12_f;
      } else {
        multiModeQuad_ROS_B.rtb_Switch_idx_3 =
          multiModeQuad_ROS_P.Constant2_Value[0];
        multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_P.Constant2_Value[1];
      }

      // End of Switch: '<S69>/Switch'

      // Product: '<S69>/Product'
      multiModeQuad_ROS_B.rtb_Switch_idx_3 /= multiModeQuad_ROS_B.K12;

      // Gain: '<S58>/Gain1' incorporates:
      //   Merge: '<S20>/Merge'
      //   Product: '<S58>/Product'
      //   Sum: '<S68>/Add'

      multiModeQuad_ROS_B.Merge[1] = (multiModeQuad_ROS_B.VectorConcatenate_m[1]
        + multiModeQuad_ROS_B.VectorConcatenate_m[3]) *
        multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Gain1_Gain_k;

      // Gain: '<S58>/Gain3' incorporates:
      //   Merge: '<S20>/Merge'
      //   Product: '<S58>/Product'
      //   Sum: '<S67>/Add'

      multiModeQuad_ROS_B.Merge[3] = (multiModeQuad_ROS_B.VectorConcatenate_m[5]
        + multiModeQuad_ROS_B.VectorConcatenate_m[7]) *
        multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Gain3_Gain;

      // Gain: '<S58>/Gain4' incorporates:
      //   Fcn: '<S43>/Fcn13'
      //   Merge: '<S20>/Merge'
      //   Product: '<S58>/Product'
      //   Sum: '<S66>/Add'

      multiModeQuad_ROS_B.Merge[0] = (multiModeQuad_ROS_B.VectorConcatenate_m[2]
        - (-multiModeQuad_ROS_B.Selector_tmp)) *
        multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Gain4_Gain;

      // Gain: '<S58>/Gain' incorporates:
      //   Merge: '<S20>/Merge'

      multiModeQuad_ROS_B.Merge[2] = multiModeQuad_ROS_P.Gain_Gain_g *
        multiModeQuad_ROS_B.K12_f;

      // End of Outputs for SubSystem: '<S53>/Maximum Value at DCM(2,2)'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S53>/Maximum Value at DCM(3,3)' incorporates:
      //   ActionPort: '<S59>/Action Port'

      // Sqrt: '<S59>/sqrt' incorporates:
      //   Constant: '<S75>/Constant'
      //   Sum: '<S75>/Add'

      multiModeQuad_ROS_B.K12_f = sqrt
        (((multiModeQuad_ROS_B.VectorConcatenate_m[8] -
           multiModeQuad_ROS_B.VectorConcatenate_m[0]) -
          multiModeQuad_ROS_B.VectorConcatenate_m[4]) +
         multiModeQuad_ROS_P.Constant_Value_pu);

      // Switch: '<S74>/Switch' incorporates:
      //   Constant: '<S74>/Constant1'
      //   Constant: '<S74>/Constant2'

      if (multiModeQuad_ROS_B.K12_f != 0.0) {
        multiModeQuad_ROS_B.rtb_Switch_idx_3 =
          multiModeQuad_ROS_P.Constant1_Value_h;
        multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_B.K12_f;
      } else {
        multiModeQuad_ROS_B.rtb_Switch_idx_3 =
          multiModeQuad_ROS_P.Constant2_Value_k[0];
        multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_P.Constant2_Value_k[1];
      }

      // End of Switch: '<S74>/Switch'

      // Product: '<S74>/Product'
      multiModeQuad_ROS_B.rtb_Switch_idx_3 /= multiModeQuad_ROS_B.K12;

      // Gain: '<S59>/Gain1' incorporates:
      //   Fcn: '<S43>/Fcn13'
      //   Merge: '<S20>/Merge'
      //   Product: '<S59>/Product'
      //   Sum: '<S71>/Add'

      multiModeQuad_ROS_B.Merge[1] = (multiModeQuad_ROS_B.VectorConcatenate_m[2]
        + -multiModeQuad_ROS_B.Selector_tmp) *
        multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Gain1_Gain_l;

      // Gain: '<S59>/Gain2' incorporates:
      //   Merge: '<S20>/Merge'
      //   Product: '<S59>/Product'
      //   Sum: '<S72>/Add'

      multiModeQuad_ROS_B.Merge[2] = (multiModeQuad_ROS_B.VectorConcatenate_m[5]
        + multiModeQuad_ROS_B.VectorConcatenate_m[7]) *
        multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Gain2_Gain;

      // Gain: '<S59>/Gain3' incorporates:
      //   Merge: '<S20>/Merge'
      //   Product: '<S59>/Product'
      //   Sum: '<S73>/Add'

      multiModeQuad_ROS_B.Merge[0] = (multiModeQuad_ROS_B.VectorConcatenate_m[3]
        - multiModeQuad_ROS_B.VectorConcatenate_m[1]) *
        multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Gain3_Gain_a;

      // Gain: '<S59>/Gain' incorporates:
      //   Merge: '<S20>/Merge'

      multiModeQuad_ROS_B.Merge[3] = multiModeQuad_ROS_P.Gain_Gain_d *
        multiModeQuad_ROS_B.K12_f;

      // End of Outputs for SubSystem: '<S53>/Maximum Value at DCM(3,3)'
      break;

     case 2:
      // Outputs for IfAction SubSystem: '<S53>/Maximum Value at DCM(1,1)' incorporates:
      //   ActionPort: '<S57>/Action Port'

      // Sqrt: '<S57>/sqrt' incorporates:
      //   Constant: '<S65>/Constant'
      //   Sum: '<S65>/Add'

      multiModeQuad_ROS_B.K12_f = sqrt
        (((multiModeQuad_ROS_B.VectorConcatenate_m[0] -
           multiModeQuad_ROS_B.VectorConcatenate_m[4]) -
          multiModeQuad_ROS_B.VectorConcatenate_m[8]) +
         multiModeQuad_ROS_P.Constant_Value_fv);

      // Switch: '<S64>/Switch' incorporates:
      //   Constant: '<S64>/Constant1'
      //   Constant: '<S64>/Constant2'

      if (multiModeQuad_ROS_B.K12_f != 0.0) {
        multiModeQuad_ROS_B.rtb_Switch_idx_3 =
          multiModeQuad_ROS_P.Constant1_Value_b;
        multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_B.K12_f;
      } else {
        multiModeQuad_ROS_B.rtb_Switch_idx_3 =
          multiModeQuad_ROS_P.Constant2_Value_o[0];
        multiModeQuad_ROS_B.K12 = multiModeQuad_ROS_P.Constant2_Value_o[1];
      }

      // End of Switch: '<S64>/Switch'

      // Product: '<S64>/Product'
      multiModeQuad_ROS_B.rtb_Switch_idx_3 /= multiModeQuad_ROS_B.K12;

      // Gain: '<S57>/Gain1' incorporates:
      //   Merge: '<S20>/Merge'
      //   Product: '<S57>/Product'
      //   Sum: '<S63>/Add'

      multiModeQuad_ROS_B.Merge[2] = (multiModeQuad_ROS_B.VectorConcatenate_m[1]
        + multiModeQuad_ROS_B.VectorConcatenate_m[3]) *
        multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Gain1_Gain_b;

      // Gain: '<S57>/Gain2' incorporates:
      //   Fcn: '<S43>/Fcn13'
      //   Merge: '<S20>/Merge'
      //   Product: '<S57>/Product'
      //   Sum: '<S61>/Add'

      multiModeQuad_ROS_B.Merge[3] = (multiModeQuad_ROS_B.VectorConcatenate_m[2]
        + -multiModeQuad_ROS_B.Selector_tmp) *
        multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Gain2_Gain_d;

      // Gain: '<S57>/Gain3' incorporates:
      //   Merge: '<S20>/Merge'
      //   Product: '<S57>/Product'
      //   Sum: '<S62>/Add'

      multiModeQuad_ROS_B.Merge[0] = (multiModeQuad_ROS_B.VectorConcatenate_m[7]
        - multiModeQuad_ROS_B.VectorConcatenate_m[5]) *
        multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Gain3_Gain_az;

      // Gain: '<S57>/Gain' incorporates:
      //   Merge: '<S20>/Merge'

      multiModeQuad_ROS_B.Merge[1] = multiModeQuad_ROS_P.Gain_Gain_i *
        multiModeQuad_ROS_B.K12_f;

      // End of Outputs for SubSystem: '<S53>/Maximum Value at DCM(1,1)'
      break;
    }

    // End of If: '<S53>/Find Maximum Diagonal Value'
    // End of Outputs for SubSystem: '<S20>/Negative Trace'
    break;
  }

  // End of If: '<S20>/If'

  // MATLAB Function: '<Root>/time to sec & nsec' incorporates:
  //   Clock: '<Root>/Clock'

  multiModeQuad_ROS_timetosecnsec(multiModeQuad_ROS_M->Timing.t[0],
    &multiModeQuad_ROS_B.rtb_Switch_idx_3, &multiModeQuad_ROS_B.K12_f);

  // BusAssignment: '<Root>/Bus Assignment' incorporates:
  //   Constant: '<S1>/Constant'
  //   Gain: '<S34>/Gain'
  //   Integrator: '<S18>/xe,ye,ze'
  //   MATLABSystem: '<Root>/Get Parameter'
  //   MATLABSystem: '<Root>/Get Parameter1'
  //   MATLABSystem: '<Root>/Get Parameter2'
  //   Sum: '<Root>/Sum'

  multiModeQuad_ROS_B.BusAssignment = multiModeQuad_ROS_P.Constant_Value_f;
  multiModeQuad_ROS_B.BusAssignment.Header.Stamp.Sec =
    multiModeQuad_ROS_B.rtb_Switch_idx_3;
  multiModeQuad_ROS_B.BusAssignment.Header.Stamp.Nsec =
    multiModeQuad_ROS_B.K12_f;
  multiModeQuad_ROS_B.BusAssignment.Pose.Position.X =
    multiModeQuad_ROS_B.rtb_Switch_idx_0 + multiModeQuad_ROS_X.xeyeze_CSTATE[0];
  multiModeQuad_ROS_B.BusAssignment.Pose.Position.Y =
    multiModeQuad_ROS_B.rtb_Switch_idx_1 + multiModeQuad_ROS_X.xeyeze_CSTATE[1];
  multiModeQuad_ROS_B.BusAssignment.Pose.Position.Z =
    multiModeQuad_ROS_P.Gain_Gain_m * multiModeQuad_ROS_X.xeyeze_CSTATE[2] +
    multiModeQuad_ROS_B.rtb_sincos_o2_idx_0;
  multiModeQuad_ROS_B.BusAssignment.Pose.Orientation.W =
    multiModeQuad_ROS_B.Merge[0];
  multiModeQuad_ROS_B.BusAssignment.Pose.Orientation.X =
    multiModeQuad_ROS_B.Merge[1];
  multiModeQuad_ROS_B.BusAssignment.Pose.Orientation.Y =
    multiModeQuad_ROS_B.Merge[2];
  multiModeQuad_ROS_B.BusAssignment.Pose.Orientation.Z =
    multiModeQuad_ROS_B.Merge[3];

  // Outputs for Atomic SubSystem: '<Root>/Publish'
  // MATLABSystem: '<S4>/SinkBlock'
  Pub_multiModeQuad_ROS_438.publish(&multiModeQuad_ROS_B.BusAssignment);

  // End of Outputs for SubSystem: '<Root>/Publish'

  // MATLAB Function: '<Root>/time to sec & nsec1' incorporates:
  //   Clock: '<Root>/Clock1'

  multiModeQuad_ROS_timetosecnsec(multiModeQuad_ROS_M->Timing.t[0],
    &multiModeQuad_ROS_B.rtb_Switch_idx_3, &multiModeQuad_ROS_B.K12_f);

  // BusAssignment: '<Root>/Bus Assignment1' incorporates:
  //   Constant: '<S2>/Constant'
  //   Integrator: '<S18>/p,q,r '

  multiModeQuad_ROS_B.BusAssignment1 = multiModeQuad_ROS_P.Constant_Value;
  multiModeQuad_ROS_B.BusAssignment1.Header.Stamp.Sec =
    multiModeQuad_ROS_B.rtb_Switch_idx_3;
  multiModeQuad_ROS_B.BusAssignment1.Header.Stamp.Nsec =
    multiModeQuad_ROS_B.K12_f;
  multiModeQuad_ROS_B.BusAssignment1.Orientation.W = multiModeQuad_ROS_B.Merge[0];
  multiModeQuad_ROS_B.BusAssignment1.Orientation.X = multiModeQuad_ROS_B.Merge[1];
  multiModeQuad_ROS_B.BusAssignment1.Orientation.Y = multiModeQuad_ROS_B.Merge[2];
  multiModeQuad_ROS_B.BusAssignment1.Orientation.Z = multiModeQuad_ROS_B.Merge[3];
  multiModeQuad_ROS_B.BusAssignment1.AngularVelocity.X =
    multiModeQuad_ROS_X.pqr_CSTATE[0];
  multiModeQuad_ROS_B.BusAssignment1.AngularVelocity.Y =
    multiModeQuad_ROS_X.pqr_CSTATE[1];
  multiModeQuad_ROS_B.BusAssignment1.AngularVelocity.Z =
    multiModeQuad_ROS_X.pqr_CSTATE[2];
  multiModeQuad_ROS_B.BusAssignment1.LinearAcceleration.X =
    multiModeQuad_ROS_B.Sum[0];
  multiModeQuad_ROS_B.BusAssignment1.LinearAcceleration.Y =
    multiModeQuad_ROS_B.Sum[1];
  multiModeQuad_ROS_B.BusAssignment1.LinearAcceleration.Z =
    multiModeQuad_ROS_B.Sum[2];

  // Outputs for Atomic SubSystem: '<Root>/Publish1'
  // MATLABSystem: '<S5>/SinkBlock'
  Pub_multiModeQuad_ROS_477.publish(&multiModeQuad_ROS_B.BusAssignment1);

  // End of Outputs for SubSystem: '<Root>/Publish1'
  if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
    // Gain: '<S365>/Filter Coefficient' incorporates:
    //   DiscreteIntegrator: '<S357>/Filter'
    //   Gain: '<S356>/Derivative Gain'
    //   Sum: '<S357>/SumD'

    multiModeQuad_ROS_B.FilterCoefficient_g =
      (multiModeQuad_ROS_P.PIDangularroll_D * multiModeQuad_ROS_B.sincos_o1[0] -
       multiModeQuad_ROS_DW.Filter_DSTATE) *
      multiModeQuad_ROS_P.PIDangularroll_N;

    // Gain: '<S360>/Proportional Gain' incorporates:
    //   DiscreteIntegrator: '<S362>/Integrator'
    //   Sum: '<S371>/Sum'

    multiModeQuad_ROS_B.ProportionalGain = ((multiModeQuad_ROS_B.sincos_o1[0] +
      multiModeQuad_ROS_DW.Integrator_DSTATE) +
      multiModeQuad_ROS_B.FilterCoefficient_g) *
      multiModeQuad_ROS_P.PIDangularroll_P;

    // Gain: '<S317>/Filter Coefficient' incorporates:
    //   DiscreteIntegrator: '<S309>/Filter'
    //   Gain: '<S308>/Derivative Gain'
    //   Sum: '<S309>/SumD'

    multiModeQuad_ROS_B.FilterCoefficient_a =
      (multiModeQuad_ROS_P.PIDangularpitch_D * multiModeQuad_ROS_B.sincos_o1[1]
       - multiModeQuad_ROS_DW.Filter_DSTATE_a) *
      multiModeQuad_ROS_P.PIDangularpitch_N;

    // Gain: '<S312>/Proportional Gain' incorporates:
    //   DiscreteIntegrator: '<S314>/Integrator'
    //   Sum: '<S323>/Sum'

    multiModeQuad_ROS_B.ProportionalGain_o = ((multiModeQuad_ROS_B.sincos_o1[1]
      + multiModeQuad_ROS_DW.Integrator_DSTATE_f) +
      multiModeQuad_ROS_B.FilterCoefficient_a) *
      multiModeQuad_ROS_P.PIDangularpitch_P;

    // Gain: '<S269>/Filter Coefficient' incorporates:
    //   DiscreteIntegrator: '<S261>/Filter'
    //   Gain: '<S260>/Derivative Gain'
    //   Sum: '<S261>/SumD'

    multiModeQuad_ROS_B.FilterCoefficient_i =
      (multiModeQuad_ROS_P.PIDangulayaw_D * multiModeQuad_ROS_B.sincos_o1[2] -
       multiModeQuad_ROS_DW.Filter_DSTATE_c) *
      multiModeQuad_ROS_P.PIDangulayaw_N;

    // Gain: '<S264>/Proportional Gain' incorporates:
    //   DiscreteIntegrator: '<S266>/Integrator'
    //   Sum: '<S275>/Sum'

    multiModeQuad_ROS_B.ProportionalGain_l = ((multiModeQuad_ROS_B.sincos_o1[2]
      + multiModeQuad_ROS_DW.Integrator_DSTATE_l) +
      multiModeQuad_ROS_B.FilterCoefficient_i) *
      multiModeQuad_ROS_P.PIDangulayaw_P;

    // Outputs for Atomic SubSystem: '<Root>/Sub setpoint thrust'
    // MATLABSystem: '<S8>/SourceBlock' incorporates:
    //   Inport: '<S16>/In1'

    tmp = Sub_multiModeQuad_ROS_500.getLatestMessage
      (&multiModeQuad_ROS_B.b_varargout_2_j);

    // Outputs for Enabled SubSystem: '<S8>/Enabled Subsystem' incorporates:
    //   EnablePort: '<S16>/Enable'

    if (tmp) {
      multiModeQuad_ROS_B.In1_ol = multiModeQuad_ROS_B.b_varargout_2_j;
    }

    // End of MATLABSystem: '<S8>/SourceBlock'
    // End of Outputs for SubSystem: '<S8>/Enabled Subsystem'
    // End of Outputs for SubSystem: '<Root>/Sub setpoint thrust'

    // Gain: '<S10>/Gain'
    multiModeQuad_ROS_B.Gain = multiModeQuad_ROS_P.Gain_Gain_dx *
      multiModeQuad_ROS_B.In1_ol.Data;
  }

  // Gain: '<S221>/Filter Coefficient' incorporates:
  //   Gain: '<S212>/Derivative Gain'
  //   Integrator: '<S213>/Filter'
  //   Sum: '<S213>/SumD'

  multiModeQuad_ROS_B.FilterCoefficient_j = (multiModeQuad_ROS_P.PIDVelocityz_D *
    multiModeQuad_ROS_B.rtb_Sum2_idx_2 - multiModeQuad_ROS_X.Filter_CSTATE_h) *
    multiModeQuad_ROS_P.PIDVelocityz_N;

  // Switch: '<S10>/Switch2' incorporates:
  //   Constant: '<S10>/Hover throttle'
  //   Gain: '<S223>/Proportional Gain'
  //   Integrator: '<S218>/Integrator'
  //   Sum: '<S10>/Sum3'
  //   Sum: '<S227>/Sum'

  if (multiModeQuad_ROS_B.In1_e.Data >= multiModeQuad_ROS_P.Switch2_Threshold) {
    multiModeQuad_ROS_B.rtb_Switch_idx_1 = multiModeQuad_ROS_B.Gain;
  } else {
    multiModeQuad_ROS_B.rtb_Switch_idx_1 = ((multiModeQuad_ROS_P.PIDVelocityz_P *
      multiModeQuad_ROS_B.rtb_Sum2_idx_2 +
      multiModeQuad_ROS_X.Integrator_CSTATE_h) +
      multiModeQuad_ROS_B.FilterCoefficient_j) +
      multiModeQuad_ROS_P.Hoverthrottle_Value;
  }

  // End of Switch: '<S10>/Switch2'

  // MATLAB Function: '<S10>/Mixer'
  multiModeQuad_ROS_B.rtb_Switch_idx_0 = multiModeQuad_ROS_B.rtb_Switch_idx_1 -
    multiModeQuad_ROS_B.ProportionalGain;
  multiModeQuad_ROS_B.rtb_Switch_idx_3 = (multiModeQuad_ROS_B.rtb_Switch_idx_0 -
    multiModeQuad_ROS_B.ProportionalGain_o) -
    multiModeQuad_ROS_B.ProportionalGain_l;

  // Saturate: '<S10>/Saturation4'
  if (multiModeQuad_ROS_B.rtb_Switch_idx_3 >
      multiModeQuad_ROS_P.Saturation4_UpperSat) {
    // Saturate: '<S10>/Saturation7'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_P.Saturation4_UpperSat;
  } else if (multiModeQuad_ROS_B.rtb_Switch_idx_3 <
             multiModeQuad_ROS_P.Saturation4_LowerSat) {
    // Saturate: '<S10>/Saturation7'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_P.Saturation4_LowerSat;
  }

  // End of Saturate: '<S10>/Saturation4'

  // Gain: '<S10>/Komega' incorporates:
  //   Math: '<S10>/Square'

  multiModeQuad_ROS_B.Komega = multiModeQuad_ROS_B.rtb_Switch_idx_3 *
    multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Komega_Gain;

  // MATLAB Function: '<S10>/Mixer'
  multiModeQuad_ROS_B.rtb_Switch_idx_3 = (multiModeQuad_ROS_B.rtb_Switch_idx_0 +
    multiModeQuad_ROS_B.ProportionalGain_o) +
    multiModeQuad_ROS_B.ProportionalGain_l;

  // Saturate: '<S10>/Saturation5'
  if (multiModeQuad_ROS_B.rtb_Switch_idx_3 >
      multiModeQuad_ROS_P.Saturation5_UpperSat) {
    // Saturate: '<S10>/Saturation7'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_P.Saturation5_UpperSat;
  } else if (multiModeQuad_ROS_B.rtb_Switch_idx_3 <
             multiModeQuad_ROS_P.Saturation5_LowerSat) {
    // Saturate: '<S10>/Saturation7'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_P.Saturation5_LowerSat;
  }

  // End of Saturate: '<S10>/Saturation5'

  // Gain: '<S10>/Komega1' incorporates:
  //   Math: '<S10>/Square1'

  multiModeQuad_ROS_B.Komega1 = multiModeQuad_ROS_B.rtb_Switch_idx_3 *
    multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Komega1_Gain;

  // MATLAB Function: '<S10>/Mixer'
  multiModeQuad_ROS_B.rtb_Switch_idx_0 = multiModeQuad_ROS_B.rtb_Switch_idx_1 +
    multiModeQuad_ROS_B.ProportionalGain;
  multiModeQuad_ROS_B.rtb_Switch_idx_3 = (multiModeQuad_ROS_B.rtb_Switch_idx_0 -
    multiModeQuad_ROS_B.ProportionalGain_o) +
    multiModeQuad_ROS_B.ProportionalGain_l;

  // Saturate: '<S10>/Saturation6'
  if (multiModeQuad_ROS_B.rtb_Switch_idx_3 >
      multiModeQuad_ROS_P.Saturation6_UpperSat) {
    // Saturate: '<S10>/Saturation7'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_P.Saturation6_UpperSat;
  } else if (multiModeQuad_ROS_B.rtb_Switch_idx_3 <
             multiModeQuad_ROS_P.Saturation6_LowerSat) {
    // Saturate: '<S10>/Saturation7'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_P.Saturation6_LowerSat;
  }

  // End of Saturate: '<S10>/Saturation6'

  // Gain: '<S10>/Komega2' incorporates:
  //   Math: '<S10>/Square2'

  multiModeQuad_ROS_B.Komega2 = multiModeQuad_ROS_B.rtb_Switch_idx_3 *
    multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Komega2_Gain;

  // MATLAB Function: '<S10>/Mixer'
  multiModeQuad_ROS_B.rtb_Switch_idx_3 = (multiModeQuad_ROS_B.rtb_Switch_idx_0 +
    multiModeQuad_ROS_B.ProportionalGain_o) -
    multiModeQuad_ROS_B.ProportionalGain_l;

  // Saturate: '<S10>/Saturation7'
  if (multiModeQuad_ROS_B.rtb_Switch_idx_3 >
      multiModeQuad_ROS_P.Saturation7_UpperSat) {
    // Saturate: '<S10>/Saturation7'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_P.Saturation7_UpperSat;
  } else if (multiModeQuad_ROS_B.rtb_Switch_idx_3 <
             multiModeQuad_ROS_P.Saturation7_LowerSat) {
    // Saturate: '<S10>/Saturation7'
    multiModeQuad_ROS_B.rtb_Switch_idx_3 =
      multiModeQuad_ROS_P.Saturation7_LowerSat;
  }

  // End of Saturate: '<S10>/Saturation7'

  // Gain: '<S10>/Komega3' incorporates:
  //   Math: '<S10>/Square3'

  multiModeQuad_ROS_B.Komega3 = multiModeQuad_ROS_B.rtb_Switch_idx_3 *
    multiModeQuad_ROS_B.rtb_Switch_idx_3 * multiModeQuad_ROS_P.Komega3_Gain;

  // Gain: '<S119>/Integral Gain'
  multiModeQuad_ROS_B.IntegralGain = multiModeQuad_ROS_P.PIDVelocityx_I *
    multiModeQuad_ROS_B.rtb_Sum2_idx_0;

  // Gain: '<S167>/Integral Gain'
  multiModeQuad_ROS_B.IntegralGain_i = multiModeQuad_ROS_P.PIDVelocityy_I *
    multiModeQuad_ROS_B.rtb_Sum2_idx_1;

  // Gain: '<S215>/Integral Gain'
  multiModeQuad_ROS_B.IntegralGain_g = multiModeQuad_ROS_P.PIDVelocityz_I *
    multiModeQuad_ROS_B.rtb_Sum2_idx_2;
  if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
    int8_T rtAction;

    // Gain: '<S263>/Integral Gain'
    multiModeQuad_ROS_B.IntegralGain_g1 = multiModeQuad_ROS_P.PIDangulayaw_I *
      multiModeQuad_ROS_B.sincos_o1[2];

    // Gain: '<S311>/Integral Gain'
    multiModeQuad_ROS_B.IntegralGain_i_m = multiModeQuad_ROS_P.PIDangularpitch_I
      * multiModeQuad_ROS_B.sincos_o1[1];

    // Gain: '<S359>/Integral Gain'
    multiModeQuad_ROS_B.IntegralGain_j = multiModeQuad_ROS_P.PIDangularroll_I *
      multiModeQuad_ROS_B.sincos_o1[0];

    // If: '<S55>/If1' incorporates:
    //   Constant: '<S55>/Constant'

    rtAction = -1;
    if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
      if (multiModeQuad_ROS_P.DirectionCosineMatrixtoQuaterni != 1.0) {
        rtAction = 0;
      }

      multiModeQuad_ROS_DW.If1_ActiveSubsystem = rtAction;
    } else {
      rtAction = multiModeQuad_ROS_DW.If1_ActiveSubsystem;
    }

    if (rtAction == 0) {
      // Outputs for IfAction SubSystem: '<S55>/If Warning//Error' incorporates:
      //   ActionPort: '<S79>/if'

      // Bias: '<S82>/Bias1' incorporates:
      //   Concatenate: '<S45>/Vector Concatenate'
      //   Math: '<S82>/Math Function'
      //   Product: '<S82>/Product'

      for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 3;
           multiModeQuad_ROS_B.c_i++) {
        for (multiModeQuad_ROS_B.c_j = 0; multiModeQuad_ROS_B.c_j < 3;
             multiModeQuad_ROS_B.c_j++) {
          multiModeQuad_ROS_B.sgn = 3 * multiModeQuad_ROS_B.c_j +
            multiModeQuad_ROS_B.c_i;
          multiModeQuad_ROS_B.Product_tmp_c[multiModeQuad_ROS_B.sgn] =
            ((multiModeQuad_ROS_B.VectorConcatenate_m[3 *
              multiModeQuad_ROS_B.c_j + 1] *
              multiModeQuad_ROS_B.Product_tmp[multiModeQuad_ROS_B.c_i + 3] +
              multiModeQuad_ROS_B.VectorConcatenate_m[3 *
              multiModeQuad_ROS_B.c_j] *
              multiModeQuad_ROS_B.Product_tmp[multiModeQuad_ROS_B.c_i]) +
             multiModeQuad_ROS_B.VectorConcatenate_m[3 * multiModeQuad_ROS_B.c_j
             + 2] * multiModeQuad_ROS_B.Product_tmp[multiModeQuad_ROS_B.c_i + 6])
            + multiModeQuad_ROS_P.Bias1_Bias[multiModeQuad_ROS_B.sgn];
        }
      }

      // End of Bias: '<S82>/Bias1'

      // RelationalOperator: '<S88>/Compare' incorporates:
      //   Abs: '<S82>/Abs2'
      //   Constant: '<S88>/Constant'

      for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 9;
           multiModeQuad_ROS_B.c_i++) {
        multiModeQuad_ROS_B.Compare[multiModeQuad_ROS_B.c_i] = (fabs
          (multiModeQuad_ROS_B.Product_tmp_c[multiModeQuad_ROS_B.c_i]) >
          multiModeQuad_ROS_P.DirectionCosineMatrixtoQuater_h);
      }

      // End of RelationalOperator: '<S88>/Compare'

      // Logic: '<S82>/Logical Operator1' incorporates:
      //   RelationalOperator: '<S88>/Compare'

      tmp = multiModeQuad_ROS_B.Compare[0];
      for (multiModeQuad_ROS_B.c_i = 0; multiModeQuad_ROS_B.c_i < 8;
           multiModeQuad_ROS_B.c_i++) {
        tmp = (tmp || multiModeQuad_ROS_B.Compare[multiModeQuad_ROS_B.c_i + 1]);
      }

      // If: '<S79>/If' incorporates:
      //   Abs: '<S83>/Abs1'
      //   Bias: '<S83>/Bias'
      //   Concatenate: '<S45>/Vector Concatenate'
      //   Constant: '<S90>/Constant'
      //   Fcn: '<S43>/Fcn13'
      //   Logic: '<S82>/Logical Operator1'
      //   Product: '<S89>/Product'
      //   Product: '<S89>/Product1'
      //   Product: '<S89>/Product2'
      //   Product: '<S89>/Product3'
      //   Product: '<S89>/Product4'
      //   Product: '<S89>/Product5'
      //   RelationalOperator: '<S90>/Compare'
      //   Reshape: '<S89>/Reshape'
      //   Sum: '<S89>/Sum'

      if (fabs((((((multiModeQuad_ROS_B.VectorConcatenate_m[0] *
                    multiModeQuad_ROS_B.VectorConcatenate_m[4] *
                    multiModeQuad_ROS_B.VectorConcatenate_m[8] -
                    multiModeQuad_ROS_B.VectorConcatenate_m[0] *
                    multiModeQuad_ROS_B.VectorConcatenate_m[5] *
                    multiModeQuad_ROS_B.VectorConcatenate_m[7]) -
                   multiModeQuad_ROS_B.VectorConcatenate_m[1] *
                   multiModeQuad_ROS_B.VectorConcatenate_m[3] *
                   multiModeQuad_ROS_B.VectorConcatenate_m[8]) +
                  multiModeQuad_ROS_B.VectorConcatenate_m[2] *
                  multiModeQuad_ROS_B.VectorConcatenate_m[3] *
                  multiModeQuad_ROS_B.VectorConcatenate_m[7]) +
                 multiModeQuad_ROS_B.VectorConcatenate_m[1] *
                 multiModeQuad_ROS_B.VectorConcatenate_m[5] *
                 -multiModeQuad_ROS_B.Selector_tmp) -
                multiModeQuad_ROS_B.VectorConcatenate_m[2] *
                multiModeQuad_ROS_B.VectorConcatenate_m[4] *
                -multiModeQuad_ROS_B.Selector_tmp) +
               multiModeQuad_ROS_P.Bias_Bias) >
          multiModeQuad_ROS_P.DirectionCosineMatrixtoQuater_h) {
        // Outputs for IfAction SubSystem: '<S79>/If Not Proper' incorporates:
        //   ActionPort: '<S81>/Action Port'

        // If: '<S81>/If' incorporates:
        //   Constant: '<S81>/Constant'

        if (multiModeQuad_ROS_P.DirectionCosineMatrixtoQuaterni == 2.0) {
          // Outputs for IfAction SubSystem: '<S81>/Warning' incorporates:
          //   ActionPort: '<S87>/Action Port'

          // Assertion: '<S87>/Assertion' incorporates:
          //   Constant: '<S81>/Constant1'

          utAssert(multiModeQuad_ROS_P.Constant1_Value_d != 0.0);

          // End of Outputs for SubSystem: '<S81>/Warning'
        } else if (multiModeQuad_ROS_P.DirectionCosineMatrixtoQuaterni == 3.0) {
          // Outputs for IfAction SubSystem: '<S81>/Error' incorporates:
          //   ActionPort: '<S86>/Action Port'

          // Assertion: '<S86>/Assertion' incorporates:
          //   Constant: '<S81>/Constant1'

          utAssert(multiModeQuad_ROS_P.Constant1_Value_d != 0.0);

          // End of Outputs for SubSystem: '<S81>/Error'
        }

        // End of If: '<S81>/If'
        // End of Outputs for SubSystem: '<S79>/If Not Proper'
      } else if (tmp) {
        // Outputs for IfAction SubSystem: '<S79>/Else If Not Orthogonal' incorporates:
        //   ActionPort: '<S80>/Action Port'

        // If: '<S80>/If' incorporates:
        //   Constant: '<S80>/Constant'

        if (multiModeQuad_ROS_P.DirectionCosineMatrixtoQuaterni == 2.0) {
          // Outputs for IfAction SubSystem: '<S80>/Warning' incorporates:
          //   ActionPort: '<S85>/Action Port'

          // Assertion: '<S85>/Assertion' incorporates:
          //   Constant: '<S80>/Constant1'

          utAssert(multiModeQuad_ROS_P.Constant1_Value_j != 0.0);

          // End of Outputs for SubSystem: '<S80>/Warning'
        } else if (multiModeQuad_ROS_P.DirectionCosineMatrixtoQuaterni == 3.0) {
          // Outputs for IfAction SubSystem: '<S80>/Error' incorporates:
          //   ActionPort: '<S84>/Action Port'

          // Assertion: '<S84>/Assertion' incorporates:
          //   Constant: '<S80>/Constant1'

          utAssert(multiModeQuad_ROS_P.Constant1_Value_j != 0.0);

          // End of Outputs for SubSystem: '<S80>/Error'
        }

        // End of If: '<S80>/If'
        // End of Outputs for SubSystem: '<S79>/Else If Not Orthogonal'
      }

      // End of If: '<S79>/If'
      // End of Outputs for SubSystem: '<S55>/If Warning//Error'
    }
  }

  if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
    if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
      // Update for DiscreteIntegrator: '<S357>/Filter'
      multiModeQuad_ROS_DW.Filter_DSTATE += multiModeQuad_ROS_P.Filter_gainval *
        multiModeQuad_ROS_B.FilterCoefficient_g;

      // Update for DiscreteIntegrator: '<S362>/Integrator'
      multiModeQuad_ROS_DW.Integrator_DSTATE +=
        multiModeQuad_ROS_P.Integrator_gainval *
        multiModeQuad_ROS_B.IntegralGain_j;

      // Update for DiscreteIntegrator: '<S314>/Integrator'
      multiModeQuad_ROS_DW.Integrator_DSTATE_f +=
        multiModeQuad_ROS_P.Integrator_gainval_e *
        multiModeQuad_ROS_B.IntegralGain_i_m;

      // Update for DiscreteIntegrator: '<S309>/Filter'
      multiModeQuad_ROS_DW.Filter_DSTATE_a +=
        multiModeQuad_ROS_P.Filter_gainval_j *
        multiModeQuad_ROS_B.FilterCoefficient_a;

      // Update for DiscreteIntegrator: '<S266>/Integrator'
      multiModeQuad_ROS_DW.Integrator_DSTATE_l +=
        multiModeQuad_ROS_P.Integrator_gainval_a *
        multiModeQuad_ROS_B.IntegralGain_g1;

      // Update for DiscreteIntegrator: '<S261>/Filter'
      multiModeQuad_ROS_DW.Filter_DSTATE_c +=
        multiModeQuad_ROS_P.Filter_gainval_c *
        multiModeQuad_ROS_B.FilterCoefficient_i;
    }
  }                                    // end MajorTimeStep

  if (rtmIsMajorTimeStep(multiModeQuad_ROS_M)) {
    rt_ertODEUpdateContinuousStates(&multiModeQuad_ROS_M->solverInfo);

    // Update absolute time for base rate
    // The "clockTick0" counts the number of times the code of this task has
    //  been executed. The absolute time is the multiplication of "clockTick0"
    //  and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
    //  overflow during the application lifespan selected.

    ++multiModeQuad_ROS_M->Timing.clockTick0;
    multiModeQuad_ROS_M->Timing.t[0] = rtsiGetSolverStopTime
      (&multiModeQuad_ROS_M->solverInfo);

    {
      // Update absolute timer for sample time: [0.005s, 0.0s]
      // The "clockTick1" counts the number of times the code of this task has
      //  been executed. The resolution of this integer timer is 0.005, which is the step size
      //  of the task. Size of "clockTick1" ensures timer will not overflow during the
      //  application lifespan selected.

      multiModeQuad_ROS_M->Timing.clockTick1++;
    }
  }                                    // end MajorTimeStep
}

// Derivatives for root system: '<Root>'
void multiModeQuad_ROS_derivatives(void)
{
  XDot_multiModeQuad_ROS_T *_rtXdot;
  _rtXdot = ((XDot_multiModeQuad_ROS_T *) multiModeQuad_ROS_M->derivs);

  // Derivatives for Integrator: '<S122>/Integrator'
  _rtXdot->Integrator_CSTATE = multiModeQuad_ROS_B.IntegralGain;

  // Derivatives for Integrator: '<S117>/Filter'
  _rtXdot->Filter_CSTATE = multiModeQuad_ROS_B.FilterCoefficient;

  // Derivatives for Integrator: '<S170>/Integrator'
  _rtXdot->Integrator_CSTATE_o = multiModeQuad_ROS_B.IntegralGain_i;

  // Derivatives for Integrator: '<S165>/Filter'
  _rtXdot->Filter_CSTATE_k = multiModeQuad_ROS_B.FilterCoefficient_m;

  // Derivatives for TransferFcn: '<S10>/Transfer Fcn2'
  _rtXdot->TransferFcn2_CSTATE = 0.0;
  _rtXdot->TransferFcn2_CSTATE += multiModeQuad_ROS_P.TransferFcn2_A *
    multiModeQuad_ROS_X.TransferFcn2_CSTATE;
  _rtXdot->TransferFcn2_CSTATE += multiModeQuad_ROS_B.Komega;

  // Derivatives for TransferFcn: '<S10>/Transfer Fcn1'
  _rtXdot->TransferFcn1_CSTATE = 0.0;
  _rtXdot->TransferFcn1_CSTATE += multiModeQuad_ROS_P.TransferFcn1_A *
    multiModeQuad_ROS_X.TransferFcn1_CSTATE;
  _rtXdot->TransferFcn1_CSTATE += multiModeQuad_ROS_B.Komega1;

  // Derivatives for TransferFcn: '<S10>/Transfer Fcn3'
  _rtXdot->TransferFcn3_CSTATE = 0.0;
  _rtXdot->TransferFcn3_CSTATE += multiModeQuad_ROS_P.TransferFcn3_A *
    multiModeQuad_ROS_X.TransferFcn3_CSTATE;
  _rtXdot->TransferFcn3_CSTATE += multiModeQuad_ROS_B.Komega2;

  // Derivatives for TransferFcn: '<S10>/Transfer Fcn4'
  _rtXdot->TransferFcn4_CSTATE = 0.0;
  _rtXdot->TransferFcn4_CSTATE += multiModeQuad_ROS_P.TransferFcn4_A *
    multiModeQuad_ROS_X.TransferFcn4_CSTATE;
  _rtXdot->TransferFcn4_CSTATE += multiModeQuad_ROS_B.Komega3;

  // Derivatives for Integrator: '<S35>/phi theta psi'
  _rtXdot->phithetapsi_CSTATE[0] =
    multiModeQuad_ROS_B.TmpSignalConversionAtphithetaps[0];

  // Derivatives for Integrator: '<S18>/ub,vb,wb'
  _rtXdot->ubvbwb_CSTATE[0] = multiModeQuad_ROS_B.Sum[0];

  // Derivatives for Integrator: '<S18>/p,q,r '
  _rtXdot->pqr_CSTATE[0] = multiModeQuad_ROS_B.Product2[0];

  // Derivatives for Integrator: '<S18>/xe,ye,ze'
  _rtXdot->xeyeze_CSTATE[0] = multiModeQuad_ROS_B.Product[0];

  // Derivatives for Integrator: '<S35>/phi theta psi'
  _rtXdot->phithetapsi_CSTATE[1] =
    multiModeQuad_ROS_B.TmpSignalConversionAtphithetaps[1];

  // Derivatives for Integrator: '<S18>/ub,vb,wb'
  _rtXdot->ubvbwb_CSTATE[1] = multiModeQuad_ROS_B.Sum[1];

  // Derivatives for Integrator: '<S18>/p,q,r '
  _rtXdot->pqr_CSTATE[1] = multiModeQuad_ROS_B.Product2[1];

  // Derivatives for Integrator: '<S18>/xe,ye,ze'
  _rtXdot->xeyeze_CSTATE[1] = multiModeQuad_ROS_B.Product[1];

  // Derivatives for Integrator: '<S35>/phi theta psi'
  _rtXdot->phithetapsi_CSTATE[2] =
    multiModeQuad_ROS_B.TmpSignalConversionAtphithetaps[2];

  // Derivatives for Integrator: '<S18>/ub,vb,wb'
  _rtXdot->ubvbwb_CSTATE[2] = multiModeQuad_ROS_B.Sum[2];

  // Derivatives for Integrator: '<S18>/p,q,r '
  _rtXdot->pqr_CSTATE[2] = multiModeQuad_ROS_B.Product2[2];

  // Derivatives for Integrator: '<S18>/xe,ye,ze'
  _rtXdot->xeyeze_CSTATE[2] = multiModeQuad_ROS_B.Product[2];

  // Derivatives for Integrator: '<S218>/Integrator'
  _rtXdot->Integrator_CSTATE_h = multiModeQuad_ROS_B.IntegralGain_g;

  // Derivatives for Integrator: '<S213>/Filter'
  _rtXdot->Filter_CSTATE_h = multiModeQuad_ROS_B.FilterCoefficient_j;
}

// Model initialize function
void multiModeQuad_ROS_initialize(void)
{
  // Registration code

  // initialize non-finites
  rt_InitInfAndNaN(sizeof(real_T));

  // non-finite (run-time) assignments
  multiModeQuad_ROS_P.Saturation4_UpperSat = rtInf;
  multiModeQuad_ROS_P.Saturation5_UpperSat = rtInf;
  multiModeQuad_ROS_P.Saturation6_UpperSat = rtInf;
  multiModeQuad_ROS_P.Saturation7_UpperSat = rtInf;

  {
    // Setup solver object
    rtsiSetSimTimeStepPtr(&multiModeQuad_ROS_M->solverInfo,
                          &multiModeQuad_ROS_M->Timing.simTimeStep);
    rtsiSetTPtr(&multiModeQuad_ROS_M->solverInfo, &rtmGetTPtr
                (multiModeQuad_ROS_M));
    rtsiSetStepSizePtr(&multiModeQuad_ROS_M->solverInfo,
                       &multiModeQuad_ROS_M->Timing.stepSize0);
    rtsiSetdXPtr(&multiModeQuad_ROS_M->solverInfo, &multiModeQuad_ROS_M->derivs);
    rtsiSetContStatesPtr(&multiModeQuad_ROS_M->solverInfo, (real_T **)
                         &multiModeQuad_ROS_M->contStates);
    rtsiSetNumContStatesPtr(&multiModeQuad_ROS_M->solverInfo,
      &multiModeQuad_ROS_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&multiModeQuad_ROS_M->solverInfo,
      &multiModeQuad_ROS_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&multiModeQuad_ROS_M->solverInfo,
      &multiModeQuad_ROS_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&multiModeQuad_ROS_M->solverInfo,
      &multiModeQuad_ROS_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&multiModeQuad_ROS_M->solverInfo, (&rtmGetErrorStatus
      (multiModeQuad_ROS_M)));
    rtsiSetRTModelPtr(&multiModeQuad_ROS_M->solverInfo, multiModeQuad_ROS_M);
  }

  rtsiSetSimTimeStep(&multiModeQuad_ROS_M->solverInfo, MAJOR_TIME_STEP);
  multiModeQuad_ROS_M->intgData.y = multiModeQuad_ROS_M->odeY;
  multiModeQuad_ROS_M->intgData.f[0] = multiModeQuad_ROS_M->odeF[0];
  multiModeQuad_ROS_M->intgData.f[1] = multiModeQuad_ROS_M->odeF[1];
  multiModeQuad_ROS_M->intgData.f[2] = multiModeQuad_ROS_M->odeF[2];
  multiModeQuad_ROS_M->contStates = ((X_multiModeQuad_ROS_T *)
    &multiModeQuad_ROS_X);
  multiModeQuad_ROS_M->periodicContStateIndices = ((int_T*)
    multiModeQuad_ROS_PeriodicIndX);
  multiModeQuad_ROS_M->periodicContStateRanges = ((real_T*)
    multiModeQuad_ROS_PeriodicRngX);
  rtsiSetSolverData(&multiModeQuad_ROS_M->solverInfo, static_cast<void *>
                    (&multiModeQuad_ROS_M->intgData));
  rtsiSetSolverName(&multiModeQuad_ROS_M->solverInfo,"ode3");
  rtmSetTPtr(multiModeQuad_ROS_M, &multiModeQuad_ROS_M->Timing.tArray[0]);
  multiModeQuad_ROS_M->Timing.stepSize0 = 0.005;

  {
    int32_T i;
    char_T b_zeroDelimTopic_0[26];
    char_T b_zeroDelimTopic[12];
    char_T b_zeroDelimTopic_2[8];
    char_T b_zeroDelimTopic_3[7];
    char_T b_zeroDelimTopic_1[5];
    static const char_T tmp[11] = { 'f', 'l', 'i', 'g', 'h', 't', '_', 'm', 'o',
      'd', 'e' };

    static const char_T tmp_0[26] = { 's', 'e', 't', 'p', 'o', 'i', 'n', 't',
      '_', 'a', 't', 't', 'i', 't', 'u', 'd', 'e', '/', 'a', 't', 't', 'i', 't',
      'u', 'd', 'e' };

    static const char_T tmp_1[35] = { 's', 'e', 't', 'p', 'o', 'i', 'n', 't',
      '_', 'v', 'e', 'l', 'o', 'c', 'i', 't', 'y', '/', 'c', 'm', 'd', '_', 'v',
      'e', 'l', '_', 'u', 'n', 's', 't', 'a', 'm', 'p', 'e', 'd' };

    static const char_T tmp_2[25] = { 's', 'e', 't', 'p', 'o', 'i', 'n', 't',
      '_', 'a', 't', 't', 'i', 't', 'u', 'd', 'e', '/', 'c', 'm', 'd', '_', 'v',
      'e', 'l' };

    static const char_T tmp_3[7] = { 'i', 'm', 'u', 'D', 'a', 't', 'a' };

    static const char_T tmp_4[6] = { 't', 'h', 'r', 'u', 's', 't' };

    static const char_T tmp_5[11] = { 'i', 'n', 'i', 't', '_', 'p', 'o', 's',
      'e', '_', 'x' };

    static const char_T tmp_6[11] = { 'i', 'n', 'i', 't', '_', 'p', 'o', 's',
      'e', '_', 'y' };

    static const char_T tmp_7[11] = { 'i', 'n', 'i', 't', '_', 'p', 'o', 's',
      'e', '_', 'z' };

    // Start for If: '<S20>/If'
    multiModeQuad_ROS_DW.If_ActiveSubsystem = -1;

    // Start for If: '<S55>/If1'
    multiModeQuad_ROS_DW.If1_ActiveSubsystem = -1;

    // InitializeConditions for Integrator: '<S122>/Integrator'
    multiModeQuad_ROS_X.Integrator_CSTATE =
      multiModeQuad_ROS_P.PIDVelocityx_InitialCondition_o;

    // InitializeConditions for Integrator: '<S117>/Filter'
    multiModeQuad_ROS_X.Filter_CSTATE =
      multiModeQuad_ROS_P.PIDVelocityx_InitialConditionFo;

    // InitializeConditions for Integrator: '<S170>/Integrator'
    multiModeQuad_ROS_X.Integrator_CSTATE_o =
      multiModeQuad_ROS_P.PIDVelocityy_InitialCondition_i;

    // InitializeConditions for Integrator: '<S165>/Filter'
    multiModeQuad_ROS_X.Filter_CSTATE_k =
      multiModeQuad_ROS_P.PIDVelocityy_InitialConditionFo;

    // InitializeConditions for TransferFcn: '<S10>/Transfer Fcn2'
    multiModeQuad_ROS_X.TransferFcn2_CSTATE = 0.0;

    // InitializeConditions for TransferFcn: '<S10>/Transfer Fcn1'
    multiModeQuad_ROS_X.TransferFcn1_CSTATE = 0.0;

    // InitializeConditions for TransferFcn: '<S10>/Transfer Fcn3'
    multiModeQuad_ROS_X.TransferFcn3_CSTATE = 0.0;

    // InitializeConditions for TransferFcn: '<S10>/Transfer Fcn4'
    multiModeQuad_ROS_X.TransferFcn4_CSTATE = 0.0;

    // InitializeConditions for Integrator: '<S35>/phi theta psi'
    multiModeQuad_ROS_X.phithetapsi_CSTATE[0] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_eul_0[0];

    // InitializeConditions for Integrator: '<S18>/ub,vb,wb'
    multiModeQuad_ROS_X.ubvbwb_CSTATE[0] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_Vm_0[0];

    // InitializeConditions for Integrator: '<S18>/p,q,r '
    multiModeQuad_ROS_X.pqr_CSTATE[0] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_pm_0[0];

    // InitializeConditions for Integrator: '<S18>/xe,ye,ze'
    multiModeQuad_ROS_X.xeyeze_CSTATE[0] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_xme_0[0];

    // InitializeConditions for Integrator: '<S35>/phi theta psi'
    multiModeQuad_ROS_X.phithetapsi_CSTATE[1] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_eul_0[1];

    // InitializeConditions for Integrator: '<S18>/ub,vb,wb'
    multiModeQuad_ROS_X.ubvbwb_CSTATE[1] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_Vm_0[1];

    // InitializeConditions for Integrator: '<S18>/p,q,r '
    multiModeQuad_ROS_X.pqr_CSTATE[1] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_pm_0[1];

    // InitializeConditions for Integrator: '<S18>/xe,ye,ze'
    multiModeQuad_ROS_X.xeyeze_CSTATE[1] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_xme_0[1];

    // InitializeConditions for Integrator: '<S35>/phi theta psi'
    multiModeQuad_ROS_X.phithetapsi_CSTATE[2] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_eul_0[2];

    // InitializeConditions for Integrator: '<S18>/ub,vb,wb'
    multiModeQuad_ROS_X.ubvbwb_CSTATE[2] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_Vm_0[2];

    // InitializeConditions for Integrator: '<S18>/p,q,r '
    multiModeQuad_ROS_X.pqr_CSTATE[2] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_pm_0[2];

    // InitializeConditions for Integrator: '<S18>/xe,ye,ze'
    multiModeQuad_ROS_X.xeyeze_CSTATE[2] =
      multiModeQuad_ROS_P.uDOFEulerAngles2_xme_0[2];

    // InitializeConditions for DiscreteIntegrator: '<S357>/Filter'
    multiModeQuad_ROS_DW.Filter_DSTATE =
      multiModeQuad_ROS_P.PIDangularroll_InitialCondition;

    // InitializeConditions for DiscreteIntegrator: '<S362>/Integrator'
    multiModeQuad_ROS_DW.Integrator_DSTATE =
      multiModeQuad_ROS_P.PIDangularroll_InitialConditi_o;

    // InitializeConditions for DiscreteIntegrator: '<S314>/Integrator'
    multiModeQuad_ROS_DW.Integrator_DSTATE_f =
      multiModeQuad_ROS_P.PIDangularpitch_InitialCondit_b;

    // InitializeConditions for DiscreteIntegrator: '<S309>/Filter'
    multiModeQuad_ROS_DW.Filter_DSTATE_a =
      multiModeQuad_ROS_P.PIDangularpitch_InitialConditio;

    // InitializeConditions for DiscreteIntegrator: '<S266>/Integrator'
    multiModeQuad_ROS_DW.Integrator_DSTATE_l =
      multiModeQuad_ROS_P.PIDangulayaw_InitialCondition_m;

    // InitializeConditions for DiscreteIntegrator: '<S261>/Filter'
    multiModeQuad_ROS_DW.Filter_DSTATE_c =
      multiModeQuad_ROS_P.PIDangulayaw_InitialConditionFo;

    // InitializeConditions for Integrator: '<S218>/Integrator'
    multiModeQuad_ROS_X.Integrator_CSTATE_h =
      multiModeQuad_ROS_P.PIDVelocityz_InitialCondition_g;

    // InitializeConditions for Integrator: '<S213>/Filter'
    multiModeQuad_ROS_X.Filter_CSTATE_h =
      multiModeQuad_ROS_P.PIDVelocityz_InitialConditionFo;

    // SystemInitialize for Atomic SubSystem: '<Root>/Flight mode'
    // SystemInitialize for Enabled SubSystem: '<S3>/Enabled Subsystem'
    // SystemInitialize for Outport: '<S13>/Out1' incorporates:
    //   Inport: '<S13>/In1'

    multiModeQuad_ROS_B.In1_e = multiModeQuad_ROS_P.Out1_Y0_p;

    // End of SystemInitialize for SubSystem: '<S3>/Enabled Subsystem'

    // Start for MATLABSystem: '<S3>/SourceBlock'
    multiModeQuad_ROS_DW.obj_l.matlabCodegenIsDeleted = false;
    multiModeQuad_ROS_DW.obj_l.isInitialized = 1;
    for (i = 0; i < 11; i++) {
      b_zeroDelimTopic[i] = tmp[i];
    }

    b_zeroDelimTopic[11] = '\x00';
    Sub_multiModeQuad_ROS_472.createSubscriber(&b_zeroDelimTopic[0], 51);
    multiModeQuad_ROS_DW.obj_l.isSetupComplete = true;

    // End of Start for MATLABSystem: '<S3>/SourceBlock'
    // End of SystemInitialize for SubSystem: '<Root>/Flight mode'

    // SystemInitialize for Atomic SubSystem: '<Root>/Sub setpoint attitude'
    // SystemInitialize for Enabled SubSystem: '<S6>/Enabled Subsystem'
    // SystemInitialize for Outport: '<S14>/Out1' incorporates:
    //   Inport: '<S14>/In1'

    multiModeQuad_ROS_B.In1 = multiModeQuad_ROS_P.Out1_Y0;

    // End of SystemInitialize for SubSystem: '<S6>/Enabled Subsystem'

    // Start for MATLABSystem: '<S6>/SourceBlock'
    multiModeQuad_ROS_DW.obj_ap.matlabCodegenIsDeleted = false;
    multiModeQuad_ROS_DW.obj_ap.isInitialized = 1;
    for (i = 0; i < 26; i++) {
      multiModeQuad_ROS_B.b_zeroDelimTopic_b[i] = tmp_0[i];
    }

    multiModeQuad_ROS_B.b_zeroDelimTopic_b[26] = '\x00';
    Sub_multiModeQuad_ROS_496.createSubscriber
      (&multiModeQuad_ROS_B.b_zeroDelimTopic_b[0], 51);
    multiModeQuad_ROS_DW.obj_ap.isSetupComplete = true;

    // End of Start for MATLABSystem: '<S6>/SourceBlock'
    // End of SystemInitialize for SubSystem: '<Root>/Sub setpoint attitude'

    // SystemInitialize for Atomic SubSystem: '<Root>/Sub setpoint velocity'
    // SystemInitialize for Enabled SubSystem: '<S9>/Enabled Subsystem'
    // SystemInitialize for Outport: '<S17>/Out1' incorporates:
    //   Inport: '<S17>/In1'

    multiModeQuad_ROS_B.In1_h = multiModeQuad_ROS_P.Out1_Y0_g;

    // End of SystemInitialize for SubSystem: '<S9>/Enabled Subsystem'

    // Start for MATLABSystem: '<S9>/SourceBlock'
    multiModeQuad_ROS_DW.obj_c.matlabCodegenIsDeleted = false;
    multiModeQuad_ROS_DW.obj_c.isInitialized = 1;
    for (i = 0; i < 35; i++) {
      multiModeQuad_ROS_B.b_zeroDelimTopic[i] = tmp_1[i];
    }

    multiModeQuad_ROS_B.b_zeroDelimTopic[35] = '\x00';
    Sub_multiModeQuad_ROS_426.createSubscriber
      (&multiModeQuad_ROS_B.b_zeroDelimTopic[0], 51);
    multiModeQuad_ROS_DW.obj_c.isSetupComplete = true;

    // End of Start for MATLABSystem: '<S9>/SourceBlock'
    // End of SystemInitialize for SubSystem: '<Root>/Sub setpoint velocity'

    // SystemInitialize for Atomic SubSystem: '<Root>/Sub setpoint rate'
    // SystemInitialize for Enabled SubSystem: '<S7>/Enabled Subsystem'
    // SystemInitialize for Outport: '<S15>/Out1' incorporates:
    //   Inport: '<S15>/In1'

    multiModeQuad_ROS_B.In1_o = multiModeQuad_ROS_P.Out1_Y0_k;

    // End of SystemInitialize for SubSystem: '<S7>/Enabled Subsystem'

    // Start for MATLABSystem: '<S7>/SourceBlock'
    multiModeQuad_ROS_DW.obj_b.matlabCodegenIsDeleted = false;
    multiModeQuad_ROS_DW.obj_b.isInitialized = 1;
    for (i = 0; i < 25; i++) {
      b_zeroDelimTopic_0[i] = tmp_2[i];
    }

    b_zeroDelimTopic_0[25] = '\x00';
    Sub_multiModeQuad_ROS_497.createSubscriber(&b_zeroDelimTopic_0[0], 51);
    multiModeQuad_ROS_DW.obj_b.isSetupComplete = true;

    // End of Start for MATLABSystem: '<S7>/SourceBlock'
    // End of SystemInitialize for SubSystem: '<Root>/Sub setpoint rate'

    // SystemInitialize for IfAction SubSystem: '<S20>/Negative Trace'
    // Start for If: '<S53>/Find Maximum Diagonal Value'
    multiModeQuad_ROS_DW.FindMaximumDiagonalValue_Active = -1;

    // End of SystemInitialize for SubSystem: '<S20>/Negative Trace'

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish'
    // Start for MATLABSystem: '<S4>/SinkBlock'
    multiModeQuad_ROS_DW.obj_m.matlabCodegenIsDeleted = false;
    multiModeQuad_ROS_DW.obj_m.isInitialized = 1;

    // End of SystemInitialize for SubSystem: '<Root>/Publish'

    // SystemInitialize for Merge: '<S20>/Merge'
    multiModeQuad_ROS_B.Merge[0] = multiModeQuad_ROS_P.Merge_InitialOutput[0];

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish'
    // Start for MATLABSystem: '<S4>/SinkBlock'
    b_zeroDelimTopic_1[0] = 'p';

    // End of SystemInitialize for SubSystem: '<Root>/Publish'

    // SystemInitialize for Merge: '<S20>/Merge'
    multiModeQuad_ROS_B.Merge[1] = multiModeQuad_ROS_P.Merge_InitialOutput[1];

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish'
    // Start for MATLABSystem: '<S4>/SinkBlock'
    b_zeroDelimTopic_1[1] = 'o';

    // End of SystemInitialize for SubSystem: '<Root>/Publish'

    // SystemInitialize for Merge: '<S20>/Merge'
    multiModeQuad_ROS_B.Merge[2] = multiModeQuad_ROS_P.Merge_InitialOutput[2];

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish'
    // Start for MATLABSystem: '<S4>/SinkBlock'
    b_zeroDelimTopic_1[2] = 's';

    // End of SystemInitialize for SubSystem: '<Root>/Publish'

    // SystemInitialize for Merge: '<S20>/Merge'
    multiModeQuad_ROS_B.Merge[3] = multiModeQuad_ROS_P.Merge_InitialOutput[3];

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish'
    // Start for MATLABSystem: '<S4>/SinkBlock'
    b_zeroDelimTopic_1[3] = 'e';
    b_zeroDelimTopic_1[4] = '\x00';
    Pub_multiModeQuad_ROS_438.createPublisher(&b_zeroDelimTopic_1[0], 1);
    multiModeQuad_ROS_DW.obj_m.isSetupComplete = true;

    // End of SystemInitialize for SubSystem: '<Root>/Publish'

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish1'
    // Start for MATLABSystem: '<S5>/SinkBlock'
    multiModeQuad_ROS_DW.obj_k.matlabCodegenIsDeleted = false;
    multiModeQuad_ROS_DW.obj_k.isInitialized = 1;
    for (i = 0; i < 7; i++) {
      b_zeroDelimTopic_2[i] = tmp_3[i];
    }

    b_zeroDelimTopic_2[7] = '\x00';
    Pub_multiModeQuad_ROS_477.createPublisher(&b_zeroDelimTopic_2[0], 1);
    multiModeQuad_ROS_DW.obj_k.isSetupComplete = true;

    // End of Start for MATLABSystem: '<S5>/SinkBlock'
    // End of SystemInitialize for SubSystem: '<Root>/Publish1'

    // SystemInitialize for Atomic SubSystem: '<Root>/Sub setpoint thrust'
    // SystemInitialize for Enabled SubSystem: '<S8>/Enabled Subsystem'
    // SystemInitialize for Outport: '<S16>/Out1' incorporates:
    //   Inport: '<S16>/In1'

    multiModeQuad_ROS_B.In1_ol = multiModeQuad_ROS_P.Out1_Y0_kv;

    // End of SystemInitialize for SubSystem: '<S8>/Enabled Subsystem'

    // Start for MATLABSystem: '<S8>/SourceBlock'
    multiModeQuad_ROS_DW.obj_ps.matlabCodegenIsDeleted = false;
    multiModeQuad_ROS_DW.obj_ps.isInitialized = 1;
    for (i = 0; i < 6; i++) {
      b_zeroDelimTopic_3[i] = tmp_4[i];
    }

    b_zeroDelimTopic_3[6] = '\x00';
    Sub_multiModeQuad_ROS_500.createSubscriber(&b_zeroDelimTopic_3[0], 51);
    multiModeQuad_ROS_DW.obj_ps.isSetupComplete = true;

    // End of Start for MATLABSystem: '<S8>/SourceBlock'
    // End of SystemInitialize for SubSystem: '<Root>/Sub setpoint thrust'

    // Start for MATLABSystem: '<Root>/Get Parameter'
    multiModeQuad_ROS_DW.obj_p.matlabCodegenIsDeleted = false;
    multiModeQuad_ROS_DW.obj_p.isInitialized = 1;
    for (i = 0; i < 11; i++) {
      b_zeroDelimTopic[i] = tmp_5[i];
    }

    b_zeroDelimTopic[11] = '\x00';
    ParamGet_multiModeQuad_ROS_459.initialize(&b_zeroDelimTopic[0]);
    ParamGet_multiModeQuad_ROS_459.initialize_error_codes(0, 1, 2, 3);
    ParamGet_multiModeQuad_ROS_459.set_initial_value(0.0);
    multiModeQuad_ROS_DW.obj_p.isSetupComplete = true;

    // End of Start for MATLABSystem: '<Root>/Get Parameter'

    // Start for MATLABSystem: '<Root>/Get Parameter1'
    multiModeQuad_ROS_DW.obj_a.matlabCodegenIsDeleted = false;
    multiModeQuad_ROS_DW.obj_a.isInitialized = 1;
    for (i = 0; i < 11; i++) {
      b_zeroDelimTopic[i] = tmp_6[i];
    }

    b_zeroDelimTopic[11] = '\x00';
    ParamGet_multiModeQuad_ROS_460.initialize(&b_zeroDelimTopic[0]);
    ParamGet_multiModeQuad_ROS_460.initialize_error_codes(0, 1, 2, 3);
    ParamGet_multiModeQuad_ROS_460.set_initial_value(0.0);
    multiModeQuad_ROS_DW.obj_a.isSetupComplete = true;

    // End of Start for MATLABSystem: '<Root>/Get Parameter1'

    // Start for MATLABSystem: '<Root>/Get Parameter2'
    multiModeQuad_ROS_DW.obj.matlabCodegenIsDeleted = false;
    multiModeQuad_ROS_DW.obj.isInitialized = 1;
    for (i = 0; i < 11; i++) {
      b_zeroDelimTopic[i] = tmp_7[i];
    }

    b_zeroDelimTopic[11] = '\x00';
    ParamGet_multiModeQuad_ROS_461.initialize(&b_zeroDelimTopic[0]);
    ParamGet_multiModeQuad_ROS_461.initialize_error_codes(0, 1, 2, 3);
    ParamGet_multiModeQuad_ROS_461.set_initial_value(0.0);
    multiModeQuad_ROS_DW.obj.isSetupComplete = true;

    // End of Start for MATLABSystem: '<Root>/Get Parameter2'

    // InitializeConditions for root-level periodic continuous states
    {
      int_T rootPeriodicContStateIndices[3] = { 0, 1, 2 };

      real_T rootPeriodicContStateRanges[6] = { -3.1415926535897931,
        3.1415926535897931, -3.1415926535897931, 3.1415926535897931,
        -3.1415926535897931, 3.1415926535897931 };

      (void) memcpy((void*)multiModeQuad_ROS_PeriodicIndX,
                    rootPeriodicContStateIndices,
                    3*sizeof(int_T));
      (void) memcpy((void*)multiModeQuad_ROS_PeriodicRngX,
                    rootPeriodicContStateRanges,
                    6*sizeof(real_T));
    }
  }
}

// Model terminate function
void multiModeQuad_ROS_terminate(void)
{
  // Terminate for Atomic SubSystem: '<Root>/Flight mode'
  // Terminate for MATLABSystem: '<S3>/SourceBlock'
  if (!multiModeQuad_ROS_DW.obj_l.matlabCodegenIsDeleted) {
    multiModeQuad_ROS_DW.obj_l.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<S3>/SourceBlock'
  // End of Terminate for SubSystem: '<Root>/Flight mode'

  // Terminate for Atomic SubSystem: '<Root>/Sub setpoint attitude'
  // Terminate for MATLABSystem: '<S6>/SourceBlock'
  if (!multiModeQuad_ROS_DW.obj_ap.matlabCodegenIsDeleted) {
    multiModeQuad_ROS_DW.obj_ap.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<S6>/SourceBlock'
  // End of Terminate for SubSystem: '<Root>/Sub setpoint attitude'

  // Terminate for Atomic SubSystem: '<Root>/Sub setpoint velocity'
  // Terminate for MATLABSystem: '<S9>/SourceBlock'
  if (!multiModeQuad_ROS_DW.obj_c.matlabCodegenIsDeleted) {
    multiModeQuad_ROS_DW.obj_c.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<S9>/SourceBlock'
  // End of Terminate for SubSystem: '<Root>/Sub setpoint velocity'

  // Terminate for Atomic SubSystem: '<Root>/Sub setpoint rate'
  // Terminate for MATLABSystem: '<S7>/SourceBlock'
  if (!multiModeQuad_ROS_DW.obj_b.matlabCodegenIsDeleted) {
    multiModeQuad_ROS_DW.obj_b.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<S7>/SourceBlock'
  // End of Terminate for SubSystem: '<Root>/Sub setpoint rate'

  // Terminate for MATLABSystem: '<Root>/Get Parameter'
  if (!multiModeQuad_ROS_DW.obj_p.matlabCodegenIsDeleted) {
    multiModeQuad_ROS_DW.obj_p.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<Root>/Get Parameter'

  // Terminate for MATLABSystem: '<Root>/Get Parameter1'
  if (!multiModeQuad_ROS_DW.obj_a.matlabCodegenIsDeleted) {
    multiModeQuad_ROS_DW.obj_a.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<Root>/Get Parameter1'

  // Terminate for MATLABSystem: '<Root>/Get Parameter2'
  if (!multiModeQuad_ROS_DW.obj.matlabCodegenIsDeleted) {
    multiModeQuad_ROS_DW.obj.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<Root>/Get Parameter2'

  // Terminate for Atomic SubSystem: '<Root>/Publish'
  // Terminate for MATLABSystem: '<S4>/SinkBlock'
  if (!multiModeQuad_ROS_DW.obj_m.matlabCodegenIsDeleted) {
    multiModeQuad_ROS_DW.obj_m.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<S4>/SinkBlock'
  // End of Terminate for SubSystem: '<Root>/Publish'

  // Terminate for Atomic SubSystem: '<Root>/Publish1'
  // Terminate for MATLABSystem: '<S5>/SinkBlock'
  if (!multiModeQuad_ROS_DW.obj_k.matlabCodegenIsDeleted) {
    multiModeQuad_ROS_DW.obj_k.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<S5>/SinkBlock'
  // End of Terminate for SubSystem: '<Root>/Publish1'

  // Terminate for Atomic SubSystem: '<Root>/Sub setpoint thrust'
  // Terminate for MATLABSystem: '<S8>/SourceBlock'
  if (!multiModeQuad_ROS_DW.obj_ps.matlabCodegenIsDeleted) {
    multiModeQuad_ROS_DW.obj_ps.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<S8>/SourceBlock'
  // End of Terminate for SubSystem: '<Root>/Sub setpoint thrust'
}

//
// File trailer for generated code.
//
// [EOF]
//
