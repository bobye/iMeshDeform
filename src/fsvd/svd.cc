#include "fsvd.hh"
#include <cmath>
#include <algorithm>
#include <mkl.h>

#define USE_SCALAR_IMPLEMENTATION

#define COMPUTE_V_AS_MATRIX
#define COMPUTE_U_AS_MATRIX

#include "Singular_Value_Decomposition_Preamble.hpp"

void fastsvd(double* const A, float * const U, float * const V, float * const Sigma) {
#include "Singular_Value_Decomposition_Kernel_Declarations.hpp"
  ENABLE_SCALAR_IMPLEMENTATION(Sa11.f=(float) A[0];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa21.f=(float) A[1];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa31.f=(float) A[2];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa12.f=(float) A[3];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa22.f=(float) A[4];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa32.f=(float) A[5];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa13.f=(float) A[6];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa23.f=(float) A[7];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa33.f=(float) A[8];)
#include "Singular_Value_Decomposition_Main_Kernel_Body.hpp"

  Sigma[0]=Sa11.f;Sigma[1]=Sa22.f;Sigma[2]=Sa33.f;
  U[0]=Su11.f;U[1]=Su21.f;U[2]=Su31.f;
  U[3]=Su12.f;U[4]=Su22.f;U[5]=Su32.f;
  U[6]=Su13.f;U[7]=Su23.f;U[8]=Su33.f;
  V[0]=Sv11.f;V[1]=Sv21.f;V[2]=Sv31.f;
  V[3]=Sv12.f;V[4]=Sv22.f;V[5]=Sv32.f;
  V[6]=Sv13.f;V[7]=Sv23.f;V[8]=Sv33.f;
}


void dfastsvd(double* const A, const int n) {
  for (int i=0; i < 9*n; i+=9) {
    float U[9], V[9], S[3];
    double *R = A+i;
    double norm = cblas_dnrm2(9, R, 1);
    cblas_dscal(9, 1/norm, R, 1);
    
    fastsvd(R, U, V, S);
    R[0] = V[0]*U[0] + V[3]*U[3] + V[6]*U[6];
    R[1] = V[1]*U[0] + V[4]*U[3] + V[7]*U[6];
    R[2] = V[2]*U[0] + V[5]*U[3] + V[8]*U[6];
    R[3] = V[0]*U[1] + V[3]*U[4] + V[6]*U[7];
    R[4] = V[1]*U[1] + V[4]*U[4] + V[7]*U[7];
    R[5] = V[2]*U[1] + V[5]*U[4] + V[8]*U[7];
    R[6] = V[0]*U[2] + V[3]*U[5] + V[6]*U[8];
    R[7] = V[1]*U[2] + V[4]*U[5] + V[7]*U[8];
    R[8] = V[2]*U[2] + V[5]*U[5] + V[8]*U[8];

  }
}
