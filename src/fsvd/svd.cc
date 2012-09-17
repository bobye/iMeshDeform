#include "fsvd.hh"
#include <cmath>
#include <algorithm>
#define USE_SCALAR_IMPLEMENTATION

#define COMPUTE_V_AS_MATRIX
#define COMPUTE_U_AS_MATRIX

#include "Singular_Value_Decomposition_Preamble.hpp"

void fastsvd(float* const A, float * const U, float * const V, float * const Sigma) {
#include "Singular_Value_Decomposition_Kernel_Declarations.hpp"
  ENABLE_SCALAR_IMPLEMENTATION(Sa11.f=A[0];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa21.f=A[1];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa31.f=A[2];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa12.f=A[3];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa22.f=A[4];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa32.f=A[5];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa13.f=A[6];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa23.f=A[7];)
  ENABLE_SCALAR_IMPLEMENTATION(Sa33.f=A[8];)
#include "Singular_Value_Decomposition_Main_Kernel_Body.hpp"

  Sigma[0]=Sa11.f;Sigma[1]=Sa22.f;Sigma[2]=Sa33.f;
  U[0]=Su11.f;U[1]=Su21.f;U[2]=Su31.f;
  U[3]=Su12.f;U[4]=Su22.f;U[5]=Su32.f;
  U[6]=Su13.f;U[7]=Su23.f;U[8]=Su33.f;
  V[0]=Sv11.f;V[1]=Sv21.f;V[2]=Sv31.f;
  V[3]=Sv12.f;V[4]=Sv22.f;V[5]=Sv32.f;
  V[6]=Sv13.f;V[7]=Sv23.f;V[8]=Sv33.f;
}
