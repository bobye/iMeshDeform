#define USE_SCALAR_IMPLEMENTATION
#define USE_SSE_IMPLEMENTATION

#define COMPUTE_V_AS_MATRIX
#define COMPUTE_U_AS_MATRIX

#include "Singular_Value_Decomposition_Preamble.hpp"
template <class T>
inline void fsvd_unit(
        T* a11,T* a21,T* a31,T* a12,T* a22,T* a32,T* a13,T* a23,T* a33,
        T* u11,T* u21,T* u31,T* u12,T* u22,T* u32,T* u13,T* u23,T* u33,
        T* v11,T* v21,T* v31,T* v12,T* v22,T* v32,T* v13,T* v23,T* v33,
        T* sigma1,T* sigma2,T* sigma3) {
#include "Singular_Value_Decomposition_Kernel_Declarations.hpp"
  ENABLE_SCALAR_IMPLEMENTATION(Sa11.f=*a11;)
    ENABLE_SCALAR_IMPLEMENTATION(Sa21.f=*a21;)
    ENABLE_SCALAR_IMPLEMENTATION(Sa31.f=*a31;)
    ENABLE_SCALAR_IMPLEMENTATION(Sa12.f=*a12;)
    ENABLE_SCALAR_IMPLEMENTATION(Sa22.f=*a22;)
    ENABLE_SCALAR_IMPLEMENTATION(Sa32.f=*a32;)
    ENABLE_SCALAR_IMPLEMENTATION(Sa13.f=*a13;)
    ENABLE_SCALAR_IMPLEMENTATION(Sa23.f=*a23;)
    ENABLE_SCALAR_IMPLEMENTATION(Sa33.f=*a33;)

#include "Singular_Value_Decomposition_Main_Kernel_Body.hpp"

    ENABLE_SCALAR_IMPLEMENTATION(*sigma1=Sa11.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*sigma2=Sa22.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*sigma3=Sa33.f;)

    ENABLE_SCALAR_IMPLEMENTATION(*u11=Su11.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*u21=Su21.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*u31=Su31.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*u12=Su12.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*u22=Su22.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*u32=Su32.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*u13=Su13.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*u23=Su23.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*u33=Su33.f;)

    ENABLE_SCALAR_IMPLEMENTATION(*v11=Sv11.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*v21=Sv21.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*v31=Sv31.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*v12=Sv12.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*v22=Sv22.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*v32=Sv32.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*v13=Sv13.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*v23=Sv23.f;)
    ENABLE_SCALAR_IMPLEMENTATION(*v33=Sv33.f;)

    //    printf("%.3f %.3f %.3f\n", *sigma1, *sigma2, *sigma3);
    //printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", *u11, *u21, *u31, *u12, *u22, *u32, *u13, *u23, *u33);
    //printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", *v11, *v21, *v31, *v12, *v22, *v32, *v13, *v23, *v33);
};
#undef USE_SCALAR_IMPLEMENTATION // #undef USE_SSE_IMPLEMENTATION
/***************************************************************************/


#include <algorithm>
#include "PTHREAD_QUEUE.h"
//#include "Singular_Value_Decomposition_Helper.h"

float *a11,*a21,*a31,*a12,*a22,*a32,*a13,*a23,*a33;
float *u11,*u21,*u31,*u12,*u22,*u32,*u13,*u23,*u33;
float *v11,*v21,*v31,*v12,*v22,*v32,*v13,*v23,*v33;
float *sigma1,*sigma2,*sigma3;

int size, number_of_threads;//the number of size should be 8x

using namespace PhysBAM;


extern PTHREAD_QUEUE* pthread_queue;


template<class T>
class Singular_Value_Decomposition_Size_Specific_Helper
{
    T* const a11,* const a21,* const a31,* const a12,* const a22,* const a32,* const a13,* const a23,* const a33;
    T* const u11,* const u21,* const u31,* const u12,* const u22,* const u32,* const u13,* const u23,* const u33;
    T* const v11,* const v21,* const v31,* const v12,* const v22,* const v32,* const v13,* const v23,* const v33;
    T* const sigma1,* const sigma2,* const sigma3;

  int size;
  int number_of_partitions;

public:
  explicit Singular_Value_Decomposition_Size_Specific_Helper (
        const int size_input, const int number_of_partitions_input,
        T* const a11_input,T* const a21_input,T* const a31_input,
        T* const a12_input,T* const a22_input,T* const a32_input,
        T* const a13_input,T* const a23_input,T* const a33_input,
        T* const u11_input,T* const u21_input,T* const u31_input,
        T* const u12_input,T* const u22_input,T* const u32_input,
        T* const u13_input,T* const u23_input,T* const u33_input,
        T* const v11_input,T* const v21_input,T* const v31_input,
        T* const v12_input,T* const v22_input,T* const v32_input,
        T* const v13_input,T* const v23_input,T* const v33_input,
        T* const sigma1_input,T* const sigma2_input,T* const sigma3_input)
        : size(size_input), number_of_partitions(number_of_partitions_input),
        a11(a11_input),a21(a21_input),a31(a31_input),a12(a12_input),a22(a22_input),a32(a32_input),a13(a13_input),a23(a23_input),a33(a33_input),
        u11(u11_input),u21(u21_input),u31(u31_input),u12(u12_input),u22(u22_input),u32(u32_input),u13(u13_input),u23(u23_input),u33(u33_input),
        v11(v11_input),v21(v21_input),v31(v31_input),v12(v12_input),v22(v22_input),v32(v32_input),v13(v13_input),v23(v23_input),v33(v33_input),
        sigma1(sigma1_input),sigma2(sigma2_input),sigma3(sigma3_input)
    {}

    void Run()
  {Run_Index_Range(0, size);}
  
//#####################################################################
    static void Allocate_Data(
        T*& a11,T*& a21,T*& a31,T*& a12,T*& a22,T*& a32,T*& a13,T*& a23,T*& a33,
        T*& u11,T*& u21,T*& u31,T*& u12,T*& u22,T*& u32,T*& u13,T*& u23,T*& u33,
        T*& v11,T*& v21,T*& v31,T*& v12,T*& v22,T*& v32,T*& v13,T*& v23,T*& v33,
        T*& sigma1,T*& sigma2,T*& sigma3);
    static void Initialize_Data(
        T*& a11,T*& a21,T*& a31,T*& a12,T*& a22,T*& a32,T*& a13,T*& a23,T*& a33,
        T*& u11,T*& u21,T*& u31,T*& u12,T*& u22,T*& u32,T*& u13,T*& u23,T*& u33,
        T*& v11,T*& v21,T*& v31,T*& v12,T*& v22,T*& v32,T*& v13,T*& v23,T*& v33,
        T*& sigma1,T*& sigma2,T*& sigma3);
    void Run_Parallel();
    void Run_Index_Range(const int imin, const int imax_plus_one);
//#####################################################################
};


template<class T>
struct Singular_Value_Decomposition_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Singular_Value_Decomposition_Size_Specific_Helper<T>* const obj;
    const int imin,imax_plus_one;
    Singular_Value_Decomposition_Size_Specific_Thread_Helper(
        Singular_Value_Decomposition_Size_Specific_Helper<T>* const obj_input,const int imin_input,const int imax_plus_one_input)
        :obj(obj_input),imin(imin_input),imax_plus_one(imax_plus_one_input) {}
    void Run(){obj->Run_Index_Range(imin,imax_plus_one);}
};

template<class T> void Singular_Value_Decomposition_Size_Specific_Helper<T>::
Run_Parallel()
{
    for(int partition=0;partition<number_of_partitions;partition++){
        int imin=(size/number_of_partitions)*partition+std::min(size%number_of_partitions,partition);
        int imax_plus_one=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1);
 	Singular_Value_Decomposition_Size_Specific_Thread_Helper<T>* task=new Singular_Value_Decomposition_Size_Specific_Thread_Helper<T>(this,imin,imax_plus_one);
 	pthread_queue->Queue(task);
    }
    pthread_queue->Wait();    
}

//#####################################################################
// Function Run_Index_Range
//#####################################################################
template<class T> void Singular_Value_Decomposition_Size_Specific_Helper<T>::
Run_Index_Range(const int imin,const int imax_plus_one)
{   

#include "Singular_Value_Decomposition_Kernel_Declarations.hpp"

#ifdef USE_SSE_IMPLEMENTATION
#define STEP_SIZE 4
#endif
#ifdef USE_AVX_IMPLEMENTATION
#define STEP_SIZE 8
#endif
#ifdef USE_SCALAR_IMPLEMENTATION
#define STEP_SIZE 1
#endif

    for(int index=imin;index<imax_plus_one;index+=STEP_SIZE){
        
    ENABLE_SCALAR_IMPLEMENTATION(Sa11.f=a11[index];)                                      ENABLE_SSE_IMPLEMENTATION(Va11=_mm_loadu_ps(a11+index);)                                  ENABLE_AVX_IMPLEMENTATION(Va11=_mm256_loadu_ps(a11+index);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa21.f=a21[index];)                                      ENABLE_SSE_IMPLEMENTATION(Va21=_mm_loadu_ps(a21+index);)                                  ENABLE_AVX_IMPLEMENTATION(Va21=_mm256_loadu_ps(a21+index);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa31.f=a31[index];)                                      ENABLE_SSE_IMPLEMENTATION(Va31=_mm_loadu_ps(a31+index);)                                  ENABLE_AVX_IMPLEMENTATION(Va31=_mm256_loadu_ps(a31+index);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa12.f=a12[index];)                                      ENABLE_SSE_IMPLEMENTATION(Va12=_mm_loadu_ps(a12+index);)                                  ENABLE_AVX_IMPLEMENTATION(Va12=_mm256_loadu_ps(a12+index);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa22.f=a22[index];)                                      ENABLE_SSE_IMPLEMENTATION(Va22=_mm_loadu_ps(a22+index);)                                  ENABLE_AVX_IMPLEMENTATION(Va22=_mm256_loadu_ps(a22+index);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa32.f=a32[index];)                                      ENABLE_SSE_IMPLEMENTATION(Va32=_mm_loadu_ps(a32+index);)                                  ENABLE_AVX_IMPLEMENTATION(Va32=_mm256_loadu_ps(a32+index);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa13.f=a13[index];)                                      ENABLE_SSE_IMPLEMENTATION(Va13=_mm_loadu_ps(a13+index);)                                  ENABLE_AVX_IMPLEMENTATION(Va13=_mm256_loadu_ps(a13+index);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa23.f=a23[index];)                                      ENABLE_SSE_IMPLEMENTATION(Va23=_mm_loadu_ps(a23+index);)                                  ENABLE_AVX_IMPLEMENTATION(Va23=_mm256_loadu_ps(a23+index);)
    ENABLE_SCALAR_IMPLEMENTATION(Sa33.f=a33[index];)                                      ENABLE_SSE_IMPLEMENTATION(Va33=_mm_loadu_ps(a33+index);)                                  ENABLE_AVX_IMPLEMENTATION(Va33=_mm256_loadu_ps(a33+index);)

#include "Singular_Value_Decomposition_Main_Kernel_Body.hpp"

    ENABLE_SCALAR_IMPLEMENTATION(u11[index]=Su11.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(u11+index,Vu11);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(u11+index,Vu11);)
    ENABLE_SCALAR_IMPLEMENTATION(u21[index]=Su21.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(u21+index,Vu21);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(u21+index,Vu21);)
    ENABLE_SCALAR_IMPLEMENTATION(u31[index]=Su31.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(u31+index,Vu31);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(u31+index,Vu31);)
    ENABLE_SCALAR_IMPLEMENTATION(u12[index]=Su12.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(u12+index,Vu12);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(u12+index,Vu12);)
    ENABLE_SCALAR_IMPLEMENTATION(u22[index]=Su22.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(u22+index,Vu22);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(u22+index,Vu22);)
    ENABLE_SCALAR_IMPLEMENTATION(u32[index]=Su32.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(u32+index,Vu32);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(u32+index,Vu32);)
    ENABLE_SCALAR_IMPLEMENTATION(u13[index]=Su13.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(u13+index,Vu13);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(u13+index,Vu13);)
    ENABLE_SCALAR_IMPLEMENTATION(u23[index]=Su23.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(u23+index,Vu23);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(u23+index,Vu23);)
    ENABLE_SCALAR_IMPLEMENTATION(u33[index]=Su33.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(u33+index,Vu33);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(u33+index,Vu33);)

    ENABLE_SCALAR_IMPLEMENTATION(v11[index]=Sv11.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(v11+index,Vv11);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(v11+index,Vv11);)
    ENABLE_SCALAR_IMPLEMENTATION(v21[index]=Sv21.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(v21+index,Vv21);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(v21+index,Vv21);)
    ENABLE_SCALAR_IMPLEMENTATION(v31[index]=Sv31.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(v31+index,Vv31);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(v31+index,Vv31);)
    ENABLE_SCALAR_IMPLEMENTATION(v12[index]=Sv12.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(v12+index,Vv12);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(v12+index,Vv12);)
    ENABLE_SCALAR_IMPLEMENTATION(v22[index]=Sv22.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(v22+index,Vv22);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(v22+index,Vv22);)
    ENABLE_SCALAR_IMPLEMENTATION(v32[index]=Sv32.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(v32+index,Vv32);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(v32+index,Vv32);)
    ENABLE_SCALAR_IMPLEMENTATION(v13[index]=Sv13.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(v13+index,Vv13);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(v13+index,Vv13);)
    ENABLE_SCALAR_IMPLEMENTATION(v23[index]=Sv23.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(v23+index,Vv23);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(v23+index,Vv23);)
    ENABLE_SCALAR_IMPLEMENTATION(v33[index]=Sv33.f;)                                      ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(v33+index,Vv33);)                                 ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(v33+index,Vv33);)

    ENABLE_SCALAR_IMPLEMENTATION(sigma1[index]=Sa11.f;)                                   ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(sigma1+index,Va11);)                              ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(sigma1+index,Va11);)
    ENABLE_SCALAR_IMPLEMENTATION(sigma2[index]=Sa22.f;)                                   ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(sigma2+index,Va22);)                              ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(sigma2+index,Va22);)
    ENABLE_SCALAR_IMPLEMENTATION(sigma3[index]=Sa33.f;)                                   ENABLE_SSE_IMPLEMENTATION(_mm_storeu_ps(sigma3+index,Va33);)                              ENABLE_AVX_IMPLEMENTATION(_mm256_storeu_ps(sigma3+index,Va33);)

    }

#undef STEP_SIZE

}

//#####################################################################

inline void init_svd(_SS_SCALAR* const A11, 
		     _SS_SCALAR* const A12,
		     _SS_SCALAR* const A13,
		     _SS_SCALAR* const A21,
		     _SS_SCALAR* const A22,
		     _SS_SCALAR* const A23,
		     _SS_SCALAR* const A31,
		     _SS_SCALAR* const A32,
		     _SS_SCALAR* const A33,
		     const int size_t, const int number_of_threads_t) {
  size = size_t;
  number_of_threads = number_of_threads_t;

  //printf("Using %d threads\n",number_of_threads);
  pthread_queue=new PhysBAM::PTHREAD_QUEUE(number_of_threads);  

    a11=new float[size];
    a21=new float[size];
    a31=new float[size];
    a12=new float[size];
    a22=new float[size];
    a32=new float[size];
    a13=new float[size];
    a23=new float[size];
    a33=new float[size];

    u11=new float[size];
    u21=new float[size];
    u31=new float[size];
    u12=new float[size];
    u22=new float[size];
    u32=new float[size];
    u13=new float[size];
    u23=new float[size];
    u33=new float[size];

    v11=new float[size];
    v21=new float[size];
    v31=new float[size];
    v12=new float[size];
    v22=new float[size];
    v32=new float[size];
    v13=new float[size];
    v23=new float[size];
    v33=new float[size];

    sigma1=new float[size];
    sigma2=new float[size];
    sigma3=new float[size];

}
inline void destroy_svd() {// to implement
}

inline void proj_rot(_SS_SCALAR* const A11, 
		     _SS_SCALAR* const A12,
		     _SS_SCALAR* const A13,
		     _SS_SCALAR* const A21,
		     _SS_SCALAR* const A22,
		     _SS_SCALAR* const A23,
		     _SS_SCALAR* const A31,
		     _SS_SCALAR* const A32,
		     _SS_SCALAR* const A33,
		     const int size) {

    for (int i=0; i<size; ++i) 
    { a11[i]=A11[i]; a21[i]=A12[i]; a31[i]=A13[i]; a12[i]=A21[i]; a22[i]=A22[i]; a32[i]=A23[i]; a13[i]=A31[i]; a23[i]=A32[i]; a33[i]=A33[i];}

    
    for (int i=0; i<size; ++i) {
      float one_over_frobenius_norm=(float)1./sqrt(
        (double)a11[i]*(double)a11[i]+(double)a12[i]*(double)a12[i]+(double)a13[i]*(double)a13[i]+
        (double)a21[i]*(double)a21[i]+(double)a22[i]*(double)a22[i]+(double)a23[i]*(double)a23[i]+
        (double)a31[i]*(double)a31[i]+(double)a32[i]*(double)a32[i]+(double)a33[i]*(double)a33[i]);

      a11[i]*=one_over_frobenius_norm;
      a12[i]*=one_over_frobenius_norm;
      a13[i]*=one_over_frobenius_norm;
      a21[i]*=one_over_frobenius_norm;
      a22[i]*=one_over_frobenius_norm;
      a23[i]*=one_over_frobenius_norm;
      a31[i]*=one_over_frobenius_norm;
      a32[i]*=one_over_frobenius_norm;
      a33[i]*=one_over_frobenius_norm;
    }

    int size4 = 4*number_of_threads*(size/(4*number_of_threads));
    int sizeu = size%(4*number_of_threads);
    if (size4!=0) {
    Singular_Value_Decomposition_Size_Specific_Helper<float> test(size4, number_of_threads,
        a11,a21,a31,a12,a22,a32,a13,a23,a33,
        u11,u21,u31,u12,u22,u32,u13,u23,u33,
        v11,v21,v31,v12,v22,v32,v13,v23,v33,
        sigma1,sigma2,sigma3);
    test.Run_Parallel();
    }

    for (int i=0; i<sizeu; ++i) {
      int j = size4+i;
      fsvd_unit<float>(
        a11+j,a21+j,a31+j,a12+j,a22+j,a32+j,a13+j,a23+j,a33+j,
        u11+j,u21+j,u31+j,u12+j,u22+j,u32+j,u13+j,u23+j,u33+j,
        v11+j,v21+j,v31+j,v12+j,v22+j,v32+j,v13+j,v23+j,v33+j,
        sigma1+j,sigma2+j,sigma3+j);
    }


    for (int i=0; i<size; ++i) {
      A11[i] = v11[i]*u11[i] + v12[i]*u12[i] + v13[i]*u13[i];
      A12[i] = v21[i]*u11[i] + v22[i]*u12[i] + v23[i]*u13[i];
      A13[i] = v31[i]*u11[i] + v32[i]*u12[i] + v33[i]*u13[i];

      A21[i] = v11[i]*u21[i] + v12[i]*u22[i] + v13[i]*u23[i];
      A22[i] = v21[i]*u21[i] + v22[i]*u22[i] + v23[i]*u23[i];
      A23[i] = v31[i]*u21[i] + v32[i]*u22[i] + v33[i]*u23[i];

      A31[i] = v11[i]*u31[i] + v12[i]*u32[i] + v13[i]*u33[i];
      A32[i] = v21[i]*u31[i] + v22[i]*u32[i] + v23[i]*u33[i];
      A33[i] = v31[i]*u31[i] + v32[i]*u32[i] + v33[i]*u33[i];
    }
    /*
    R[0] = V[0]*U[0] + V[3]*U[3] + V[6]*U[6];
    R[1] = V[1]*U[0] + V[4]*U[3] + V[7]*U[6];
    R[2] = V[2]*U[0] + V[5]*U[3] + V[8]*U[6];
    R[3] = V[0]*U[1] + V[3]*U[4] + V[6]*U[7];
    R[4] = V[1]*U[1] + V[4]*U[4] + V[7]*U[7];
    R[5] = V[2]*U[1] + V[5]*U[4] + V[8]*U[7];
    R[6] = V[0]*U[2] + V[3]*U[5] + V[6]*U[8];
    R[7] = V[1]*U[2] + V[4]*U[5] + V[7]*U[8];
    R[8] = V[2]*U[2] + V[5]*U[5] + V[8]*U[8];
    */
}

inline void init_svd(_SS_SCALAR* A, const int n, const int npthread) {
  init_svd(A, A+n, A+2*n, A+3*n, A+4*n, A+5*n, A+6*n, A+7*n, A+8*n, n, npthread);
}

inline void proj_rot(_SS_SCALAR* A, const int n) {
  proj_rot(A, A+n, A+2*n, A+3*n, A+4*n, A+5*n, A+6*n, A+7*n, A+8*n, n);
}
