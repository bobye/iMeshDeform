#include "fsvd.hh"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

int main(int argc,char* argv[])
{
    if(argc!=2){printf("Must specify an integer random seed as argument\n");exit(1);}
    srand(atoi(argv[1]));

    float test_A[9], norm_inverse=0., U[9],V[9],S[3];

    for (int i=0; i<9; ++i) 
      {
	test_A[i]=2.*(float)rand()/(float)RAND_MAX-1.;
	norm_inverse += test_A[i]*test_A[i];
      }
    norm_inverse = 1/std::sqrt(norm_inverse);
    for (int i=0; i<9; ++i) test_A[i]*=norm_inverse;
    
    fastsvd(test_A,U,V,S);
    for (int i=0; i<9; ++i) {
      std::cout << " " << U[i];
      if (i%3==2) std::cout << std::endl; 
    }
    std::cout << std::endl;
    for (int i=0; i<9; ++i) {
      std::cout << " " << V[i];
      if (i%3==2) std::cout << std::endl; 
    }
    std::cout << std::endl;
    for (int i=0; i<3; ++i) {
      std::cout << " " << S[i];
    }
    std::cout << std::endl;

    return 0;
}
//#####################################################################
