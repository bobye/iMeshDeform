iMeshDeform
===========
Variational Reduced Deformable Model

Prerequisites:

 - Linux
 - [Intel C++ Compiler for Linux](http://software.intel.com/en-us/non-commercial-software-development)
 - CMake

Libraries:

 - OpenGL + GLUT
 - PETSc with external sparse direct solver, e.g. mumps/umfpack
 - SLEPc
 - Intel MKL [doc](http://software.intel.com/en-us/articles/intel-math-kernel-library-documentation)
 - [trimesh2++](https://github.com/bobye/trimesh2plus)

Env variables:

	source /opt/intel/bin/compilervars.sh intel64
	export SLEPC_DIR=/home/bobye/pub/slepc/slepc-3.3-p3
	export PETSC_DIR=/home/bobye/pub/petsc/petsc-3.3-p5
	export PETSC_ARCH=linux-gnu-c-debug
	export TRIMESH2_DIR=/home/bobye/pub/trimesh2++
	export TRIMESH2_ARCH=Linux64
	export MKLPATH=${MKLROOT}/lib/intel64
	PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:${PKG_CONFIG_PATH}
	export PKG_CONFIG_PATH


How to run:
	
	# run main program
	./bin/iMeshDeform test2.off test2_cluster1.id test2_rotation.id test2_rigid.id -pc_factor_mat_solver_package umfpack

	# run pre-configuration
	./bin/preconfig test2.off
