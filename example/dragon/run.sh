#!/bin/sh

# run main program
../../bin/iMeshDeform \
	dragon.off \
	dragon_cluster30.id \
	dragon_cluster30.id \
	dragon_empty.id \
	-pc_factor_mat_solver_package umfpack

