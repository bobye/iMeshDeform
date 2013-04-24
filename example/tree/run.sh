#!/bin/sh

# run main program
../../bin/iMeshDeform \
	tree.off \
	tree_c60.id \
	tree_c60.id \
	tree_empty.id \
	-pc_factor_mat_solver_package umfpack

