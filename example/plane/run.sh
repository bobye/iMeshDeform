#!/bin/sh

# run main program
../../bin/iMeshDeform plane.off plane_linear25_cluster.id plane_rot25_cluster.id plane_affine.id -pc_factor_mat_solver_package umfpack

