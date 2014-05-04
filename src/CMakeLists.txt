ADD_EXECUTABLE(iVertSelect main.cc)
TARGET_LINK_LIBRARIES(iVertSelect subspace ${LIBS})

#INCLUDE_DIRECTORIES(${SLEPC_INCLUDES})
#ADD_EXECUTABLE(eigen_solver eigen_solver.cc)
#TARGET_LINK_LIBRARIES(eigen_solver ${SOLVER_LIBS})
