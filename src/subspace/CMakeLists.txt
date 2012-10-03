SET(SUBSPACE_LIB_SRC gui.cc Object.cc VertSelect.cc subspace.cc)

ADD_LIBRARY(subspace STATIC ${SUBSPACE_LIB_SRC})
TARGET_LINK_LIBRARIES(subspace ${LIBS})
