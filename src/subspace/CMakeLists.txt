SET(SUBSPACE_LIB_SRC gui.cc Object.cc Select.cc subspace.cc PTHREAD_QUEUE.cpp)

ADD_LIBRARY(subspace STATIC ${SUBSPACE_LIB_SRC})
TARGET_LINK_LIBRARIES(subspace ${LIBS})
