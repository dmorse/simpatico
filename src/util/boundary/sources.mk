
util_boundary_SRCS=\
    $(SRC_DIR)/util/boundary/OrthoRegion.cpp \
    $(SRC_DIR)/util/boundary/OrthorhombicBoundary.cpp \
    $(SRC_DIR)/util/boundary/MonoclinicBoundary.cpp 
    #$(SRC_DIR)/util/boundary/MonoclinicBoundaryMI.cpp

util_boundary_OBJS=$(util_boundary_SRCS:.cpp=.o)

