
simp_boundary_=\
    simp/boundary/OrthoRegion.cpp \
    simp/boundary/OrthorhombicBoundary.cpp \
    simp/boundary/MonoclinicBoundary.cpp 
    #simp/boundary/MonoclinicBoundaryMI.cpp

simp_boundary_SRCS=$(addprefix $(SRC_DIR)/, $(simp_boundary_))
simp_boundary_OBJS=$(addprefix $(BLD_DIR)/, $(simp_boundary_:.cpp=.o))

