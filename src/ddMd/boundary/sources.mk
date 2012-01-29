ddMd_boundary_SRCS=$(SRC_DIR)/ddMd/boundary/CubicBoundary.cpp \
    $(SRC_DIR)/ddMd/boundary/LatticeSystem.cpp \
    $(SRC_DIR)/ddMd/boundary/OrthoBoundaryBase.cpp \
    $(SRC_DIR)/ddMd/boundary/OrthoRegion.cpp \
    $(SRC_DIR)/ddMd/boundary/OrthorhombicBoundary.cpp \
    $(SRC_DIR)/ddMd/boundary/TetragonalBoundary.cpp 

ddMd_boundary_OBJS=$(ddMd_boundary_SRCS:.cpp=.o)

