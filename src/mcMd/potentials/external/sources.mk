mcMd_potentials_external_SRCS=\
    $(SRC_DIR)/mcMd/potentials/external/BoxExternal.cpp \
    $(SRC_DIR)/mcMd/potentials/external/OrthoBoxExternal.cpp \
    $(SRC_DIR)/mcMd/potentials/external/SlitExternal.cpp \
    $(SRC_DIR)/mcMd/potentials/external/TanhCosineExternal.cpp \
    $(SRC_DIR)/mcMd/potentials/external/ExternalFactory.cpp 

mcMd_potentials_external_OBJS=$(mcMd_potentials_external_SRCS:.cpp=.o)

