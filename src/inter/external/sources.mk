inter_external_SRCS=\
    $(SRC_DIR)/inter/external/BoxExternal.cpp \
    $(SRC_DIR)/inter/external/OrthoBoxExternal.cpp \
    $(SRC_DIR)/inter/external/SlitExternal.cpp \
    $(SRC_DIR)/inter/external/LamellarOrderingExternal.cpp \
    $(SRC_DIR)/inter/external/LocalLamellarOrderingExternal.cpp \
    $(SRC_DIR)/inter/external/PeriodicExternal.cpp 

inter_external_OBJS=$(inter_external_SRCS:.cpp=.o)

