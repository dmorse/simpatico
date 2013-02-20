inter_external_SRCS=\
    $(SRC_DIR)/inter/external/BoxExternal.cpp \
    $(SRC_DIR)/inter/external/OrthoBoxExternal.cpp \
    $(SRC_DIR)/inter/external/SlitExternal.cpp \
    $(SRC_DIR)/inter/external/TanhCosineExternal.cpp \
    $(SRC_DIR)/inter/external/OrderingExternal.cpp

inter_external_OBJS=$(inter_external_SRCS:.cpp=.o)

