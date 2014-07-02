inter_external_=\
    inter/external/BoxExternal.cpp \
    inter/external/OrthoBoxExternal.cpp \
    inter/external/SlitExternal.cpp \
    inter/external/LamellarOrderingExternal.cpp \
    inter/external/LocalLamellarOrderingExternal.cpp \
    inter/external/SimplePeriodicExternal.cpp \
    inter/external/GeneralPeriodicExternal.cpp \
    inter/external/NucleationExternal.cpp \
    inter/external/PeriodicExternal.cpp

inter_external_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_external_))
inter_external_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_external_:.cpp=.o))

