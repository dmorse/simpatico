simp_interaction_external_=\
    simp/interaction/external/BoxExternal.cpp \
    simp/interaction/external/OrthoBoxExternal.cpp \
    simp/interaction/external/SlitExternal.cpp \
    simp/interaction/external/LamellarOrderingExternal.cpp \
    simp/interaction/external/LocalLamellarOrderingExternal.cpp \
    simp/interaction/external/SimplePeriodicExternal.cpp \
    simp/interaction/external/GeneralPeriodicExternal.cpp \
    simp/interaction/external/NucleationExternal.cpp \
    simp/interaction/external/PeriodicExternal.cpp \
    simp/interaction/external/SphericalTabulatedExternal.cpp 

simp_interaction_external_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_interaction_external_))
simp_interaction_external_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_interaction_external_:.cpp=.o))

