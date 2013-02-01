util_ensembles_SRCS=\
    $(SRC_DIR)/util/ensembles/EnergyEnsemble.cpp \
    $(SRC_DIR)/util/ensembles/BoundaryEnsemble.cpp \
    $(SRC_DIR)/util/ensembles/SpeciesEnsemble.cpp 

util_ensembles_OBJS=$(util_ensembles_SRCS:.cpp=.o)

