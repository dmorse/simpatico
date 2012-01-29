mcMd_ensembles_SRCS=$(SRC_DIR)/mcMd/ensembles/BoundaryEnsemble.cpp \
    $(SRC_DIR)/mcMd/ensembles/EnergyEnsemble.cpp \
    $(SRC_DIR)/mcMd/ensembles/SpeciesEnsemble.cpp 

mcMd_ensembles_OBJS=$(mcMd_ensembles_SRCS:.cpp=.o)

