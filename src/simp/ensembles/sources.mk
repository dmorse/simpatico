simp_ensembles_=\
    simp/ensembles/EnergyEnsemble.cpp \
    simp/ensembles/BoundaryEnsemble.cpp \
    simp/ensembles/SpeciesEnsemble.cpp 

simp_ensembles_SRCS=$(addprefix $(SRC_DIR)/, $(simp_ensembles_))
simp_ensembles_OBJS=$(addprefix $(BLD_DIR)/, $(simp_ensembles_:.cpp=.o))

