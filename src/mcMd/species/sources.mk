mcMd_species_= \
    mcMd/species/SpeciesMutator.cpp \
    mcMd/species/SpeciesFactory.cpp \
    mcMd/species/SpeciesManager.cpp 

ifdef SIMP_BOND
mcMd_species_+= \
    mcMd/species/HomopolymerSG.cpp \
    mcMd/species/LinearSG.cpp
endif

mcMd_species_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_species_))
mcMd_species_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_species_:.cpp=.o))

