
mcMd_species_= \
    mcMd/species/Species.cpp \
    mcMd/species/SpeciesMutator.cpp \
    mcMd/species/SpeciesFactory.cpp \
    mcMd/species/SpeciesManager.cpp \
    mcMd/species/Point.cpp \
    mcMd/species/species_queries.cpp 

ifdef INTER_BOND
mcMd_species_+= \
    mcMd/species/Linear.cpp \
    mcMd/species/Homopolymer.cpp \
    mcMd/species/Diblock.cpp \
    mcMd/species/Multiblock.cpp \
    mcMd/species/HomopolymerSG.cpp \
    mcMd/species/Ring.cpp \
    mcMd/species/HomoRing.cpp 
endif

mcMd_species_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_species_))
mcMd_species_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_species_:.cpp=.o))

