
mcMd_species_=mcMd/species/Diblock.cpp \
    mcMd/species/Homopolymer.cpp \
    mcMd/species/HomopolymerSG.cpp \
    mcMd/species/HomoRing.cpp mcMd/species/Linear.cpp \
    mcMd/species/Point.cpp mcMd/species/Ring.cpp \
    mcMd/species/Species.cpp \
    mcMd/species/species_queries.cpp \
    mcMd/species/SpeciesFactory.cpp \
    mcMd/species/SpeciesManager.cpp \
    mcMd/species/SpeciesMutator.cpp 

mcMd_species_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_species_))
mcMd_species_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_species_:.cpp=.o))

