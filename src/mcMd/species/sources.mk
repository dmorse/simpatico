
mcMd_species_SRCS=$(SRC_DIR)/mcMd/species/Diblock.cpp \
    $(SRC_DIR)/mcMd/species/Homopolymer.cpp \
    $(SRC_DIR)/mcMd/species/HomopolymerSG.cpp \
    $(SRC_DIR)/mcMd/species/HomoRing.cpp $(SRC_DIR)/mcMd/species/Linear.cpp \
    $(SRC_DIR)/mcMd/species/Point.cpp $(SRC_DIR)/mcMd/species/Ring.cpp \
    $(SRC_DIR)/mcMd/species/Species.cpp \
    $(SRC_DIR)/mcMd/species/species_queries.cpp \
    $(SRC_DIR)/mcMd/species/SpeciesFactory.cpp \
    $(SRC_DIR)/mcMd/species/SpeciesManager.cpp \
    $(SRC_DIR)/mcMd/species/SpeciesMutator.cpp 

mcMd_species_OBJS=$(mcMd_species_SRCS:.cpp=.o)

