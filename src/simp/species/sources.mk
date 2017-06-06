
simp_species_= \
    simp/species/Species.cpp \
    simp/species/Point.cpp 

ifdef SIMP_BOND
simp_species_+= \
    simp/species/Linear.cpp \
    simp/species/Homopolymer.cpp \
    simp/species/Diblock.cpp \
    simp/species/Multiblock.cpp \
    simp/species/Ring.cpp \
    simp/species/HomoRing.cpp 
endif

ifdef UTIL_MPI
simp_species_+= simp/species/SpeciesGroup.cpp
endif

simp_species_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_species_))
simp_species_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_species_:.cpp=.o))

