mcMd_colvars_= \
    mcMd/colvars/ColVar.cpp \
    mcMd/colvars/MdColVar.cpp \
    mcMd/colvars/CvPotential.cpp \
    mcMd/colvars/CvHarmonicPotential.cpp 

mcMd_colvars_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_colvars_))
mcMd_colvars_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_colvars_:.cpp=.o))

