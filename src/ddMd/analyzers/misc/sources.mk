ddMd_analyzers_misc_=\
     ddMd/analyzers/misc/OrderParamNucleation.cpp

ifdef SIMP_BOND
ddMd_analyzers_misc_+=\
     ddMd/analyzers/misc/BondTensorAutoCorr.cpp
endif

ddMd_analyzers_misc_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_analyzers_misc_))
ddMd_analyzers_misc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_analyzers_misc_:.cpp=.o))

