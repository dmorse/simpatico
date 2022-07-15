mcMd_analyzers_mdSystem_=\
    mcMd/analyzers/mdSystem/MdAnalyzerFactory.cpp \
    mcMd/analyzers/mdSystem/MdEnergyAnalyzer.cpp \
    mcMd/analyzers/mdSystem/MdKineticEnergyAverage.cpp \
    mcMd/analyzers/mdSystem/MdPairEnergyCoefficients.cpp \
    mcMd/analyzers/mdSystem/MdPressureAnalyzer.cpp \
    mcMd/analyzers/mdSystem/MdStressAutoCorr.cpp 

ifdef SIMP_BOND
mcMd_analyzers_mdSystem_+=\
    mcMd/analyzers/mdSystem/MdIntraBondStressAutoCorr.cpp \
    mcMd/analyzers/mdSystem/MdIntraBondTensorAutoCorr.cpp 
endif

mcMd_analyzers_mdSystem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_mdSystem_))
mcMd_analyzers_mdSystem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_mdSystem_:.cpp=.o))

