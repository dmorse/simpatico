mcMd_analyzers_mdSystem_=\
    mcMd/analyzers/mdSystem/MdAnalyzerFactory.cpp \
    mcMd/analyzers/mdSystem/MdEnergyOutput.cpp \
    mcMd/analyzers/mdSystem/MdKineticEnergyAverage.cpp \
    mcMd/analyzers/mdSystem/MdPairEnergyCoefficients.cpp \
    mcMd/analyzers/mdSystem/MdPotentialEnergyAverage.cpp \
    mcMd/analyzers/mdSystem/MdPressureAverage.cpp \
    mcMd/analyzers/mdSystem/MdStressAutoCorrelation.cpp \
    mcMd/analyzers/mdSystem/MdVirialStressTensorAverage.cpp \
    mcMd/analyzers/mdSystem/MdIntraBondStressAutoCorr.cpp \
    mcMd/analyzers/mdSystem/MdIntraBondTensorAutoCorr.cpp 

mcMd_analyzers_mdSystem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_mdSystem_))
mcMd_analyzers_mdSystem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_mdSystem_:.cpp=.o))

