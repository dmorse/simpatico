mcMd_analyzers_mcSystem_=\
    mcMd/analyzers/mcSystem/McBondEnergyAverage.cpp \
    mcMd/analyzers/mcSystem/McAnalyzerFactory.cpp \
    mcMd/analyzers/mcSystem/McEnergyAverage.cpp \
    mcMd/analyzers/mcSystem/McEnergyOutput.cpp \
    mcMd/analyzers/mcSystem/McPairEnergyAverage.cpp \
    mcMd/analyzers/mcSystem/McPressureAverage.cpp \
    mcMd/analyzers/mcSystem/McVirialStressTensorAverage.cpp \
    mcMd/analyzers/mcSystem/McStressAutoCorr.cpp \
    mcMd/analyzers/mcSystem/McIntraBondStressAutoCorr.cpp \
    mcMd/analyzers/mcSystem/McIntraBondTensorAutoCorr.cpp \
    mcMd/analyzers/mcSystem/McNVTChemicalPotential.cpp 

ifdef INTER_EXTERNAL
mcMd_analyzers_mcSystem_+=\
    mcMd/analyzers/mcSystem/McExternalEnergyAverage.cpp 
endif

mcMd_analyzers_mcSystem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_mcSystem_))
mcMd_analyzers_mcSystem_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_analyzers_mcSystem_:.cpp=.o))

