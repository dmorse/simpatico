mcMd_analyzers_mcSystem_=\
    mcMd/analyzers/mcSystem/McEnergyAnalyzer.cpp \
    mcMd/analyzers/mcSystem/McEnergyAverage.cpp \
    mcMd/analyzers/mcSystem/McEnergyOutput.cpp \
    mcMd/analyzers/mcSystem/McPairEnergyAverage.cpp \
    mcMd/analyzers/mcSystem/McPressureAverage.cpp \
    mcMd/analyzers/mcSystem/McVirialStressTensorAverage.cpp \
    mcMd/analyzers/mcSystem/McAnalyzerFactory.cpp 

ifndef SIMP_NOPAIR
mcMd_analyzers_mcSystem_+=\
    mcMd/analyzers/mcSystem/McMuExchange.cpp 
endif

ifdef SIMP_BOND
mcMd_analyzers_mcSystem_+=\
    mcMd/analyzers/mcSystem/McBondEnergyAverage.cpp \
    mcMd/analyzers/mcSystem/McIntraBondStressAutoCorr.cpp \
    mcMd/analyzers/mcSystem/McIntraBondTensorAutoCorr.cpp \
    mcMd/analyzers/mcSystem/McNVTChemicalPotential.cpp 
endif

ifdef SIMP_EXTERNAL
mcMd_analyzers_mcSystem_+=\
    mcMd/analyzers/mcSystem/McExternalEnergyAverage.cpp 
endif

mcMd_analyzers_mcSystem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_mcSystem_))
mcMd_analyzers_mcSystem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_mcSystem_:.cpp=.o))

