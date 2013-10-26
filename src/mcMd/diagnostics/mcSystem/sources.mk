mcMd_diagnostics_mcSystem_=\
    mcMd/diagnostics/mcSystem/McBondEnergyAverage.cpp \
    mcMd/diagnostics/mcSystem/McDiagnosticFactory.cpp \
    mcMd/diagnostics/mcSystem/McEnergyAverage.cpp \
    mcMd/diagnostics/mcSystem/McEnergyOutput.cpp \
    mcMd/diagnostics/mcSystem/McPairEnergyAverage.cpp \
    mcMd/diagnostics/mcSystem/McPressureAverage.cpp \
    mcMd/diagnostics/mcSystem/McStressAutoCorr.cpp \
    mcMd/diagnostics/mcSystem/McIntraBondStressAutoCorr.cpp \
    mcMd/diagnostics/mcSystem/McIntraBondTensorAutoCorr.cpp \
    mcMd/diagnostics/mcSystem/McNVTChemicalPotential.cpp 

ifdef INTER_EXTERNAL
mcMd_diagnostics_mcSystem_+=\
    mcMd/diagnostics/mcSystem/McExternalEnergyAverage.cpp 
endif

mcMd_diagnostics_mcSystem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_diagnostics_mcSystem_))
mcMd_diagnostics_mcSystem_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_diagnostics_mcSystem_:.cpp=.o))

