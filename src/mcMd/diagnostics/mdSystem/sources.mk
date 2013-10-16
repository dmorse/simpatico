mcMd_diagnostics_mdSystem_=\
    mcMd/diagnostics/mdSystem/MdDiagnosticFactory.cpp \
    mcMd/diagnostics/mdSystem/MdEnergyOutput.cpp \
    mcMd/diagnostics/mdSystem/MdKineticEnergyAverage.cpp \
    mcMd/diagnostics/mdSystem/MdPairEnergyCoefficients.cpp \
    mcMd/diagnostics/mdSystem/MdPotentialEnergyAverage.cpp \
    mcMd/diagnostics/mdSystem/MdPressureAverage.cpp \
    mcMd/diagnostics/mdSystem/MdIntraBondStressAutoCorr.cpp \
    mcMd/diagnostics/mdSystem/MdIntraBondTensorAutoCorr.cpp 

mcMd_diagnostics_mdSystem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_diagnostics_mdSystem_))
mcMd_diagnostics_mdSystem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_diagnostics_mdSystem_:.cpp=.o))

