ddMd_diagnostics_=\
     ddMd/diagnostics/Diagnostic.cpp\
     ddMd/diagnostics/DiagnosticManager.cpp\
     ddMd/diagnostics/DiagnosticFactory.cpp\
     ddMd/diagnostics/WriteConfig.cpp\
     ddMd/diagnostics/OutputEnergy.cpp\
     ddMd/diagnostics/OutputPressure.cpp\
     ddMd/diagnostics/OutputBoxdim.cpp \
     ddMd/diagnostics/OutputTemperature.cpp\
     ddMd/diagnostics/OutputPairEnergies.cpp\
     ddMd/diagnostics/StructureFactor.cpp\
     ddMd/diagnostics/StructureFactorGrid.cpp\
     ddMd/diagnostics/OrderParamNucleation.cpp

ddMd_diagnostics_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_diagnostics_))
ddMd_diagnostics_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_diagnostics_:.cpp=.o))

