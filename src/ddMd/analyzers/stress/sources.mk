ddMd_analyzers_stress_=\
     ddMd/analyzers/stress/PressureAnalyzer.cpp\
     ddMd/analyzers/stress/OutputPressure.cpp\
     ddMd/analyzers/stress/OutputBoxdim.cpp \
     ddMd/analyzers/stress/StressAnalyzer.cpp\
     ddMd/analyzers/stress/VirialStressAnalyzer.cpp\
     ddMd/analyzers/stress/OutputStressTensor.cpp\
     ddMd/analyzers/stress/VirialStressTensorAverage.cpp\
     ddMd/analyzers/stress/VirialStressTensor.cpp\
     ddMd/analyzers/stress/StressAutoCorr.cpp\
     ddMd/analyzers/stress/StressAutoCorrelation.cpp

ddMd_analyzers_stress_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_analyzers_stress_))
ddMd_analyzers_stress_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_analyzers_stress_:.cpp=.o))

