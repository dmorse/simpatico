include $(SRC_DIR)/ddMd/analyzers/trajectory/sources.mk
include $(SRC_DIR)/ddMd/analyzers/energy/sources.mk
include $(SRC_DIR)/ddMd/analyzers/stress/sources.mk
include $(SRC_DIR)/ddMd/analyzers/scattering/sources.mk
include $(SRC_DIR)/ddMd/analyzers/misc/sources.mk

ddMd_analyzers_=\
 $(ddMd_analyzers_trajectory_)\
 $(ddMd_analyzers_energy_)\
 $(ddMd_analyzers_stress_)\
 $(ddMd_analyzers_scattering_)\
 $(ddMd_analyzers_misc_)

ddMd_analyzers_+=\
     ddMd/analyzers/Analyzer.cpp\
     ddMd/analyzers/AnalyzerManager.cpp\
     ddMd/analyzers/AnalyzerFactory.cpp\
     ddMd/analyzers/AverageAnalyzer.cpp\
     ddMd/analyzers/TensorAverageAnalyzer.cpp\
     ddMd/analyzers/SymmTensorAverageAnalyzer.cpp

ddMd_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_analyzers_))
ddMd_analyzers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_analyzers_:.cpp=.o))

