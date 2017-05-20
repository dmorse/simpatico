ddMd_analyzers_energy_=\
     ddMd/analyzers/energy/LogEnergy.cpp\
     ddMd/analyzers/energy/OutputEnergy.cpp\
     ddMd/analyzers/energy/KineticEnergyAnalyzer.cpp\
     ddMd/analyzers/energy/EnergyAnalyzer.cpp\
     ddMd/analyzers/energy/OutputTemperature.cpp\
     ddMd/analyzers/energy/PairEnergyAverage.cpp\
     ddMd/analyzers/energy/PairEnergyAnalyzer.cpp\
     ddMd/analyzers/energy/OutputPairEnergies.cpp

ifdef SIMP_EXTERNAL
ddMd_analyzers_energy_+=\
     ddMd/analyzers/energy/ExternalEnergyAverage.cpp\
     ddMd/analyzers/energy/ExternalEnergyAnalyzer.cpp
endif

ddMd_analyzers_energy_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_analyzers_energy_))
ddMd_analyzers_energy_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_analyzers_energy_:.cpp=.o))

