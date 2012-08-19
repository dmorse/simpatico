ddMd_integrators_SRCS=$(SRC_DIR)/ddMd/integrators/Integrator.cpp \
                     $(SRC_DIR)/ddMd/integrators/TwoStepIntegrator.cpp \
                     $(SRC_DIR)/ddMd/integrators/NveIntegrator.cpp \
                     $(SRC_DIR)/ddMd/integrators/NvtIntegrator.cpp \
                     $(SRC_DIR)/ddMd/integrators/IntegratorFactory.cpp

ddMd_integrators_OBJS=$(ddMd_integrators_SRCS:.cpp=.o)
