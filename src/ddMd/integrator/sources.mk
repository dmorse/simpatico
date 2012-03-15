ddMd_integrator_SRCS=$(SRC_DIR)/ddMd/integrator/Integrator.cpp \
                     $(SRC_DIR)/ddMd/integrator/IntegratorFactory.cpp \
                     $(SRC_DIR)/ddMd/integrator/NveIntegrator.cpp \
                     $(SRC_DIR)/ddMd/integrator/NvtIntegrator.cpp 

ddMd_integrator_OBJS=$(ddMd_integrator_SRCS:.cpp=.o)
