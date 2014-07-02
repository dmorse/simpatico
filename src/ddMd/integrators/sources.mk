ddMd_integrators_=\
   ddMd/integrators/Integrator.cpp \
   ddMd/integrators/TwoStepIntegrator.cpp \
   ddMd/integrators/NveIntegrator.cpp \
   ddMd/integrators/NvtIntegrator.cpp \
   ddMd/integrators/NptIntegrator.cpp \
   ddMd/integrators/NphIntegrator.cpp \
   ddMd/integrators/IntegratorFactory.cpp

ddMd_integrators_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_integrators_))
ddMd_integrators_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_integrators_:.cpp=.o))

