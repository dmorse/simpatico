ddMd_tests_communicate_= \
    ddMd/tests/communicate/BufferTest.cpp \
    ddMd/tests/communicate/DomainTest.cpp \
    ddMd/tests/communicate/AtomDistributorTest.cpp \
    ddMd/tests/communicate/GroupDistributorTest.cpp \
    ddMd/tests/communicate/AtomCollectorTest.cpp \
    ddMd/tests/communicate/BondCollectorTest.cpp \
    ddMd/tests/communicate/PlanTest.cpp \
    ddMd/tests/communicate/ExchangerTest.cpp \
    ddMd/tests/communicate/ExchangerForceTest.cpp \
    ddMd/tests/communicate/Test.cpp 

ddMd_tests_communicate_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_communicate_))
ddMd_tests_communicate_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_tests_communicate_:.cpp=.o))

