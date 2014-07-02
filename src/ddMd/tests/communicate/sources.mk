ddMd_tests_communicate_= \
    ddMd/tests/communicate/BufferTest.cc \
    ddMd/tests/communicate/DomainTest.cc \
    ddMd/tests/communicate/AtomDistributorTest.cc \
    ddMd/tests/communicate/GroupDistributorTest.cc \
    ddMd/tests/communicate/AtomCollectorTest.cc \
    ddMd/tests/communicate/BondCollectorTest.cc \
    ddMd/tests/communicate/PlanTest.cc \
    ddMd/tests/communicate/ExchangerTest.cc \
    ddMd/tests/communicate/ExchangerForceTest.cc \
    ddMd/tests/communicate/Test.cc 

ddMd_tests_communicate_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_communicate_))
ddMd_tests_communicate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_communicate_:.cc=.o))

