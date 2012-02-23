ddMd_communicate_tests_SRCS=$(SRC_DIR)/ddMd/communicate/tests/BufferTest.cc \
    $(SRC_DIR)/ddMd/communicate/tests/DomainTest.cc \
    $(SRC_DIR)/ddMd/communicate/tests/AtomDistributorTest.cc \
    $(SRC_DIR)/ddMd/communicate/tests/BondDistributorTest.cc \
    $(SRC_DIR)/ddMd/communicate/tests/PlanTest.cc \
    $(SRC_DIR)/ddMd/communicate/tests/ExchangerTest.cc \
    $(SRC_DIR)/ddMd/communicate/tests/CollectorTest.cc \
    $(SRC_DIR)/ddMd/communicate/tests/Test.cc 

ddMd_communicate_tests_OBJS=$(ddMd_communicate_tests_SRCS:.cc=.o)

