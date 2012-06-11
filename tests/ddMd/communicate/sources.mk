tests_ddMd_communicate_SRCS=$(TESTS_DIR)/ddMd/communicate/BufferTest.cc \
    $(TESTS_DIR)/ddMd/communicate/DomainTest.cc \
    $(TESTS_DIR)/ddMd/communicate/AtomDistributorTest.cc \
    $(TESTS_DIR)/ddMd/communicate/BondDistributorTest.cc \
    $(TESTS_DIR)/ddMd/communicate/GroupDistributorTest.cc \
    $(TESTS_DIR)/ddMd/communicate/AtomCollectorTest.cc \
    $(TESTS_DIR)/ddMd/communicate/BondCollectorTest.cc \
    $(TESTS_DIR)/ddMd/communicate/PlanTest.cc \
    $(TESTS_DIR)/ddMd/communicate/ExchangerTest.cc \
    $(TESTS_DIR)/ddMd/communicate/ExchangerForceTest.cc \
    $(TESTS_DIR)/ddMd/communicate/Test.cc 

tests_ddMd_communicate_OBJS=$(tests_ddMd_communicate_SRCS:.cc=.o)

