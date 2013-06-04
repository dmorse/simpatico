SRC_DIR_REL =../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/inter/defines.mk
include $(SRC_DIR)/ddMd/defines.mk
include $(TESTS_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk


include sources.mk

all: $(tests_ddMd_communicate_OBJS) BufferTest DomainTest PlanTest \
        AtomDistributorTest GroupDistributorTest \
        ExchangerTest ExchangerForceTest \
        AtomCollectorTest BondCollectorTest Test

clean:
	rm -f $(tests_ddMd_communicate_OBJS) \
              BufferTest DomainTest PlanTest \
              AtomDistributorTest GroupDistributorTest \
              ExchangerTest ExchangerForceTest \
              AtomCollectorTest BondCollectorTest Test \
              $(tests_ddMd_communicate_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_ddMd_communicate_OBJS:.o=.d)

BufferTest: BufferTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o BufferTest BufferTest.o $(ddMd_LIB) $(util_LIB)

DomainTest: DomainTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o DomainTest DomainTest.o $(ddMd_LIB) $(util_LIB)

PlanTest: PlanTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o PlanTest PlanTest.o $(ddMd_LIB) $(util_LIB)

AtomDistributorTest: AtomDistributorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o AtomDistributorTest AtomDistributorTest.o $(ddMd_LIB) $(util_LIB)

GroupDistributorTest: GroupDistributorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o GroupDistributorTest GroupDistributorTest.o $(ddMd_LIB) $(util_LIB)

ExchangerTest: ExchangerTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o ExchangerTest ExchangerTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)

ExchangerForceTest: ExchangerForceTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o ExchangerForceTest ExchangerForceTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)

AtomCollectorTest: AtomCollectorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o AtomCollectorTest AtomCollectorTest.o $(ddMd_LIB) $(util_LIB)

BondCollectorTest: BondCollectorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
        -o BondCollectorTest BondCollectorTest.o $(ddMd_LIB) $(util_LIB)

Test: Test.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
        -o Test Test.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include $(tests_ddMd_communicate_OBJS:.o=.d)

