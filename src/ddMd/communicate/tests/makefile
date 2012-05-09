SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include sources.mk

all: $(ddMd_communicate_tests_OBJS) BufferTest DomainTest PlanTest \
     AtomDistributorTest BondDistributorTest ExchangerTest ExchangerForceTest \
     AtomCollectorTest BondCollectorTest Test

clean:
	rm -f $(ddMd_communicate_tests_OBJS) BufferTest DomainTest \
        PlanTest AtomDistributorTest BondDistributorTest ExchangerTest \
        ExchangerForceTest AtomCollectorTest BondCollectorTest Test \
        $(ddMd_communicate_tests_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_communicate_tests_OBJS:.o=.d)

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

BondDistributorTest: BondDistributorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
                -o BondDistributorTest BondDistributorTest.o $(ddMd_LIB) $(util_LIB)

GroupDistributorTest: GroupDistributorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o GroupDistributorTest GroupDistributorTest.o $(ddMd_LIB) $(util_LIB)

ExchangerTest: ExchangerTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o ExchangerTest ExchangerTest.o $(ddMd_LIB) $(util_LIB)

ExchangerForceTest: ExchangerForceTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o ExchangerForceTest ExchangerForceTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)

AtomCollectorTest: AtomCollectorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o AtomCollectorTest AtomCollectorTest.o $(ddMd_LIB) $(util_LIB)

BondCollectorTest: BondCollectorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
        -o BondCollectorTest BondCollectorTest.o $(ddMd_LIB) $(util_LIB)

Test: Test.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
        -o Test Test.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include $(ddMd_communicate_tests_OBJS:.o=.d)

