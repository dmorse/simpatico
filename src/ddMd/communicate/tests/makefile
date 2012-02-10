SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include sources.mk

all: $(ddMd_communicate_tests_OBJS) BufferTest DomainTest PlanTest \
     AtomDistributorTest BondDistributorTest ExchangerTest CollectorTest Test

clean:
	rm -f $(ddMd_communicate_tests_OBJS) BufferTest DomainTest ExchangerTest Test CollectorTest AtomDistributorTest $(ddMd_communicate_tests_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_communicate_tests_OBJS:.o=.d)

BufferTest: BufferTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o BufferTest BufferTest.o \
	       $(ddMd_LIB) $(util_LIB)

DomainTest: DomainTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o DomainTest DomainTest.o \
	       $(ddMd_LIB) $(util_LIB)

PlanTest: PlanTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o PlanTest PlanTest.o \
	       $(ddMd_LIB) $(util_LIB)

AtomDistributorTest: AtomDistributorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AtomDistributorTest AtomDistributorTest.o \
	       $(ddMd_LIB) $(util_LIB)

BondDistributorTest: BondDistributorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o BondDistributorTest BondDistributorTest.o \
	       $(ddMd_LIB) $(util_LIB)

ExchangerTest: ExchangerTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o ExchangerTest ExchangerTest.o \
	       $(ddMd_LIB) $(util_LIB)

CollectorTest: CollectorTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o CollectorTest CollectorTest.o \
	       $(ddMd_LIB) $(util_LIB)

Test: Test.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(ddMd_LIB) $(util_LIB)

-include $(ddMd_communicate_tests_OBJS:.o=.d)

