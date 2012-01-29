
SRC=../../src
CNV=.

include $(SRC)/makefile.global
include $(SRC)/makefile.inc

clean:	
	rm *.o

# ==========================================================================
# Test directory macros
CNV_PROG=$(CNV)/convert
 
$(CNV_PROG): $(CNV_PROG).cpp $(SRC_ALL) $(OBJ_ALL)
	$(CC) -o $(CNV_PROG) $(CNV_PROG).cpp $(OBJ_ALL)
