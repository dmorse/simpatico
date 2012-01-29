SRC=../../../src

include $(SRC)/makefile.cc
include $(SRC)/mcMd/makefile.inc
RINGMAKER=RingMaker

clean:	
	-rm $(RINGMAKER)

# ==========================================================================

$(RINGMAKER): $(RINGMAKER).cpp $(OBJ_MCMD_ALL)
	$(CC) -o $(RINGMAKER) $(RINGMAKER).cpp $(OBJ_MCMD_ALL)
