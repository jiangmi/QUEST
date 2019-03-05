QUEST_DIR = ../..
include $(QUEST_DIR)/make.inc

FLAGS=$(FC_FLAGS) -I$(QUEST_DIR)/src

.SUFFIXES : .mod .o .F90

test:   $(QUEST_DIR)/$(DQMCLIB) main.F90
	$(FC) $(FLAGS) -o main.e main.F90 $(QUEST_DIR)/$(DQMCLIB) $(QUEST_DIR)/$(SPRNGLIB)  $(LIB)

clean:
	rm -f main.e
