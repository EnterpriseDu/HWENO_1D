CC = gcc
CmpFLAGS = -g -c
LnkFLAGS = ../file_io-new/file_io.a  ../Riemann_solver/Riemann_solver.a  -lm

INCS = $(patsubst %, -I%, $(filter-out ../../.git%, $(shell find -type d)))
INCS +=  -I../file_io-new/inc/  -I../Riemann_solver/inc/
VPATH += $(patsubst  -I%,%,$(INCS))
export VPATH


#find all the directories in './src', but drop './src' itself
#be careful above the order of the list, solvers should go before others
#here since solvers begin with capital letters, we can achieve this by
#alphabetically sort the list
LIST = $(sort $(filter-out %src,  $(shell find $(CURDIR)/src -type d)  )  )
DIRS = $(notdir   $(foreach dir, $(notdir $(LIST)), $(wildcard ./src/$(dir)))  )
#         ^                                               ^
#         |                                               |
#         |      drop all the sub-directories, but this will add the prefix './src/'
#         |
#         |
#  so we drop the prefix again
#construct the list of archives
ARXIVS = $(foreach dir, $(DIRS), ./src/$(dir)/$(dir).a)




all :
include main.d
include riemann_sol.d

%.d : %.c
	@set -e; rm -f $@; \
	$(CC) -MM $(INCS)  $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$



all : main.o $(ARXIVS) -lm
	gcc -o riemann_sol $< $(ARXIVS) $(LnkFLAGS)
	gcc -o main $< $(ARXIVS) $(LnkFLAGS)
all : riemann_sol.o $(ARXIVS) -lm


main.o :
	gcc $< $(CmpFLAGS) $(patsubst %, -I%,$(sort $(dir $(filter %.h, $^))))
riemann_sol.o :
	gcc $< $(CmpFLAGS) $(patsubst %, -I%,$(sort $(dir $(filter %.h, $^))))

$(ARXIVS) :
	$(MAKE) --directory=$(@D)




clean :
	@for d in $(DIRS); do \
		$(MAKE) --directory=src/$$d clean; \
	done
	rm -f *.[oa]
	rm -f main


test :
	@echo ===================

.PHONY : clean test all simple $(ARXIVS)
