# for C++ define  CC = g++
CC = g++
CFLAGS = -g -Wall -Weffc++ -fPIC -m64 -std=gnu++11
LFLAGS = -L$(ROOTSYS)/lib 
#LFLAGS = -L$(ROOTSYS)/lib -L$(ALICE_ROOT)/lib
INC =	-I$(ROOTSYS)/include -I./
#INC =	-I$(ROOTSYS)/include  -I$(ALICE_ROOT)/include -I./
TGT =	libTrackO2.so
DICT=	trackO2Dict.cxx
DICTO=	trackO2Dict.o

SRC = 	Track.cxx TrackIO.cxx

HDR =	$(SRC:.cxx=.h) Utils.h Constants.h

OBJ = 	$(SRC:.cxx=.o)


.PHONY: depend clean

all: 	$(TGT)
	@echo creating libTrackO2.so

$(TGT):	$(OBJ) $(DICTO)
	$(CC) $(CFLAGS)  -shared -o $(TGT) $(OBJ) $(DICTO) `root-config --ldflags` $(LFLAGS)

# pull in dependency info for *existing* .o files
-include $(OBJ:.o=.d)

%.o : %.cxx
	$(CC) $(CFLAGS) $(INC) -c $<  -o $@
	$(CC) -MM $(CFLAGS) $(INC) -c $*.cxx > $*.d
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean:
	rm -f *.o *~ *.so *.d *Dict.{h,cxx}

$(DICT): $(HDR) trackO2LinkDef.h
	rootcint -f $@ -c $(INC) $(HDR) $^


depend: $(SRC)
	makedepend $(INC) $^

# DO NOT DELETE THIS LINE -- make depend needs it
