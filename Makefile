# Determine the samblaster build number
BUILDNUM = 24
# INTERNAL = TRUE

OBJS = samblaster.o sbhash.o

ifdef INTERNAL
PROG = samblaster$(BUILDNUM)
CCFLAGS  = -Wall -Winline -O0 -g -D BUILDNUM=$(BUILDNUM)
else
PROG = samblaster
CCFLAGS  = -Wall -O3 -D BUILDNUM=$(BUILDNUM)
endif

CC	 = gcc
CPP      = g++

CFLAGS = $(CCFLAGS) -std=gnu99
CPPFLAGS = $(CCFLAGS)

SAMBLASTER: $(PROG)

$(PROG): $(OBJS)
	$(CPP) $(LDFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

clean:
	rm -f $(PROG)  *.o *.co *.linkinfo
