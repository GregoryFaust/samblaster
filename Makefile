# Determine the samblaster build number
BUILDNUM = 25
# INTERNAL = TRUE

OBJS = samblaster.o sbhash.o

ifdef INTERNAL
PROG = samblaster$(BUILDNUM)
CCFLAGS  = -Wall -Winline -O0 -g -D BUILDNUM=$(BUILDNUM)
else
PROG = samblaster
CCFLAGS  = -Wall -Werror=literal-suffix -O3 -D BUILDNUM=$(BUILDNUM)
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
