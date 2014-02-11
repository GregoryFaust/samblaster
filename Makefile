# Determine the samblaster build number
BUILDNUM = 9
# USERMODE = TRUE

PROG = samblaster
OBJS = samblaster.o

ifdef USERMODE
CCFLAGS  = -Wall -Winline -O3 -D COMPILE_USER_MODE -D BUILDNUM=$(BUILDNUM)
else
CCFLAGS  = -Wall -Winline -O3 -g -D BUILDNUM=$(BUILDNUM)
endif

CC	 = gcc
CPP      = g++

# We need gnu++0x to get the hash table set in stl.
CFLAGS = $(CCFLAGS) -std=gnu99
CPPFLAGS = $(CCFLAGS) -std=gnu++0x

SAMBLASTER: $(PROG)

$(PROG): $(OBJS)
	$(CPP) $(LDFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

clean:
	rm -f $(PROG)  *.o *.co *.linkinfo
