ROOT=../bcl

CPPFLAGS = -I$(ROOT)/include -I.
CFLAGS   = -Wall -Wextra -pedantic -std=c99 -Werror 
LDFLAGS  = -L$(ROOT)/lib
LDLIBS   = -lbcl -lm -lfftw3

BIN=colorization

.PHONY: all
all: $(BIN)

.PHONY: clean cleanall
clean:
	$(RM) *.o *.ppm
cleanall: clean
	$(RM) $(BIN)
