# Makefile for calculating power spectrum from mock data
CC	= gcc
CFLAGS	= -Wall -O3 -std=c99

BIN = ../bin
D_UTIL	= $(BIN)

LIB_ALL	= -lm

O_UTIL	= $(D_UTIL)/util.o $(D_UTIL)/fit_param.o 

ALLCODE = create-euclid-patch test-calibration 

all:
	@ for c in $(ALLCODE); do \
		(make $$c); \
	done

$(BIN)/%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

create-euclid-patch: $(BIN)/create-euclid-patch.o
	$(CC) $(CFLAGS) -o $(BIN)/$@ $(BIN)/create-euclid-patch.o -lm

test-calibration: $(BIN)/test-calibration.o $(O_UTIL)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $(BIN)/test-calibration.o -lm $(O_UTIL) 

clean:
	-rm -f $(BIN)/*.o $(BIN)/create-euclid-patch $(BIN)/test-calibration