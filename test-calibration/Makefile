# Makefile for calculating power spectrum from mock data
CC	= gcc
CFLAGS	= -Wall -O3

D_UTIL	= ./

LIB_ALL	= -lm

O_UTIL	= $(D_UTIL)/util.o $(D_UTIL)/fit_param.o 

ALLCODE = create-euclid-patch create-euclid-survey test-calibration 

all:
	@ for c in $(ALLCODE); do \
		(make $$c); \
	done

create-euclid-patch: create-euclid-patch.o 
	$(CC)  $(CFLAGS) -o $@ create-euclid-patch.o -lm

create-euclid-survey: create-euclid-survey.o 
	$(CC)  $(CFLAGS) -o $@ create-euclid-survey.o -lm

test-calibration: test-calibration.o $(O_UTIL) 
	$(CC)  $(CFLAGS) -o $@ test-calibration.o -lm $(O_UTIL) 

clean:
	-rm -f *.o create-euclid-patch create-euclid-survey test-calibration
