FASTERAC_CFLAGS = $(shell pkg-config  --cflags libfasterac)
FASTERAC_LIBS   = $(shell pkg-config  --libs   libfasterac)

ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs) -lSpectrum

GSL_CFLAGS      = $(shell gsl-config --cflags)
GSL_LIBS        = $(shell gsl-config --libs)

CC        = g++
CFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS} ${GSL_CFLAGS} -Iinclude
LIBS      = ${FASTERAC_LIBS}   ${ROOT_LIBS}   ${GSL_LIBS} -lboost_filesystem -lboost_system -lfftw3

SRCS      = $(wildcard src/*.c++)
OBJS      = $(patsubst src/%.c++,obj/%.o,$(SRCS))
EXE       = FullFunction

all : $(EXE)

$(EXE): $(OBJS)
	${CC} -std=c++17 ${CFLAGS} ${OBJS} ${LIBS} -o $@ ../../lib/SignalDict/SignalDict.so

obj/%.o: src/%.c++ include/%.hh
	@mkdir -p obj
	${CC} -std=c++17 -c ${CFLAGS} $< -o $@

clean :
	rm -f obj/*.o
	rm -f $(EXE)