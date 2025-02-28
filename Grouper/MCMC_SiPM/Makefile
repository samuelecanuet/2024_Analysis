FASTERAC_CFLAGS = $(shell pkg-config  --cflags libfasterac)
FASTERAC_LIBS   = $(shell pkg-config  --libs   libfasterac)

ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs) -lSpectrum

GSL_CFLAGS      = $(shell gsl-config --cflags)
GSL_LIBS        = $(shell gsl-config --libs)

CC        = g++
CFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS} ${GSL_CFLAGS} -Iinclude
LIBS      = ${FASTERAC_LIBS}   ${ROOT_LIBS} ${GSL_LIBS} -lboost_filesystem -lboost_system -lrt ../include/Detectors.hh ../include/Messenger.hh ../include/Utilities.hh ../include/SIPM_Calibration.hh

SRCS      = $(wildcard *.c++)
OBJS      = $(patsubst %.c++,%.o,$(SRCS))
EXE       = MCMC_SiPM

all : $(EXE)

$(EXE): $(OBJS)
	${CC} -std=c++17 ${CFLAGS} ${LIBS} $@.o -o $@ ../../../lib/SignalDict/SignalDict.so -lfftw3 ../Dic/AutoDict_vector_vector_std__pair_Signal__Signal____cxx.so

%.o: %.c++ %.hh
	@mkdir -p obj
	${CC} -std=c++17 -c ${CFLAGS} $< -o $@

clean :
	rm -f *.o
	rm -f $(EXE)


