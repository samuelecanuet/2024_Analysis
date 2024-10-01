FASTERAC_CFLAGS = $(shell pkg-config  --cflags libfasterac)
FASTERAC_LIBS   = $(shell pkg-config  --libs   libfasterac)

ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs) -lSpectrum

GSL_CFLAGS      = $(shell gsl-config --cflags)
GSL_LIBS        = $(shell gsl-config --libs)

GTOOLS_L = $(shell gtools_lib)

CC        = g++
CFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS} ${GSL_CFLAGS} -Iinclude -I/softs/clhep/2.4.6.2/include/ -I/softs/clhep/2.4.6.2/include/CLHEP
LIBS      = ${FASTERAC_LIBS}   ${ROOT_LIBS} ${GTOOLS_L} ${GSL_LIBS} -lboost_filesystem -lboost_system include/Detectors.hh include/Messenger.hh -lCLHEP -L/softs/clhep/2.4.6.2/lib

SRCS      = $(wildcard src/*.c++)
OBJS      = $(patsubst src/%.c++,obj/%.o,$(SRCS))
EXE       = Grouper Matcher Merger Calibration Res_electronic LifeTime Defect ChannelTime

all : $(EXE)

Grouper: obj/Grouper.o

$(EXE): $(OBJS) include/Detectors.hh include/Messenger.hh
	${CC} -std=c++17 ${CFLAGS} ${LIBS} obj/$@.o -o $@ ../../lib/SignalDict/SignalDict.so -lfftw3 ../Simulation/Merger/AutoDict_vector_CLHEP__Hep3Vector__cxx.so ../Simulation/Merger/AutoDict_vector_vector_CLHEP__Hep3Vector___cxx.so ../Simulation/Merger/AutoDict_vector_vector_int____cxx.so ../Simulation/Merger/AutoDict_vector_vector_vector_CLHEP__Hep3Vector____cxx.so ../Simulation/Merger/AutoDict_vector_vector_vector_double______cxx.so

obj/%.o: src/%.c++ include/%.hh
	@mkdir -p obj
	${CC} -std=c++17 -c ${CFLAGS} $< -o $@

clean :
	rm -f obj/*.o
	rm -f $(EXE)


