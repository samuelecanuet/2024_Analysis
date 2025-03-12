FASTERAC_CFLAGS = $(shell pkg-config  --cflags libfasterac)
FASTERAC_LIBS   = $(shell pkg-config  --libs   libfasterac)

ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs) -lSpectrum

GSL_CFLAGS      = $(shell gsl-config --cflags)
GSL_LIBS        = $(shell gsl-config --libs)


CC        = g++
CFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS} ${GSL_CFLAGS} -Iinclude  -I/softs/clhep/2.4.6.2/include/ -I/softs/clhep/2.4.6.2/include/CLHEP
LIBS      = ${FASTERAC_LIBS}   ${ROOT_LIBS} ${GSL_LIBS} -lboost_filesystem -lboost_system ../Grouper/include/Messenger.hh -lCLHEP -L/softs/clhep/2.4.6.2/lib

SRCS      = $(wildcard src/*.c++)
OBJS      = $(patsubst src/%.c++,obj/%.o,$(SRCS))
EXE       = Reader_CRADLE RMATRIX Reader2 ScatteringModel ReaderNew InterStrip

all : $(EXE)

$(EXE): $(OBJS) ../Grouper/include/Messenger.hh
	${CC} -std=c++17 ${CFLAGS} ${LIBS} obj/$@.o -o $@ ../../lib/SignalDict/SignalDict.so Merger/AutoDict_vector_CLHEP__Hep3Vector__cxx.so Merger/AutoDict_vector_vector_CLHEP__Hep3Vector___cxx.so Merger/AutoDict_vector_vector_int____cxx.so Merger/AutoDict_vector_vector_vector_CLHEP__Hep3Vector____cxx.so Merger/AutoDict_vector_vector_vector_double______cxx.so

obj/%.o: src/%.c++ include/%.hh
	@mkdir -p obj
	${CC} -std=c++17 -c ${CFLAGS} $< -o $@

clean :
	rm -f obj/*.o
	rm -f $(EXE)


