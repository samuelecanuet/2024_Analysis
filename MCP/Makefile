FASTERAC_CFLAGS = $(shell pkg-config  --cflags libfasterac)
FASTERAC_LIBS   = $(shell pkg-config  --libs   libfasterac)

ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs) -lSpectrum

GSL_CFLAGS      = $(shell gsl-config --cflags)
GSL_LIBS        = $(shell gsl-config --libs)

algib_PATH = ./home/local1/Downloads/alglib-4.05.0.cpp.gpl/alglib-cpp/
ALGLIB_SRCS = src/ap.cpp src/interpolation.cpp src/linalg.cpp

CC        = g++
CFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS} ${GSL_CFLAGS} -Iinclude
LIBS      = ${FASTERAC_LIBS}   ${ROOT_LIBS} ${GSL_LIBS} -lboost_filesystem -lboost_system ../Grouper/include/Detectors.hh ../Grouper/include/Messenger.hh ../Grouper/include/Utilities.hh

SRCS      = $(wildcard src/*.c++) $(algib_PATH + $(ALGLIB_SRCS))
OBJS      = $(patsubst src/%.c++,obj/%.o,$(SRCS))
EXE       = Calibration Grouper Grouper_fit

all : $(EXE)

$(EXE): $(OBJS)
	${CC} -std=c++17 obj/$@.o -o $@ ../Grouper/../../lib/SignalDict/SignalDict.so ../Grouper/Dic/AutoDict_vector_vector_std__pair_Signal__Signal____cxx.so ${CFLAGS} ${LIBS}

obj/%.o: src/%.c++ include/%.hh
	@mkdir -p obj
	${CC} -std=c++17 -c $< -o $@ ${CFLAGS}

clean :
	rm -f obj/*.o
	rm -f $(EXE)


