FASTERAC_CFLAGS = $(shell pkg-config  --cflags libfasterac)
FASTERAC_LIBS   = $(shell pkg-config  --libs   libfasterac)

ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs) -lSpectrum

GSL_CFLAGS      = $(shell gsl-config --cflags)
GSL_LIBS        = $(shell gsl-config --libs)

GTOOLS_L = $(shell gtools_lib)

CC        = g++
CFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS} ${GSL_CFLAGS} -Iinclude
LIBS      = ${FASTERAC_LIBS}   ${ROOT_LIBS} ${GTOOLS_L} ${GSL_LIBS} -lboost_filesystem -lboost_system include/Detectors.hh include/Messenger.hh 

SRCS      = $(wildcard src/*.c++)
OBJS      = $(patsubst src/%.c++,obj/%.o,$(SRCS))
EXE       = Grouper Matcher Merger Calibration Res_electronic LifeTime Defect SIPM_Calibration Bi207 SIPM_Matching SIPM_Dark

all : $(EXE)

Grouper: obj/Grouper.o
Calibration: obj/Calibration.o
Matcher: obj/Matcher.o
Merger: obj/Merger.o
Res_electronic: obj/Res_electronic.o
LifeTime: obj/LifeTime.o
Defect: obj/Defect.o
SIPM_Calibration: obj/SIPM_Calibration.o
SIPM_Matching: obj/SIPM_Matching.o
SIPM_Dark: obj/SIPM_Dark.o	
Bi207: obj/Bi207.o

$(EXE): $(OBJS) include/Detectors.hh include/Messenger.hh
	${CC} -std=c++17 ${CFLAGS} ${LIBS} obj/$@.o -o $@ ../../lib/SignalDict/SignalDict.so -lfftw3 

obj/%.o: src/%.c++ include/%.hh
	@mkdir -p obj
	${CC} -std=c++17 -c ${CFLAGS} $< -o $@

clean :
	rm -f obj/*.o
	rm -f $(EXE)


