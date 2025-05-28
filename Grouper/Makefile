FASTERAC_CFLAGS = $(shell pkg-config  --cflags libfasterac)
FASTERAC_LIBS   = $(shell pkg-config  --libs   libfasterac)

ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs) -lSpectrum

GSL_CFLAGS      = $(shell gsl-config --cflags)
GSL_LIBS        = $(shell gsl-config --libs)

CC        = g++
CFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS} ${GSL_CFLAGS} -Iinclude
LIBS      = ${FASTERAC_LIBS}   ${ROOT_LIBS} ${GSL_LIBS} -lboost_filesystem -lboost_system include/Detectors.hh include/Messenger.hh include/Utilities.hh

SRCS      = $(wildcard src/*.c++)
OBJS      = $(patsubst src/%.c++,obj/%.o,$(SRCS))
EXE       = Grouper Matcher Merger Calibration Res_electronic LifeTime Defect SIPM_Calibration Bi207 SIPM_Matching SIPM_Dark SiPM_Time Analyser Source_Grouper BasicAnalyser Rate Pulse SIPM_PileUp Time_Pulse_Calibration

all : $(EXE)

Grouper: obj/Grouper.o
Calibration: obj/Calibration.o
Matcher: obj/Matcher.o
Merger: obj/Merger.o
Res_electronic: obj/Res_electronic.o
LifeTime: obj/LifeTime.o
Defect: obj/Defect.o
Source_Grouper: obj/Source_Grouper.o
SIPM_Calibration: obj/SIPM_Calibration.o
SIPM_Matching: obj/SIPM_Matching.o
SIPM_Dark: obj/SIPM_Dark.o	
SiPM_Time: obj/SiPM_Time.o
Bi207: obj/Bi207.o
Analyser: obj/Analyser.o
BasicAnalyser: obj/BasicAnalyser.o
Rate: obj/Rate.o
Pulse: obj/Pulse.o
SIPM_PileUp: obj/SIPM_PileUp.o
Time_Pulse_Calibration: obj/Time_Pulse_Calibration.o

$(EXE): $(OBJS)
	${CC} -std=c++17 obj/$@.o -o $@ ../../lib/SignalDict/SignalDict.so Dic/AutoDict_vector_vector_std__pair_Signal__Signal____cxx.so ${CFLAGS} ${LIBS}

obj/%.o: src/%.c++ include/%.hh
	@mkdir -p obj
	${CC} -std=c++17 -c $< -o $@ ${CFLAGS}

clean :
	rm -f obj/*.o
	rm -f $(EXE)


