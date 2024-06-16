########################################################################
#
# Compilation for analysis
#
########################################################################

BASEDIR		:= ${shell cd ..; pwd}
ANADIR		:= ${shell expr "$(BASEDIR)" : "\(.*\)Analysis/.*"}Analysis

#-----------------------------------------------------------------------
# Color highlights
ifeq ($(findstring /sh,${SHELL}),/sh)
    RED     = "[31m"
    PURPLE  = "[35m"
    GREEN   = "[32m"
    BLUE    = "[34m"
    BOLD    = "[1m"
    UNDER   = "[4m"
    REVERSE = "[7m"
    END     = "[0m"
else
    ifeq ($(findstring /csh,${SHELL}),/csh)
        RED     = "\033[31m"
        PURPLE  = "\033[35m"
        GREEN   = "\033[32m"
        BLUE    = "\033[34m"
        BOLD    = "\033[1m"
        UNDER   = "\033[4m"
        REVERSE = "\033[7m"
        END     = "\033[0m"
    else
        RED     = ""
        PURPLE  = ""
        GREEN   = ""
        BLUE    = ""
        BOLD    = ""
        UNDER   = ""
        REVERSE = ""
        END     = ""
    endif
endif

#-----------------------------------------------------------------------
# check ROOT version
ROOTCONFIG	        := root-config
ROOTV	                := $(shell ${ROOTCONFIG} --version)
ROOT_VERSION	        := $(shell expr ${ROOTV} : '\(.*\)/.*')

# test if ROOT v5
ROOT5_TST	:= $(shell echo "${ROOT_VERSION} < 6" | bc)

# change from versions 5.xx to 6.xx
ROOTCINT	:= rootcling
ifeq ("${ROOT5_TST}","1")
      ROOTCINT        := rootcint
endif

#-----------------------------------------------------------------------
# Analysis class library
LIBNAME		:= FWisard

# Destination directories
LIBDIR          := ${BASEDIR}/lib
INCDIR          := ${BASEDIR}/include
BINDIR          := ${BASEDIR}/bin

DUMMY		:= ${shell mkdir -p obj} ${shell mkdir -p $(BINDIR)} \
		   ${shell mkdir -p $(LIBDIR)} ${shell mkdir -p $(INCDIR)}

#-----------------------------------------------------------------------
# Compiler

AR		:= ar -r
CPP		:= g++
LD		:= g++
SOFLAGS		:= -shared -lc 

#-----------------------------------------------------------------------
# Sources / Objets / Header files
EXE_SRC		:= ${shell \ls *.cxx}
EXE		:= ${EXE_SRC:%.cxx=${BINDIR}/%}
EXE_LNK		:= ${EXE_SRC:%.cxx=../../%}

# files for ROOT dictionary generation (rootcint)
DICT_SRC	:= ${LIBNAME}_dict.cpp
DICT_OBJ	:= ${DICT_SRC:%.cpp=obj/%.o}
ifeq (${ROOTCINT},"rootcling")
	DICT_PCM	:= 
else
	DICT_PCM	:= lib${LIBNAME}_rdict.pcm
endif

LIB_SRC		:= ${shell \ls *.cpp}
LIB_OBJ		:= ${LIB_SRC:%.cpp=obj/%.o} obj/${LIBNAME}_dict.o
LIB_INCL	:= ${shell \ls *.hh}
LIB_ICC		:= ${shell \ls icc/*.icc}

# Dynamic library
LIB		:= ${LIBDIR}/lib${LIBNAME}.so

# Compilation options

INCFLAGS	:= ${shell faster_process_inc}
CPPFLAGS	:= ${shell faster_process_cpp}
AUXFLAGS	:= ${shell faster_process_aux}
LDFLAGS		:= -Wl,-rpath,${LIBDIR} -L${LIBDIR} -l${LIBNAME} \
		   ${shell faster_process_lib}

DEBUG_FLAGS	:= ${DEBUG:%=-DGDEBUG=%}

#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Dependancies

all:		${LIB} ${EXE}

lib:		${LIB}

# -> main program
${BINDIR}/%:	%.cxx ${LIB}
		@echo ""
		@echo "${PURPLE}Executable: $@${END}"
		@echo ""
		${CPP} ${CPPFLAGS} ${AUXFLAGS} ${DEBUG_FLAGS} ${INCFLAGS} ${LDFLAGS} $< -o $@
		@\rm -f $(@:${BINDIR}/%=../../%)
		@ln -s $@ $(@:${BINDIR}/%=../../%)
		@echo ""

# -> dynamic library
${LIB}:		${LIB_OBJ}
		@echo ""
		@echo "${PURPLE}-> Library $@${END}"
		${LD} ${SOFLAGS} ${LIB_OBJ} -o ${LIB}
		@\cp -f $(LIB_INCL) $(INCDIR)
		@\cp -Rf icc $(INCDIR)
    ifeq ("${ROOTCINT}","rootcling")
		@\cp -fv $(DICT_PCM) $(LIBDIR)
		@\cp -fv lib$(LIBNAME).rootmap $(LIBDIR)
    endif
		@echo ""

# -> object files
obj/%.o:	%.cpp ${LIB_INCL} ${LIB_ICC} ${FASTER_PROCESS_LIB}/libFasterProcess.so
		@echo ""
		@echo "${GREEN}-> Compilation ${CPP} of $@${END}"
		${CPP} ${CPPFLAGS} ${AUXFLAGS} ${INCFLAGS} ${DEBUG_FLAGS} $< -c -o $@

# -> ROOT dictionary source
%_dict.cpp:	${LIB_INCL} ${LIBNAME}_linkdef.h
		@echo ""
		@echo "${BLUE}-> Generating Dictionary${END}"
    ifeq ("${ROOTCINT}","rootcling")
#		@echo "PCM file  = "${DICT_PCM}
#		@echo ""
#		(2019/07) the -I${INCDIR} (final include directory)
#		          does not seem to be needed
		${ROOTCINT} -f $@ -s lib${LIBNAME}.so -rml lib${LIBNAME}.so ${INCFLAGS} ${LIB_INCL} ${LIBNAME}_linkdef.h
    else
		@echo ""
		${ROOTCINT} -f $@ -c -p ${INCFLAGS} ${LIB_INCL} ${LIBNAME}_linkdef.h
    endif

# -> ROOT dictionary source
# %.cxx:	Signal.h	LinkDef.h
# 		@echo ""
# 		@echo "-> Generating Dictionary from $<"
# 		${ROOTCINT} -f $@ -c -p ${INCFLAGS} $<

# -> compilation cleaning
clean:
		@echo ""
		@echo "Cleaning up..."
		@set nonomatch; \rm -Rfv obj ${LIBDIR} ${INCDIR} ${BINDIR} *.pcm *.rootmap ${EXE_LNK} SignalDict.cxx SignalDict_rdict.pcm
		@echo ""

# -> compilation cleaning
clean_root:
		@echo ""
		@echo "Cleaning up ROOT compiled macro..."
		@set nonomatch; \rm -fv ../Root/*.d ../Root/*.so
		@set nonomatch; \rm -fv ../Calibration/*.d ../Calibration/*.so
		@set nonomatch; \rm -fv ../Root/*.bck ../Calibration/*.bck
		@echo ""

# -> compilation information
info:
		@echo ""
		@echo "BASEDIR:  "${BASEDIR}
		@echo ""
		@echo "ROOTCINT: "${ROOTCINT}
		@echo "DICT_PCM: "${DICT_PCM}
		@echo ""
		@echo "EXE_LNK:  "${EXE_LNK}
		@echo ""
