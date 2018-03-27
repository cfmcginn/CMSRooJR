#-Werror
#std-c++11

CXX = g++
CXXFLAGS = -Wall -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

INCLUDE = -I $(PWD)
ROOT = `root-config --cflags --glibs`
ROOUNF=/afs/cern.ch/work/c/cmcginn/public/Packages/RooUnfold/RooUnfold/

MKDIR_BIN = mkdir -p $(PWD)/bin
MKDIR_OUTPUT = mkdir -p $(PWD)/output
MKDIR_PDF = mkdir -p $(PWD)/pdfDir

all: mkdirBin mkdirPdf mkdirOutput buildAndTest

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

buildAndTest: src/buildAndTest.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -I $(ROOUNF) -L $(ROOUNF)  -lRooUnfold -lUnfold -o bin/buildAndTest.exe src/buildAndTest.C

clean:
	rm -f ./*~
	rm -f src/*~
	rm -f include/*~
	rm -f configs/*~
	rm -f ./#*#
	rm -f src/#*#
	rm -f include/#*#
	rm -f configs/#*#
	rm -f bin/*.exe
	rm -rf  bin