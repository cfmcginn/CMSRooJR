#-Werror
#std-c++11

CXX = g++
CXXFLAGS = -Wall -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

INCLUDE = -I $(PWD)
ROOT = `root-config --cflags --glibs`
ROOUNF=/afs/cern.ch/work/c/cmcginn/public/Packages/RooUnfold/RooUnfold/

MKDIR_BIN = mkdir -p $(PWD)/bin
MKDIR_OUTPUT = mkdir -p $(PWD)/output
MKDIR_PDF = mkdir -p $(PWD)/pdfDir

#all: mkdirBin mkdirPdf mkdirOutput makeFullRAAHist_FromTree valgrindTest makeFirstRAAHist_FromTree makeFinalRAAHist_FromTree

all: mkdirBin mkdirPdf mkdirOutput bin/recreateV2V3Tree.exe bin/recreateV2V3TreeHist.exe bin/plotRecreateV2V3Tree.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

buildAndTest: src/buildAndTest.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -I $(ROOUNF) -L $(ROOUNF)  -lRooUnfold -o bin/buildAndTest.exe src/buildAndTest.C

bin/v2AndV3.exe: src/v2AndV3.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -I $(ROOUNF) -L $(ROOUNF)  -lRooUnfold -o bin/v2AndV3.exe src/v2AndV3.C

bin/recreateV2V3_PosNeg.exe: src/recreateV2V3_PosNeg.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -o bin/recreateV2V3_PosNeg.exe src/recreateV2V3_PosNeg.C

bin/recreateV2V3.exe: src/recreateV2V3.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -o bin/recreateV2V3.exe src/recreateV2V3.C

bin/recreateV2V3Tree.exe: src/recreateV2V3Tree.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -o bin/recreateV2V3Tree.exe src/recreateV2V3Tree.C

bin/recreateV2V3TreeHist.exe: src/recreateV2V3TreeHist.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -o bin/recreateV2V3TreeHist.exe src/recreateV2V3TreeHist.C

bin/plotRecreateV2V3Tree.exe: src/plotRecreateV2V3Tree.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -o bin/plotRecreateV2V3Tree.exe src/plotRecreateV2V3Tree.C

bin/makeFullRAAHist_FromTree.exe: src/makeFullRAAHist_FromTree.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -I $(ROOUNF) -L $(ROOUNF)  -lRooUnfold -lpthread -o bin/makeFullRAAHist_FromTree.exe src/makeFullRAAHist_FromTree.C

bin/makeFirstRAAHist_FromTree.exe: src/makeFirstRAAHist_FromTree.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -I $(ROOUNF) -L $(ROOUNF)  -lRooUnfold -lpthread -o bin/makeFirstRAAHist_FromTree.exe src/makeFirstRAAHist_FromTree.C

bin/makePlotValidation_FromTree.exe: src/makePlotValidation_FromTree.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -I $(ROOUNF) -L $(ROOUNF)  -lRooUnfold -lpthread -o bin/makePlotValidation_FromTree.exe src/makePlotValidation_FromTree.C

makeFirstRawRAAHist_FromTree: src/makeFirstRawRAAHist_FromTree.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -I $(ROOUNF) -L $(ROOUNF)  -lRooUnfold -lpthread -o bin/makeFirstRawRAAHist_FromTree.exe src/makeFirstRawRAAHist_FromTree.C

bin/makeFinalRAAHist_FromTree.exe: src/makeFinalRAAHist_FromTree.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -I $(ROOUNF) -L $(ROOUNF)  -lRooUnfold -lpthread -o bin/makeFinalRAAHist_FromTree.exe src/makeFinalRAAHist_FromTree.C

bin/rc.exe: src/rc.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -o bin/rc.exe src/rc.C

valgrindTest: src/valgrindTest.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o bin/valgrindTest.exe src/valgrindTest.C

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
	rm -rf bin