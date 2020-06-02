
CXXOBJS = MSAGuideTree.o MSAClusterTree.o MSAPartProbs.o MSAReadMatrix.o CalculateFeatures.o MSA.o 

TRAINOBJ = Train.o 

CXX = g++
COMMON_FLAGS = -std=c++11 -O3 -Wall -funroll-loops -I . -I /usr/include
CXXFLAGS = $(COMMON_FLAGS)

EXEC = dlpalign_main

TRAIN_EXEC = dlpalign_train

all: $(CXXOBJS)
	@$(CXX) $(CXXFLAGS) -o $(EXEC) $(CXXOBJS) $(NVCCOBJS) $(NVCCLIBS)
	@strip $(EXEC)
	@echo "Compile Finished."
	@rm -rf *.o

clean:
	@rm -rf *.o $(EXEC) $(TRAIN_EXEC)

train: $(TRAINOBJ)
	@$(CXX) $(CXXFLAGS) -o $(TRAIN_EXEC) $(TRAINOBJ) $(NVCCOBJS) $(NVCCLIBS)
	@strip $(TRAIN_EXEC)
	@echo "Compile Finished for Train."
	@rm -rf *.o
	

