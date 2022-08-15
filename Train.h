#ifndef _TRAIN_H
#define _TRAIN_H
#include "SafeVector.h"
#include "MultiSequence.h"
#include "ScoreType.h"
#include "ProbabilisticModel.h"
#include <string>

using namespace std;

class Train
{
public:
    Train(int argc, char* argv[]);
    ~Train();
private:
    SafeVector<string> ParseParams(int argc, char *argv[]);
    void PrintParameters(const VF &initDistrib,
                         const VF &gapOpen, const VF &gapExtend, const VVF &emitPairs,
                         const VF &emitSingle, const char *filename);
    void DoAlignTrain(MultiSequence *sequences,
                      const ProbabilisticModel &model, VF &initDistrib, VF &gapOpen,
                      VF &gapExtend, VVF &emitPairs, VF &emitSingle, int first_step);
    void ReadParameters();

};

#endif
