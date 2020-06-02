#ifndef _MSA_H
#define _MSA_H
#include "MSAGuideTree.h"
#include "SafeVector.h"
#include "MultiSequence.h"
#include "ScoreType.h"
#include "ProbabilisticModel.h"
#include "SparseMatrix.h"
#include <string>

using namespace std;

class MSAGuideTree;
struct TreeNode;
class MSA
{
public:
    MSA(int argc, char* argv[]);
    ~MSA();
    MSAGuideTree* getGuideTree()
    {
        return tree;
    }
    int * getSeqsWeights()
    {
        return seqsWeights;
    }
private:
    //print usage
    void printUsage();
    //do multiple sequence alignment
    void doAlign();
    //for sequence weights
    void createSeqsWeights(int seqsNum);
    void releaseSeqsWeights();

    //weights of sequences
    int * seqsWeights;
    //guide tree
    MSAGuideTree* tree;
    //output file
    string alignOutFileName;
    std::ostream* alignOutFile;
private:
    SafeVector<string> ParseParams(int argc, char *argv[]);

    SafeVector<string> PostProbsParseParams(int argc, char **argv);
    MultiSequence *doAlign(MultiSequence *sequence,
                           const ProbabilisticModel &model, int first_step, int second_step, int third_step);
    void ReadParameters();
    MultiSequence* ProcessTree(TreeNode *tree, MultiSequence *sequences,
                               const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                               const ProbabilisticModel &model);
    MultiSequence *ComputeFinalAlignment(MSAGuideTree *tree,
                                         MultiSequence *sequences,
                                         const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                         const ProbabilisticModel &model,int first_step);
    MultiSequence *AlignAlignments(MultiSequence *align1, MultiSequence *align2,
                                   const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                   const ProbabilisticModel &model);
    SafeVector<SafeVector<SparseMatrix *> > DoRelaxation(float* seqsWeights,
            MultiSequence *sequences,
            SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);
    SafeVector<SafeVector<SparseMatrix *> > DoRelaxation(MultiSequence *sequences,
            SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);

    void Relax(float weight, SparseMatrix *matXZ, SparseMatrix *matZY,VF &posterior);//weight
    void Relax1(float weight, SparseMatrix *matXZ, SparseMatrix *matZY,VF &posterior);//weight
    void Relax(SparseMatrix *matXZ, SparseMatrix *matZY,VF &posterior);//unweight
    void Relax1(SparseMatrix *matXZ, SparseMatrix *matZY,VF &posterior);//unweight

    int DoIterativeRefinement(
        const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
        const ProbabilisticModel &model, MultiSequence* &alignment);
    void DoIterativeRefinementTreeNode(
        const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
        const ProbabilisticModel &model, MultiSequence* &alignment,
        int nodeIndex);
};

#endif
