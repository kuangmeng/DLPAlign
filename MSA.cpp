#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <iomanip>

#include "MSA.h"
#include "MSAClusterTree.h"
#include "Defaults.h"
#include "CalculateFeatures.h"

string parametersInputFilename = "";

bool enableAlignOrder = false;
int numConsistencyReps = 2;
int numPreTrainingReps = 0;
int numIterativeRefinementReps = 100;

float cutoff = 0;
int first_step = 0;
bool get_features = false;
int second_step = 0;
int third_step = 0;


VF initDistrib(NumMatrixTypes);
VF gapOpen(2 * NumInsertStates);
VF gapExtend(2 * NumInsertStates);
VVF emitPairs(256, VF(256, 1e-10));//subsititution matrix for double affine pair-HMM
VF emitSingle(256, 1e-5);//subsititution matrix for double affine pair-HMM

string alphabet = alphabetDefault;

///////////////////////////////
// global pair-HMM variables
//////////////////////////////

float g_gap_open1, g_gap_open2, g_gap_ext1, g_gap_ext2;
char *aminos, *bases, matrixtype[20] = "gonnet_160";
int subst_index[26];

double sub_matrix[26][26];      //subsititution matrix for global pair-HMM
double normalized_matrix[26][26];// be used to compute sequences' similarity score
int firstread = 0;		//this makes sure that matrices are read only once

float TEMPERATURE = 5;
int MATRIXTYPE = 160;

float GAPOPEN = 0;
float GAPEXT = 0;

//argument support
typedef struct
{
    char input[30];
    int matrix;
    int N;
    float T;
    float beta;
    char opt;			//can be 'P' or 'M'
    float gapopen;
    float gapext;
} argument_decl;

argument_decl argument;

extern void init_arguments();

MSA::MSA(int argc, char* argv[])
{

    //parse program parameters
    SafeVector<string> sequenceNames = ParseParams(argc, argv);

    //initialize arguments for partition function
    init_arguments();

    ReadParameters();

    //read the input sequences
    MultiSequence *sequences = new MultiSequence();
    assert(sequences);
    for (int i = 0; i < (int) sequenceNames.size(); i++)
    {
        sequences->LoadMFA(sequenceNames[i], true);
    }
    //allocate space for sequence weights
    this->seqsWeights = new int[sequences->GetNumSequences()];
    /**************************************
    // average percent identity
    // Divergent(<=25%) levelid = 0
    // Medium(25%-40%) levelid = 1
    // Similar(40%-70%) levelid = 2
    // High Similar(>70%) levelid = 3
    ***************************************/
    if(get_features == true)
    {
        CalculateFeatures *calculate = new CalculateFeatures(sequences, ProbabilisticModel(initDistrib, gapOpen, gapExtend, emitPairs,emitSingle));
        cout << calculate->getFeatures() << endl;
        return;
    }
    // now, we can perform the alignments and write them out
    MultiSequence *alignment = doAlign(sequences,
                                       ProbabilisticModel(initDistrib, gapOpen, gapExtend, emitPairs,emitSingle), first_step, second_step, third_step);

    alignment->WriteMFA(*alignOutFile);

    //release resources
    delete alignment;
    delete[] this->seqsWeights;
    delete sequences;
}

MSA::~MSA()
{
    /*close the output file*/
    if (alignOutFileName.length() > 0)
    {
        ((std::ofstream*) alignOutFile)->close();
    }
}


/////////////////////////////////////////////////////////////////
// doAlign()
//
// First computes all pairwise posterior probability matrices.
// Then, computes new parameters if training, or a final
// alignment, otherwise.
/////////////////////////////////////////////////////////////////

extern VF *ComputePostProbs(int a, int b, string seq1, string seq2);

MultiSequence* MSA::doAlign(MultiSequence *sequences,
                            const ProbabilisticModel &model, int first_step, int second_step, int third_step)
{
    assert(sequences);

    //get the number of sequences
    const int numSeqs = sequences->GetNumSequences();
    //create distance matrix
    VVF distances(numSeqs, VF(numSeqs, 0));
    //creat sparseMatrices
    SafeVector<SafeVector<SparseMatrix *> > sparseMatrices(numSeqs,
            SafeVector<SparseMatrix *>(numSeqs, NULL));

    // do all pairwise alignments for posterior probability matrices

    for (int a = 0; a < numSeqs - 1; a++)
    {
        for (int b = a + 1; b < numSeqs; b++)
        {
            Sequence *seq1 = sequences->GetSequence(a);
            Sequence *seq2 = sequences->GetSequence(b);

            //posterior probability matrix
            VF* posterior;

            if(first_step == 0)
            {
                //cerr << first_step << "Step 1: using Double Affine Pair-HMMs." << endl;
                // compute forward and backward probabilities
                VF *forward = model.ComputeForwardMatrix(seq1, seq2);
                assert(forward);
                VF *backward = model.ComputeBackwardMatrix(seq1, seq2);
                assert(backward);
                // compute posterior probability
                posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward,*backward);
                delete forward;
                delete backward;
            }

            else if(first_step == 1){
                //cerr << first_step << "Step 1: using Partition Function." << endl;
                posterior = ::ComputePostProbs(a, b, seq1->GetString(),seq2->GetString());
            }

 //divergent use combined model
            else if(first_step == 2)
            {
               // cerr << first_step << "Step 1: using MSAProbs Model." << endl;
                //double affine pair-HMM
                // compute forward and backward probabilities
                VF *forward = model.ComputeForwardMatrix(seq1, seq2);
                assert(forward);
                VF *backward = model.ComputeBackwardMatrix(seq1, seq2);
                assert(backward);
                // compute posterior probability
                VF *double_posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward,*backward);
                assert(double_posterior);
                delete forward;
                delete backward;

                // compute posterior probability
                posterior = ::ComputePostProbs(a, b, seq1->GetString(),seq2->GetString());
                assert(posterior);
                //combined model
                //merge probalign + local + probcons
                VF::iterator ptr1 = double_posterior->begin();
                VF::iterator ptr = posterior->begin();
                for (int i = 0; i <= seq1->GetLength(); i++)
                {
                    for (int j = 0; j <= seq2->GetLength(); j++)
                    {
                        float v1 = *ptr1;
                        float v3 = *ptr;
                        *ptr = sqrt((v1*v1 + v3*v3)/2);
                        ptr1++;
                        ptr++;
                    }
                }
                delete double_posterior;
            }
            //divergent use combined model
            else
            {
               // cerr << first_step << "Step 1: using Combination Model." << endl;
                //double affine pair-HMM
                // compute forward and backward probabilities
                VF *forward = model.ComputeForwardMatrix(seq1, seq2);
                assert(forward);
                VF *backward = model.ComputeBackwardMatrix(seq1, seq2);
                assert(backward);
                // compute posterior probability
                VF *double_posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward,*backward);
                assert(double_posterior);
                delete forward;
                delete backward;

                // compute posterior probability
                VF *global_posterior = ::ComputePostProbs(a, b, seq1->GetString(),seq2->GetString());
                assert(global_posterior);
                //local pair-HMM
                // compute forward and backward probabilities
                forward = model.ComputeForwardMatrix(seq1, seq2, false);
                assert(forward);
                backward = model.ComputeBackwardMatrix(seq1, seq2, false);
                assert(backward);
                // compute posterior probability
                posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward,*backward, false);
                assert(posterior);
                delete forward;
                delete backward;
                //combined model
                //merge probalign + local + probcons
                VF::iterator ptr1 = double_posterior->begin();
                VF::iterator ptr2 = global_posterior->begin();
                VF::iterator ptr = posterior->begin();
                for (int i = 0; i <= seq1->GetLength(); i++)
                {
                    for (int j = 0; j <= seq2->GetLength(); j++)
                    {
                        float v1 = *ptr1;
                        float v2 = *ptr2;
                        float v3 = *ptr;
                        *ptr = sqrt((v1*v1 + v2*v2 + v3*v3)/3);
                        ptr1++;
                        ptr2++;
                        ptr++;
                    }
                }
                delete double_posterior;
                delete global_posterior;
            }
            assert(posterior);
            // perform the pairwise sequence alignment
            pair<SafeVector<char> *, float> alignment = model.ComputeAlignment(
                        seq1->GetLength(), seq2->GetLength(), *posterior);

            //compute expected accuracy
            distances[a][b] = distances[b][a] = 1.0f - alignment.second
                                                / min(seq1->GetLength(), seq2->GetLength());

            // compute sparse representations
            sparseMatrices[a][b] = new SparseMatrix(seq1->GetLength(),
                                                    seq2->GetLength(), *posterior);
            sparseMatrices[b][a] = NULL;

            delete posterior;
            delete alignment.first;
        }
    }

    //create the guide tree
    this->tree = new MSAClusterTree(this, distances, numSeqs);
    this->tree->create(second_step);

    // perform the consistency transformation the desired number of times
    float* fweights = new float[numSeqs];
    for (int r = 0; r < numSeqs; r++)
    {
        fweights[r] = ((float) seqsWeights[r]) / INT_MULTIPLY;
        fweights[r] *= 10;
    }
    for (int r = 0; r < numConsistencyReps; r++)
    {
        SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices;
        if (third_step == 1){
            //cerr << "Step 3: using weighted Probabilistic Comsistency Transformation." << endl;
            newSparseMatrices = DoRelaxation(fweights, sequences, sparseMatrices);
        }else{
            //cerr << "Step 3: using unweighted Probabilistic Comsistency Transformation." << endl;
            newSparseMatrices = DoRelaxation(sequences, sparseMatrices);
        }

        // now replace the old posterior matrices
        for (int i = 0; i < numSeqs; i++)
        {
            for (int j = 0; j < numSeqs; j++)
            {
                delete sparseMatrices[i][j];
                sparseMatrices[i][j] = newSparseMatrices[i][j];
            }
        }
    }
    delete[] fweights;

    //compute the final multiple sequence alignment
    MultiSequence *finalAlignment = ComputeFinalAlignment(this->tree, sequences,
                                    sparseMatrices, model, first_step);
    //destroy the guide tree
    delete this->tree;
    this->tree = 0;

    // delete sparse matrices
    for (int a = 0; a < numSeqs - 1; a++)
    {
        for (int b = a + 1; b < numSeqs; b++)
        {
            delete sparseMatrices[a][b];
            delete sparseMatrices[b][a];
        }
    }

    return finalAlignment;
}

/////////////////////////////////////////////////////////////////
// GetInteger()
//
// Attempts to parse an integer from the character string given.
// Returns true only if no parsing error occurs.
/////////////////////////////////////////////////////////////////

bool GetInteger(char *data, int *val)
{
    char *endPtr;
    long int retVal;
    assert(val);
    errno = 0;
    retVal = strtol(data, &endPtr, 0);
    if (retVal == 0 && (errno != 0 || data == endPtr))
        return false;
    if (errno != 0 && (retVal == LONG_MAX || retVal == LONG_MIN))
        return false;
    if (retVal < (long) INT_MIN || retVal > (long) INT_MAX)
        return false;
    *val = (int) retVal;
    return true;
}


/////////////////////////////////////////////////////////////////
// ReadParameters()
//
// Read initial distribution, transition, and emission
// parameters from a file.
/////////////////////////////////////////////////////////////////

void MSA::ReadParameters()
{

    ifstream data;

    emitPairs = VVF(256, VF(256, 1e-10));
    emitSingle = VF(256, 1e-5);

    // read initial state distribution and transition parameters
    if (parametersInputFilename == string(""))
    {
        if (NumInsertStates == 1)
        {
            for (int i = 0; i < NumMatrixTypes; i++)
                initDistrib[i] = initDistrib1Default[i];
            for (int i = 0; i < 2 * NumInsertStates; i++)
                gapOpen[i] = gapOpen1Default[i];
            for (int i = 0; i < 2 * NumInsertStates; i++)
                gapExtend[i] = gapExtend1Default[i];
        }
        else if (NumInsertStates == 2)
        {
            for (int i = 0; i < NumMatrixTypes; i++)
                initDistrib[i] = initDistrib2Default[i];
            for (int i = 0; i < 2 * NumInsertStates; i++)
                gapOpen[i] = gapOpen2Default[i];
            for (int i = 0; i < 2 * NumInsertStates; i++)
                gapExtend[i] = gapExtend2Default[i];
        }
        else
        {
            cerr
                    << "ERROR: No default initial distribution/parameter settings exist"
                    << endl << "       for " << NumInsertStates
                    << " pairs of insert states.  Use --paramfile." << endl;
            exit(1);
        }

        alphabet = alphabetDefault;

        for (int i = 0; i < (int) alphabet.length(); i++)
        {
            emitSingle[(unsigned char) tolower(alphabet[i])] =
                emitSingleDefault[i];
            emitSingle[(unsigned char) toupper(alphabet[i])] =
                emitSingleDefault[i];
            for (int j = 0; j <= i; j++)
            {
                emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(
                            alphabet[j])] = emitPairsDefault[i][j];
                emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(
                            alphabet[j])] = emitPairsDefault[i][j];
                emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(
                            alphabet[j])] = emitPairsDefault[i][j];
                emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(
                            alphabet[j])] = emitPairsDefault[i][j];
                emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(
                            alphabet[i])] = emitPairsDefault[i][j];
                emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(
                            alphabet[i])] = emitPairsDefault[i][j];
                emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(
                            alphabet[i])] = emitPairsDefault[i][j];
                emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(
                            alphabet[i])] = emitPairsDefault[i][j];
            }
        }
    }
    else
    {
        data.open(parametersInputFilename.c_str());
        if (data.fail())
        {
            cerr << "ERROR: Unable to read parameter file: "
                 << parametersInputFilename << endl;
            exit(1);
        }

        string line[3];
        for (int i = 0; i < 3; i++)
        {
            if (!getline(data, line[i]))
            {
                cerr
                        << "ERROR: Unable to read transition parameters from parameter file: "
                        << parametersInputFilename << endl;
                exit(1);
            }
        }
        istringstream data2;
        data2.clear();
        data2.str(line[0]);
        for (int i = 0; i < NumMatrixTypes; i++)
            data2 >> initDistrib[i];
        data2.clear();
        data2.str(line[1]);
        for (int i = 0; i < 2 * NumInsertStates; i++)
            data2 >> gapOpen[i];
        data2.clear();
        data2.str(line[2]);
        for (int i = 0; i < 2 * NumInsertStates; i++)
            data2 >> gapExtend[i];

        if (!getline(data, line[0]))
        {
            cerr << "ERROR: Unable to read alphabet from scoring matrix file: "
                 << parametersInputFilename << endl;
            exit(1);
        }

        // read alphabet as concatenation of all characters on alphabet line
        alphabet = "";
        string token;
        data2.clear();
        data2.str(line[0]);
        while (data2 >> token)
            alphabet += token;

        for (int i = 0; i < (int) alphabet.size(); i++)
        {
            for (int j = 0; j <= i; j++)
            {
                float val;
                data >> val;
                emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(
                            alphabet[j])] = val;
                emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(
                            alphabet[j])] = val;
                emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(
                            alphabet[j])] = val;
                emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(
                            alphabet[j])] = val;
                emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(
                            alphabet[i])] = val;
                emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(
                            alphabet[i])] = val;
                emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(
                            alphabet[i])] = val;
                emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(
                            alphabet[i])] = val;
            }
        }

        for (int i = 0; i < (int) alphabet.size(); i++)
        {
            float val;
            data >> val;
            emitSingle[(unsigned char) tolower(alphabet[i])] = val;
            emitSingle[(unsigned char) toupper(alphabet[i])] = val;
        }
        data.close();
    }
}

/////////////////////////////////////////////////////////////////
// ParseParams()
//
// Parse all command-line options.
/////////////////////////////////////////////////////////////////

void MSA::printUsage()
{
    cerr << "Usage:" << endl
         << "       -G           -- Get features for Classification on Posterior Probability Matrix." << endl
         << "       -s1 <0|1|2>  -- Use different Posterior Probability Matrix to do MSA." << endl
         << "       -s2 <0|1>    -- Use different Clustering Methods to generate Guide Tree." << endl
         << "       -s3 <0|1>    -- Use different Consistency Transformation Matrix." << endl << endl;
}

SafeVector<string> MSA::ParseParams(int argc, char **argv)
{
    if (argc < 2)
    {
        printUsage();
        exit(1);
    }
    SafeVector<string> sequenceNames;
    int tempInt;

    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-')
        {
            //help
            if (!strcmp(argv[i], "-help") || !strcmp(argv[i], "-?"))
            {
                printUsage();
                exit(1);
                //output file name
            }

            else if (!strcmp(argv[i], "-G")
                     || !strcmp(argv[i], "-get"))
            {
                get_features = true;
            }

            else if (!strcmp(argv[i], "-s1")
                     || !strcmp(argv[i], "-step1"))
            {
                if (!GetInteger(argv[++i], &tempInt) && i < argc - 1)
                {
                    first_step = 3;
                }
                else
                {
                        if(tempInt == 0)
                        {
                            first_step = 0;
                        }
                        else if(tempInt == 1)
                        {
                            first_step = 1;
                        }
                        else if(tempInt == 2)
                        {
                            first_step = 2;
                        }
                        else{
                            first_step = 3;
                        }
                    
                }
            }
            else if (!strcmp(argv[i], "-s2")
                     || !strcmp(argv[i], "-step2"))
            {
                if (!GetInteger(argv[++i], &tempInt) && i < argc - 1)
                {
                    second_step = 0;
                }
                else
                {
                    if (tempInt > 1 || tempInt < 0)
                    {
                        second_step = 0;
                    }
                    else
                    {
                        if(tempInt == 0)
                        {
                            second_step = 0;
                        }
                        else
                        {
                            second_step = 1;
                        }
                    }
                }
            }
            else if (!strcmp(argv[i], "-s3")
                     || !strcmp(argv[i], "-step3"))
            {
                if (!GetInteger(argv[++i], &tempInt) && i < argc - 1)
                {
                    third_step = 0;
                }
                else
                {
                    if (tempInt > 1 || tempInt < 0)
                    {
                        third_step = 0;
                    }
                    else
                    {
                        if(tempInt == 0)
                        {
                            third_step = 0;
                        }
                        else
                        {
                            third_step = 1;
                        }
                    }
                }
            }
            else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--outfile"))
            {
                if (i < argc - 1)
                {
                    alignOutFileName = argv[++i];	//get the file name
                }
                else
                {
                    cerr << "ERROR: String expected for option " << argv[i]
                         << endl;
                    exit(1);
                }
                // number of consistency transformations
            }

            // alignment order
            else if (!strcmp(argv[i], "-a")
                     || !strcmp(argv[i], "--alignment-order"))
            {
                enableAlignOrder = true;
            }
            // bad arguments
            else
            {
                cerr << "ERROR: Unrecognized option: " << argv[i] << endl;
                exit(1);
            }
        }
        else
        {
            sequenceNames.push_back(string(argv[i]));
        }
    }

    if (alignOutFileName.length() == 0)
    {
        alignOutFile = &std::cout;
    }
    else
    {
        cerr << "Open the output file " << alignOutFileName << endl;
        alignOutFile = new ofstream(alignOutFileName.c_str(),
                                    ios::binary | ios::out | ios::trunc);
    }
    return sequenceNames;
}

/////////////////////////////////////////////////////////////////
// ProcessTree()
//
// Process the tree recursively.  Returns the aligned sequences
// corresponding to a node or leaf of the tree.
/////////////////////////////////////////////////////////////////

MultiSequence* MSA::ProcessTree(TreeNode *tree, MultiSequence *sequences,
                                const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                const ProbabilisticModel &model)
{

    MultiSequence *result;

    // check if this is a node of the alignment tree
    //if (tree->GetSequenceLabel() == -1){
    if (tree->leaf == NODE)
    {
        MultiSequence *alignLeft = ProcessTree(tree->left, sequences,
                                               sparseMatrices, model);
        MultiSequence *alignRight = ProcessTree(tree->right, sequences,
                                                sparseMatrices, model);

        assert(alignLeft);
        assert(alignRight);

        result = AlignAlignments(alignLeft, alignRight, sparseMatrices, model);
        assert(result);

        delete alignLeft;
        delete alignRight;
    }

    // otherwise, this is a leaf of the alignment tree
    else
    {
        result = new MultiSequence();
        assert(result);
        result->AddSequence(sequences->GetSequence(tree->idx)->Clone());
    }

    return result;
}

/////////////////////////////////////////////////////////////////
// ComputeFinalAlignment()
//
// Compute the final alignment by calling ProcessTree(), then
// performing iterative refinement as needed.
/////////////////////////////////////////////////////////////////

MultiSequence* MSA::ComputeFinalAlignment(MSAGuideTree*tree,
        MultiSequence *sequences,
        const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
        const ProbabilisticModel &model, int first_step)
{

    MultiSequence *alignment = ProcessTree(tree->getRoot(), sequences,
                                           sparseMatrices, model);

    SafeVector<int> oldOrdering;
    int numSeqs = alignment->GetNumSequences();
    if (enableAlignOrder)
    {
        for (int i = 0; i < numSeqs; i++)
            oldOrdering.push_back(alignment->GetSequence(i)->GetSortLabel());
        alignment->SaveOrdering();
        enableAlignOrder = false;
    }

    //DoIterativeRefinement() return 1,2: this refinement unsuccessful
    if(first_step == 2 || numSeqs >= 150 ) numIterativeRefinementReps=10;
    int ineffectiveness = 0;
    for (int i = 0; i < numIterativeRefinementReps; i++)
    {
        int flag = DoIterativeRefinement(sparseMatrices, model, alignment);
        if(numSeqs > 25 && numSeqs < 150 && first_step < 2)
        {
            if(flag > 0)
            {
                if(numIterativeRefinementReps < 4*numSeqs)
                    numIterativeRefinementReps ++;
                if(flag == 1) ineffectiveness ++;
            }
            //else ineffectiveness = 0;
            if(ineffectiveness > 2*numSeqs && i>100) break;
        }
    }
    return alignment;
}

/////////////////////////////////////////////////////////////////
// AlignAlignments()
//
// Returns the alignment of two MultiSequence objects.
/////////////////////////////////////////////////////////////////

MultiSequence* MSA::AlignAlignments(MultiSequence *align1,
                                    MultiSequence *align2,
                                    const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                    const ProbabilisticModel &model)
{

#if 0
    VF *posterior = model.BuildPosterior (align1, align2, sparseMatrices, cutoff);
#else
    VF *posterior = model.BuildPosterior(getSeqsWeights(), align1, align2,
                                         sparseMatrices, cutoff);
#endif
    // compute an "accuracy" measure for the MSA before refinement

    pair<SafeVector<char> *, float> alignment;
    //perform alignment
    alignment = model.ComputeAlignment(align1->GetSequence(0)->GetLength(),
                                       align2->GetSequence(0)->GetLength(), *posterior);

    delete posterior;


    // now build final alignment
    MultiSequence *result = new MultiSequence();
    for (int i = 0; i < align1->GetNumSequences(); i++)
        result->AddSequence(
            align1->GetSequence(i)->AddGaps(alignment.first, 'X'));
    for (int i = 0; i < align2->GetNumSequences(); i++)
        result->AddSequence(
            align2->GetSequence(i)->AddGaps(alignment.first, 'Y'));
    if (!enableAlignOrder)
        result->SortByLabel();

    // free temporary alignment
    delete alignment.first;

    return result;
}


/////////////////////////////////////////////////////////////////
// DoRelaxation()
//
// Performs one round of the weighted probabilistic consistency transformation.
/////////////////////////////////////////////////////////////////

SafeVector<SafeVector<SparseMatrix *> > MSA::DoRelaxation(float* seqsWeights,
        MultiSequence *sequences,
        SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices)
{
    const int numSeqs = sequences->GetNumSequences();

    SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices(numSeqs,
            SafeVector<SparseMatrix *>(numSeqs, NULL));

    // for every pair of sequences
    for (int i = 0; i < numSeqs; i++)
    {
        float wi = seqsWeights[i];
        for (int j = i + 1; j < numSeqs; j++)
        {
            float wj = seqsWeights[j];
            Sequence *seq1 = sequences->GetSequence(i);
            Sequence *seq2 = sequences->GetSequence(j);
            // get the original posterior matrix
            VF *posteriorPtr = sparseMatrices[i][j]->GetPosterior();
            assert(posteriorPtr);
            VF &posterior = *posteriorPtr;

            const int seq1Length = seq1->GetLength();
            const int seq2Length = seq2->GetLength();

            // contribution from the summation where z = x and z = y
            float w = wi * wi * wj + wi * wj * wj;
            float sumW = w;
            for (int k = 0; k < (seq1Length + 1) * (seq2Length + 1); k++)
            {
                posterior[k] = w * posterior[k];
            }

            // contribution from all other sequences
            for (int k = 0; k < numSeqs; k++)
            {
                if (k != i && k != j)
                {
                    float wk = seqsWeights[k];
                    float w = wi * wj * wk;
                    sumW += w;
                    if (k < i)
                        Relax1(w, sparseMatrices[k][i], sparseMatrices[k][j],
                               posterior);
                    else if (k > i && k < j)
                        Relax(w, sparseMatrices[i][k], sparseMatrices[k][j],
                              posterior);
                    else
                    {
                        SparseMatrix *temp =
                            sparseMatrices[j][k]->ComputeTranspose();
                        Relax(w, sparseMatrices[i][k], temp, posterior);
                        delete temp;
                    }
                }
            }
            //cerr<<"sumW "<<sumW<<endl;
            for (int k = 0; k < (seq1Length + 1) * (seq2Length + 1); k++)
            {
                posterior[k] /= sumW;
            }
            // mask out positions not originally in the posterior matrix
            SparseMatrix *matXY = sparseMatrices[i][j];
            for (int y = 0; y <= seq2Length; y++)
                posterior[y] = 0;
            for (int x = 1; x <= seq1Length; x++)
            {
                SafeVector<PIF>::iterator XYptr = matXY->GetRowPtr(x);
                SafeVector<PIF>::iterator XYend = XYptr + matXY->GetRowSize(x);
                VF::iterator base = posterior.begin() + x * (seq2Length + 1);
                int curr = 0;
                while (XYptr != XYend)
                {

                    // zero out all cells until the first filled column
                    while (curr < XYptr->first)
                    {
                        base[curr] = 0;
                        curr++;
                    }

                    // now, skip over this column
                    curr++;
                    ++XYptr;
                }

                // zero out cells after last column
                while (curr <= seq2Length)
                {
                    base[curr] = 0;
                    curr++;
                }
            }

            // save the new posterior matrix
            newSparseMatrices[i][j] = new SparseMatrix(seq1->GetLength(),
                    seq2->GetLength(), posterior);
            newSparseMatrices[j][i] = NULL;

            delete posteriorPtr;
        }
    }

    return newSparseMatrices;
}
/////////////////////////////////////////////////////////////////
// DoRelaxation()
//
// Performs one round of the unweighted probabilistic consistency transformation.
/////////////////////////////////////////////////////////////////

SafeVector<SafeVector<SparseMatrix *> > MSA::DoRelaxation(MultiSequence *sequences,
        SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices)
{
    const int numSeqs = sequences->GetNumSequences();

    SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices(numSeqs,
            SafeVector<SparseMatrix *>(numSeqs, NULL));

    // for every pair of sequences
    for (int i = 0; i < numSeqs; i++)
    {
        for (int j = i + 1; j < numSeqs; j++)
        {
            Sequence *seq1 = sequences->GetSequence(i);
            Sequence *seq2 = sequences->GetSequence(j);
            // get the original posterior matrix
            VF *posteriorPtr = sparseMatrices[i][j]->GetPosterior();
            assert(posteriorPtr);
            VF &posterior = *posteriorPtr;

            const int seq1Length = seq1->GetLength();
            const int seq2Length = seq2->GetLength();

            for (int k = 0; k < (seq1Length + 1) * (seq2Length + 1); k++)
            {
                posterior[k] += posterior[k];
            }

            // contribution from all other sequences
            for (int k = 0; k < numSeqs; k++)
            {
                if (k != i && k != j)
                {
                    if (k < i)
                        Relax1(sparseMatrices[k][i], sparseMatrices[k][j],posterior);
                    else if (k > i && k < j)
                        Relax(sparseMatrices[i][k], sparseMatrices[k][j],posterior);
                    else
                    {
                        SparseMatrix *temp = sparseMatrices[j][k]->ComputeTranspose();
                        Relax(sparseMatrices[i][k], temp, posterior);
                        delete temp;
                    }
                }
            }
            for (int k = 0; k < (seq1Length + 1) * (seq2Length + 1); k++)
            {
                posterior[k] /= numSeqs;
            }
            // mask out positions not originally in the posterior matrix
            SparseMatrix *matXY = sparseMatrices[i][j];
            for (int y = 0; y <= seq2Length; y++)
                posterior[y] = 0;
            for (int x = 1; x <= seq1Length; x++)
            {
                SafeVector<PIF>::iterator XYptr = matXY->GetRowPtr(x);
                SafeVector<PIF>::iterator XYend = XYptr + matXY->GetRowSize(x);
                VF::iterator base = posterior.begin() + x * (seq2Length + 1);
                int curr = 0;
                while (XYptr != XYend)
                {
                    // zero out all cells until the first filled column
                    while (curr < XYptr->first)
                    {
                        base[curr] = 0;
                        curr++;
                    }
                    // now, skip over this column
                    curr++;
                    ++XYptr;
                }
                // zero out cells after last column
                while (curr <= seq2Length)
                {
                    base[curr] = 0;
                    curr++;
                }
            }
            // save the new posterior matrix
            newSparseMatrices[i][j] = new SparseMatrix(seq1->GetLength(),
                    seq2->GetLength(), posterior);
            newSparseMatrices[j][i] = NULL;
            delete posteriorPtr;

        }
    }

    return newSparseMatrices;
}


/////////////////////////////////////////////////////////////////
// Relax()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

//with weight
void MSA::Relax(float weight, SparseMatrix *matXZ, SparseMatrix *matZY,
                VF &posterior)
{

    assert(matXZ);
    assert(matZY);

    int lengthX = matXZ->GetSeq1Length();
    int lengthY = matZY->GetSeq2Length();
    assert(matXZ->GetSeq2Length() == matZY->GetSeq1Length());

    // for every x[i]
    for (int i = 1; i <= lengthX; i++)
    {
        SafeVector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
        SafeVector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

        VF::iterator base = posterior.begin() + i * (lengthY + 1);

        // iterate through all x[i]-z[k]
        while (XZptr != XZend)
        {
            SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
            SafeVector<PIF>::iterator ZYend = ZYptr
                                              + matZY->GetRowSize(XZptr->first);
            const float XZval = XZptr->second;

            // iterate through all z[k]-y[j]
            while (ZYptr != ZYend)
            {
                base[ZYptr->first] += weight * XZval * ZYptr->second;
                ZYptr++;
            }
            XZptr++;
        }
    }
}

//without weight
void MSA::Relax(SparseMatrix *matXZ, SparseMatrix *matZY,
                VF &posterior)
{

    assert(matXZ);
    assert(matZY);

    int lengthX = matXZ->GetSeq1Length();
    int lengthY = matZY->GetSeq2Length();
    assert(matXZ->GetSeq2Length() == matZY->GetSeq1Length());

    // for every x[i]
    for (int i = 1; i <= lengthX; i++)
    {
        SafeVector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
        SafeVector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

        VF::iterator base = posterior.begin() + i * (lengthY + 1);

        // iterate through all x[i]-z[k]
        while (XZptr != XZend)
        {
            SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
            SafeVector<PIF>::iterator ZYend = ZYptr
                                              + matZY->GetRowSize(XZptr->first);
            const float XZval = XZptr->second;

            // iterate through all z[k]-y[j]
            while (ZYptr != ZYend)
            {
                base[ZYptr->first] += XZval * ZYptr->second;
                ZYptr++;
            }
            XZptr++;
        }
    }
}

/////////////////////////////////////////////////////////////////
// Relax1()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

//with weight
void MSA::Relax1(float weight, SparseMatrix *matZX, SparseMatrix *matZY,VF &posterior)
{

    assert(matZX);
    assert(matZY);

    int lengthZ = matZX->GetSeq1Length();
    int lengthY = matZY->GetSeq2Length();

    // for every z[k]
    for (int k = 1; k <= lengthZ; k++)
    {
        SafeVector<PIF>::iterator ZXptr = matZX->GetRowPtr(k);
        SafeVector<PIF>::iterator ZXend = ZXptr + matZX->GetRowSize(k);

        // iterate through all z[k]-x[i]
        while (ZXptr != ZXend)
        {
            SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(k);
            SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(k);
            const float ZXval = ZXptr->second;
            VF::iterator base = posterior.begin()
                                + ZXptr->first * (lengthY + 1);

            // iterate through all z[k]-y[j]
            while (ZYptr != ZYend)
            {
                base[ZYptr->first] += weight * ZXval * ZYptr->second;
                ZYptr++;
            }
            ZXptr++;
        }
    }
}

//without weight
void MSA::Relax1(SparseMatrix *matZX, SparseMatrix *matZY,VF &posterior)
{

    assert(matZX);
    assert(matZY);

    int lengthZ = matZX->GetSeq1Length();
    int lengthY = matZY->GetSeq2Length();

    // for every z[k]
    for (int k = 1; k <= lengthZ; k++)
    {
        SafeVector<PIF>::iterator ZXptr = matZX->GetRowPtr(k);
        SafeVector<PIF>::iterator ZXend = ZXptr + matZX->GetRowSize(k);

        // iterate through all z[k]-x[i]
        while (ZXptr != ZXend)
        {
            SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(k);
            SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(k);
            const float ZXval = ZXptr->second;
            VF::iterator base = posterior.begin()
                                + ZXptr->first * (lengthY + 1);

            // iterate through all z[k]-y[j]
            while (ZYptr != ZYend)
            {
                base[ZYptr->first] += ZXval * ZYptr->second;
                ZYptr++;
            }
            ZXptr++;
        }
    }
}

/////////////////////////////////////////////////////////////////
// DoIterativeRefinement()
//
// Performs a single round of randomized partionining iterative
// refinement.
// return 0: successful refinement, 1: ineffective refinement, 2: fail to partion two profiles
/////////////////////////////////////////////////////////////////

int MSA::DoIterativeRefinement(
    const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
    const ProbabilisticModel &model, MultiSequence* &alignment)
{
    set<int> groupOne, groupTwo;
    int numSeqs = alignment->GetNumSequences();
    int i;
    // create two separate groups
    for (i = 0; i < numSeqs; i++)
    {
        int index = rand();
        if (index % 2)
        {
            groupOne.insert(i);
        }
        else
        {
            groupTwo.insert(i);
        }
    }
    if (groupOne.empty() || groupTwo.empty()) return 2;

    // project into the two groups
    MultiSequence *groupOneSeqs = alignment->Project(groupOne);
    assert(groupOneSeqs);
    MultiSequence *groupTwoSeqs = alignment->Project(groupTwo);
    assert(groupTwoSeqs);

//no weight profile-profile for refinement
#if 1
    VF *posterior = model.BuildPosterior (groupOneSeqs, groupTwoSeqs, sparseMatrices, cutoff);
#else
    VF *posterior = model.BuildPosterior(getSeqsWeights(), groupOneSeqs, groupTwoSeqs,
                                         sparseMatrices, cutoff);
#endif
    // compute an "accuracy" for the currrent MSA before refinement
    SafeVector<SafeVector<char>::iterator> oldOnePtrs(groupOne.size());
    SafeVector<SafeVector<char>::iterator> oldTwoPtrs(groupTwo.size());
    i=0;
    for (set<int>::const_iterator iter = groupOne.begin();
            iter != groupOne.end(); ++iter)
    {
        oldOnePtrs[i++] = alignment->GetSequence(*iter)->GetDataPtr();
    }
    i=0;
    for (set<int>::const_iterator iter = groupTwo.begin();
            iter != groupTwo.end(); ++iter)
    {
        oldTwoPtrs[i++] = alignment->GetSequence(*iter)->GetDataPtr();
    }

    VF &posteriorArr = *posterior;
    int oldLength = alignment->GetSequence(0)->GetLength();
    int groupOneindex=0;
    int groupTwoindex=0;
    float accuracy_before = 0;
    int j;
    for (i = 1; i <= oldLength; i++)
    {
        // check to see if there is a gap in every sequence of the set
        bool foundOne = false;
        for (j = 0; !foundOne && j < (int) groupOne.size(); j++)
            foundOne = (oldOnePtrs[j][i] != '-');
        // if not, then this column counts towards the sequence length
        if (foundOne) groupOneindex ++;
        bool foundTwo = false;
        for (j = 0; !foundTwo && j < (int) groupTwo.size(); j++)
            foundTwo = (oldTwoPtrs[j][i] != '-');
        if (foundTwo) groupTwoindex ++;
        if(foundOne && foundTwo) accuracy_before +=
                posteriorArr[groupOneindex * (groupTwoSeqs->GetSequence(0)->GetLength() + 1) + groupTwoindex];
    }

    pair<SafeVector<char> *, float> refinealignment;
    //perform alignment
    refinealignment = model.ComputeAlignment(groupOneSeqs->GetSequence(0)->GetLength(),
                      groupTwoSeqs->GetSequence(0)->GetLength(), *posterior);
    delete posterior;
    // now build final alignment
    MultiSequence *result = new MultiSequence();
    for (int i = 0; i < groupOneSeqs->GetNumSequences(); i++)
        result->AddSequence(
            groupOneSeqs->GetSequence(i)->AddGaps(refinealignment.first, 'X'));
    for (int i = 0; i < groupTwoSeqs->GetNumSequences(); i++)
        result->AddSequence(
            groupTwoSeqs->GetSequence(i)->AddGaps(refinealignment.first, 'Y'));
    // free temporary alignment
    delete refinealignment.first;
    delete alignment;
    alignment = result;
    delete groupOneSeqs;
    delete groupTwoSeqs;
    if(accuracy_before == refinealignment.second) return 1;
    else return 0;
}


void MSA::DoIterativeRefinementTreeNode(
    const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
    const ProbabilisticModel &model, MultiSequence* &alignment,
    int nodeIndex)
{
    set<int> groupOne, groupTwo;
    int numSeqs = alignment->GetNumSequences();

    vector<bool> inGroup1;
    inGroup1.resize(numSeqs);
    for (int i = 0; i < numSeqs; i++)
    {
        inGroup1[i] = false;
    }

    AlignmentOrder* orders = this->tree->getAlignOrders();
    AlignmentOrder* order = &orders[nodeIndex];
    for (int i = 0; i < order->leftNum; i++)
    {
        int si = order->leftLeafs[i];
        inGroup1[si] = true;
    }
    for (int i = 0; i < order->rightNum; i++)
    {
        int si = order->rightLeafs[i];
        inGroup1[si] = true;
    }
    // create two separate groups
    for (int i = 0; i < numSeqs; i++)
    {
        if (inGroup1[i])
        {
            groupOne.insert(i);
        }
        else
        {
            groupTwo.insert(i);
        }
    }
    if (groupOne.empty() || groupTwo.empty())
        return;

    // project into the two groups
    MultiSequence *groupOneSeqs = alignment->Project(groupOne);
    assert(groupOneSeqs);
    MultiSequence *groupTwoSeqs = alignment->Project(groupTwo);
    assert(groupTwoSeqs);
    delete alignment;

    // realign
    alignment = AlignAlignments(groupOneSeqs, groupTwoSeqs, sparseMatrices,
                                model);

    delete groupOneSeqs;
    delete groupTwoSeqs;
}

int main(int argc, char* argv[])
{
    MSA msa(argc, argv);

    return 0;
}

