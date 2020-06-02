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
#include <string>
#include <string.h>
#include "Train.h"
#include "Defaults.h"

int first_step = 1;

VF initDistrib(NumMatrixTypes);
VF gapOpen(2 * NumInsertStates);
VF gapExtend(2 * NumInsertStates);
VVF emitPairs(256, VF(256, 1e-10));
VF emitSingle(256, 1e-5);

string alphabet = alphabetDefault;

string parametersOutputFilename = "tmp.train";

Train::Train(int argc, char* argv[])
{
    //parse program parameters
    SafeVector<string> sequenceNames = ParseParams(argc, argv);
    //initialize arguments for partition function
    ReadParameters();
    // build new model for aligning
    ProbabilisticModel model (initDistrib, gapOpen, gapExtend, emitPairs, emitSingle);

// prepare to average parameters
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] = 0;
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] = 0;
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] = 0;
    for (int i = 0; i < (int) emitPairs.size(); i++)
        for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] = 0;
    for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] = 0;

    // align each file individually
    for (int i = 0; i < (int) sequenceNames.size(); i++)
    {

        VF thisInitDistrib (NumMatrixTypes);
        VF thisGapOpen (2*NumInsertStates);
        VF thisGapExtend (2*NumInsertStates);
        VVF thisEmitPairs (256, VF (256, 1e-10));
        VF thisEmitSingle (256, 1e-5);
        MultiSequence *sequences = new MultiSequence();
        assert(sequences);
        sequences->LoadMFA (sequenceNames[i], true);
        // align sequences
        DoAlignTrain (sequences, model, thisInitDistrib, thisGapOpen, thisGapExtend, thisEmitPairs, thisEmitSingle, first_step);

        // add in contribution of the derived parameters
        for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] += thisInitDistrib[i];
        for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] += thisGapOpen[i];
        for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] += thisGapExtend[i];
        for (int i = 0; i < (int) emitPairs.size(); i++)
            for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] += thisEmitPairs[i][j];
        for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] += thisEmitSingle[i];

        delete sequences;
    }

    // compute new parameters and print them out
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] /= (int) sequenceNames.size();
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] /= (int) sequenceNames.size();
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] /= (int) sequenceNames.size();
    for (int i = 0; i < (int) emitPairs.size(); i++)
        for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] /= (int) sequenceNames.size();
    for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] /= sequenceNames.size();


    PrintParameters (initDistrib, gapOpen, gapExtend, emitPairs, emitSingle,
                     parametersOutputFilename.c_str());
}

Train::~Train()
{
}

/////////////////////////////////////////////////////////////////
// PrintParameters()
/////////////////////////////////////////////////////////////////

void Train::PrintParameters (const VF &initDistrib, const VF &gapOpen,
                             const VF &gapExtend, const VVF &emitPairs, const VF &emitSingle, const char *filename)
{

    // if a file name is specified
    if (filename)
    {
        // attempt to open the file for writing
        FILE *file = fopen (filename, "w");
        if (!file)
        {
            cerr << "ERROR: Unable to write parameter file: " << filename << endl;
            exit (1);
        }

        // if successful, then write the parameters to the file
        for (int i = 0; i < NumMatrixTypes; i++) fprintf (file, "%.10f ", initDistrib[i]);
        fprintf (file, "\n");
        for (int i = 0; i < 2*NumInsertStates; i++) fprintf (file, "%.10f ", gapOpen[i]);
        fprintf (file, "\n");
        for (int i = 0; i < 2*NumInsertStates; i++) fprintf (file, "%.10f ", gapExtend[i]);
        fprintf (file, "\n");
        fprintf (file, "%s\n", alphabet.c_str());
        for (int i = 0; i < (int) alphabet.size(); i++)
        {
            for (int j = 0; j <= i; j++)
                fprintf (file, "%.10f ", emitPairs[(unsigned char) alphabet[i]][(unsigned char) alphabet[j]]);
            fprintf (file, "\n");
        }
        for (int i = 0; i < (int) alphabet.size(); i++)
            fprintf (file, "%.10f ", emitSingle[(unsigned char) alphabet[i]]);
        fprintf (file, "\n");
        fclose (file);
    }
}


/////////////////////////////////////////////////////////////////
// doAlign()
//
// First computes all pairwise posterior probability matrices.
// Then, computes new parameters if training, or a final
// alignment, otherwise.
/////////////////////////////////////////////////////////////////
void Train::DoAlignTrain(MultiSequence *sequences,
                         const ProbabilisticModel &model, VF &initDistrib, VF &gapOpen,
                         VF &gapExtend, VVF &emitPairs, VF &emitSingle, int first_step)
{
    assert(sequences);
    //get the number of sequences
    const int numSeqs = sequences->GetNumSequences();

    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] = 0;
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] = 0;
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] = 0;
    for (int i = 0; i < (int) emitPairs.size(); i++)
        for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] = 0;
    for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] = 0;

    for (int a = 0; a < numSeqs - 1; a++)
    {
        for (int b = a + 1; b < numSeqs; b++)
        {
            Sequence *seq1 = sequences->GetSequence(a);
            Sequence *seq2 = sequences->GetSequence(b);
            //posterior probability matrix
            VF* forward ;
            VF* backward;
            //Double Affine pair-HMM

            if(first_step == 1)
            {
                // compute forward and backward probabilities
                forward = model.ComputeForwardMatrix(seq1, seq2);
                assert(forward);
                backward = model.ComputeBackwardMatrix(seq1, seq2);
                assert(backward);
            }
            else
            {
                //local pair-HMM
                // compute forward and backward probabilities
                forward = model.ComputeForwardMatrix(seq1, seq2,false);
                assert(forward);
                backward = model.ComputeBackwardMatrix(seq1, seq2,false);
                assert(backward);
            }
            // compute new parameters
            VF thisInitDistrib (NumMatrixTypes);
            VF thisGapOpen (2*NumInsertStates);
            VF thisGapExtend (2*NumInsertStates);
            VVF thisEmitPairs (256, VF (256, 1e-10));
            VF thisEmitSingle (256, 1e-5);

            model.ComputeNewParameters (seq1, seq2, *forward, *backward, thisInitDistrib, thisGapOpen, thisGapExtend, thisEmitPairs, thisEmitSingle, true);

            // add in contribution of the derived parameters
            for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] += thisInitDistrib[i];
            for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] += thisGapOpen[i];
            for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] += thisGapExtend[i];
            for (int i = 0; i < (int) emitPairs.size(); i++)
                for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] += thisEmitPairs[i][j];
            for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] += thisEmitSingle[i];
            delete forward;
            delete backward;
        }
    }
    // compute new parameters
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] /= numSeqs * (numSeqs - 1) / 2;
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] /= numSeqs * (numSeqs - 1) / 2;
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] /= numSeqs * (numSeqs - 1) / 2;
    for (int i = 0; i < (int) emitPairs.size(); i++)
        for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] /= numSeqs * (numSeqs - 1) / 2;
    for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] /= numSeqs * (numSeqs - 1) / 2;
}

/////////////////////////////////////////////////////////////////
// ReadParameters()
//
// Read initial distribution, transition, and emission
// parameters from a file.
/////////////////////////////////////////////////////////////////

void Train::ReadParameters()
{
    emitPairs = VVF(256, VF(256, 1e-10));
    emitSingle = VF(256, 1e-5);
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
// ParseParams()
//
// Parse all command-line options.
/////////////////////////////////////////////////////////////////

SafeVector<string> Train::ParseParams(int argc, char **argv)
{
    if (argc < 2)
    {
        exit(1);
    }
    SafeVector<string> sequenceNames;
    int tempInt;

    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (!strcmp(argv[i], "-p")
                    || !strcmp(argv[i], "-program"))
            {
                if (!GetInteger(argv[++i], &tempInt) && i < argc - 1)
                {
                    first_step = 0;
                }
                else
                {
                    if (tempInt >= 2 || tempInt < 0)
                    {
                        first_step = 0;
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
                    }
                }
            }

            else  if (!strcmp(argv[i], "-f")
                      || !strcmp(argv[i], "-file"))
            {
                parametersOutputFilename = string(argv[++i]);
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
    return sequenceNames;
}



int main(int argc, char* argv[])
{
    Train train(argc, argv);
    return 0;
}
