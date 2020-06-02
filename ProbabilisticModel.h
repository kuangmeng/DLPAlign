/////////////////////////////////////////////////////////////////
// ProbabilisticModel.h
//
// Routines for (1) local pair-HMM posterior probability computations
// 		(2) double affine pair-HMM posterior probability computations
//              (2) chained anchoring
//              (3) maximum weight trace alignment
/////////////////////////////////////////////////////////////////

#ifndef PROBABILISTICMODEL_H
#define PROBABILISTICMODEL_H

#include <list>
#include <cmath>
#include <cstdio>
#include "SafeVector.h"
#include "ScoreType.h"
#include "SparseMatrix.h"
#include "MultiSequence.h"


using namespace std;

const int NumMatchStates = 1;
const int NumInsertStates = 2;                                             // for double affine pair-HMM
const int NumMatrixTypes = NumMatchStates + NumInsertStates * 2;

/////////////////////////////////////////////////////////////////
// ProbabilisticModel
//
// Class for storing the parameters of a probabilistic model and
// performing different computations based on those parameters.
// In particular, this class handles the computation of
// posterior probabilities that may be used in alignment.
/////////////////////////////////////////////////////////////////

class ProbabilisticModel
{

    float initialDistribution[NumMatrixTypes];               // holds the initial probabilities for each state
    float transProb[NumMatrixTypes][NumMatrixTypes];         // holds all state-to-state transition probabilities for double affine pair-HMM
    float matchProb[256][256];                               // emission probabilities for match states
    float insProb[256][NumMatrixTypes];                      // emission probabilities for insert states
    float local_transProb[3][3];				   // holds central state-to-state transition probabilities for local pair-HMM
    float random_transProb[2];				   // holds flanking state-to-state transition probabilities for local pair-HMM

public:

    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::ProbabilisticModel()
    //
    // Constructor.  Builds a new probabilistic model using the
    // given parameters.
    /////////////////////////////////////////////////////////////////

    ProbabilisticModel (const VF &initDistribMat, const VF &gapOpen, const VF &gapExtend,
                        const VVF &emitPairs, const VF &emitSingle)
    {

        /**********************************
        double affine pair-HMM:
        initDistribMat[0,1,3,4]: initialization of starting states
        gapOpen[0,2]: long and short open gap probabilities
        gapExtend[0,2]:long and short extend gap probabilities

        local pair-HMM:
        initDistribMat[2]: leave from a flanking state probability
        gapOpen[1]: central open gap probability
        gapExtend[1]:central extend gap probability
        ***********************************/

//double affine pair-HMM
        // build transition matrix
        VVF transMat (NumMatrixTypes, VF (NumMatrixTypes, 0.0f));
        transMat[0][0] = 1;
        for (int i = 0; i < NumInsertStates; i++)
        {
            transMat[0][2*i+1] = gapOpen[2*i];
            transMat[0][2*i+2] = gapOpen[2*i];
            transMat[0][0] -= (gapOpen[2*i] + gapOpen[2*i]);
            assert (transMat[0][0] > 0);
            transMat[2*i+1][2*i+1] = gapExtend[2*i];
            transMat[2*i+2][2*i+2] = gapExtend[2*i];
            transMat[2*i+1][2*i+2] = 0;
            transMat[2*i+2][2*i+1] = 0;
            transMat[2*i+1][0] = 1 - gapExtend[2*i];
            transMat[2*i+2][0] = 1 - gapExtend[2*i];
        }

        // create initial and transition probability matrices
        for (int i = 0; i < NumMatrixTypes; i++)
        {
            initialDistribution[i] = LOG (initDistribMat[i]);
            for (int j = 0; j < NumMatrixTypes; j++)
                transProb[i][j] = LOG (transMat[i][j]);
        }

        //due to local model parameters' initilization
        //need to correct initialDistribution[2]
        initialDistribution[2] = LOG (initDistribMat[1]);

        // create insertion and match probability matrices
        for (int i = 0; i < 256; i++)
        {
            for (int j = 0; j < NumMatrixTypes; j++)
                insProb[i][j] = LOG (emitSingle[i]);
            for (int j = 0; j < 256; j++)
                matchProb[i][j] = LOG (emitPairs[i][j]);
        }

//local pair-HMM
        // build transition matrix
        VVF ltransMat (3, VF (3, 0.0f));
        ltransMat[0][0] = 1;

        ltransMat[0][1] = gapOpen[1];
        ltransMat[0][2] = gapOpen[1];
        ltransMat[0][0] -= (gapOpen[1] + gapOpen[1]);
        assert (ltransMat[0][0] > 0);
        ltransMat[1][1] = gapExtend[1];
        ltransMat[2][2] = gapExtend[1];
        ltransMat[1][2] = 0;
        ltransMat[2][1] = 0;
        ltransMat[1][0] = 1 - gapExtend[1];
        ltransMat[2][0] = 1 - gapExtend[1];

        // create initial and transition probability matrices
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
                local_transProb[i][j] = LOG (ltransMat[i][j]);
        }

        // create initial and transition probability matrices
        random_transProb[0] = LOG (initDistribMat[2]);//probability to leave from a randam state
        random_transProb[1] = LOG (1-initDistribMat[2]);//probability to stay in a randam state

    }

    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::ComputeForwardMatrix()
    //
    // Computes a set of forward probability matrices for aligning
    // seq1 and seq2.
    //
    // For efficiency reasons, a single-dimensional floating-point
    // array is used here, with the following indexing scheme:
    //
    //    forward[i + NumMatrixTypes * (j * (seq2Length+1) + k)]
    //    refers to the probability of aligning through j characters
    //    of the first sequence, k characters of the second sequence,
    //    and ending in state i.
    //    flag: 1 probcons, 0 local
    /////////////////////////////////////////////////////////////////

    VF *ComputeForwardMatrix (Sequence *seq1, Sequence *seq2, bool flag=true) const
    {

        assert (seq1);
        assert (seq2);

        const int seq1Length = seq1->GetLength();
        const int seq2Length = seq2->GetLength();

        // retrieve the points to the beginning of each sequence
        SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
        SafeVector<char>::iterator iter2 = seq2->GetDataPtr();

        // create matrix
        VF *forwardPtr;
        if(flag) forwardPtr = new VF (NumMatrixTypes * (seq1Length+1) * (seq2Length+1), LOG_ZERO);
        else forwardPtr = new VF (3 * (seq1Length+1) * (seq2Length+1), LOG_ZERO);
        assert (forwardPtr);
        VF &forward = *forwardPtr;

        // initialization condition
        if(flag)
        {
            forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] =
                initialDistribution[0] + matchProb[(unsigned char) iter1[1]][(unsigned char) iter2[1]];

            for (int k = 0; k < NumInsertStates; k++)
            {
                forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] =
                    initialDistribution[2*k+1] + insProb[(unsigned char) iter1[1]][k];
                forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] =
                    initialDistribution[2*k+2] + insProb[(unsigned char) iter2[1]][k];
            }
        }

        // remember offset for each index combination
        int ij = 0;
        int i1j = -seq2Length - 1;
        int ij1 = -1;
        int i1j1 = -seq2Length - 2;

        if(flag)
        {
            ij *= NumMatrixTypes;
            i1j *= NumMatrixTypes;
            ij1 *= NumMatrixTypes;
            i1j1 *= NumMatrixTypes;
        }
        else
        {
            ij *= 3;
            i1j *= 3;
            ij1 *= 3;
            i1j1 *= 3;
        }

        // compute forward scores
        for (int i = 0; i <= seq1Length; i++)
        {
            unsigned char c1 = (i == 0) ? '~' : (unsigned char) iter1[i];
            for (int j = 0; j <= seq2Length; j++)
            {
                unsigned char c2 = (j == 0) ? '~' : (unsigned char) iter2[j];
                //local
                if(i == 1 && j == 1 && !flag) forward[0 + ij] =
                        matchProb[c1][c2] - insProb[c1][0] - insProb[c2][0] - 2*random_transProb[1];

                if (i > 1 || j > 1)
                {
                    if (i > 0 && j > 0)
                    {
                        if(flag)
                        {
                            forward[0 + ij] = forward[0 + i1j1] + transProb[0][0];
                            for (int k = 1; k < NumMatrixTypes; k++)
                                LOG_PLUS_EQUALS (forward[0 + ij], forward[k + i1j1] + transProb[k][0]);
                            forward[0 + ij] += matchProb[c1][c2];
                        }
                        //local
                        else
                        {
                            forward[0 + ij] = matchProb[c1][c2] - insProb[c1][0] - insProb[c2][0] - 2*random_transProb[1];
                            for (int k = 0; k < 3; k++)
                                LOG_PLUS_EQUALS (forward[0 + ij], matchProb[c1][c2] - insProb[c1][0] - insProb[c2][0] +
                                                 forward[k + i1j1] + local_transProb[k][0] - 2*random_transProb[1]);
                        }
                    }
                    if (i > 0)
                    {
                        if(flag)
                        {
                            for (int k = 0; k < NumInsertStates; k++)
                                forward[2*k+1 + ij] = insProb[c1][k] +
                                                      LOG_ADD (forward[0 + i1j] + transProb[0][2*k+1],
                                                               forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1]);
                        }

                        //local
                        else
                        {
                            forward[1 + ij] = LOG_ADD (forward[0 + i1j] + local_transProb[0][1] - random_transProb[1],
                                                       forward[1 + i1j] + local_transProb[1][1] - random_transProb[1]);
                        }

                    }
                    if (j > 0)
                    {
                        if(flag)
                        {
                            for (int k = 0; k < NumInsertStates; k++)
                                forward[2*k+2 + ij] = insProb[c2][k] +
                                                      LOG_ADD (forward[0 + ij1] + transProb[0][2*k+2],
                                                               forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2]);
                        }
                        //local
                        else
                        {
                            forward[2 + ij] = LOG_ADD (forward[0 + ij1] + local_transProb[0][2] - random_transProb[1],
                                                       forward[2 + ij1] + local_transProb[2][2] - random_transProb[1]);
                        }
                    }
                }
                if(flag)
                {
                    ij += NumMatrixTypes;
                    i1j += NumMatrixTypes;
                    ij1 += NumMatrixTypes;
                    i1j1 += NumMatrixTypes;
                }
                else
                {
                    ij += 3;
                    i1j += 3;
                    ij1 += 3;
                    i1j1 += 3;
                }
            }
        }

        return forwardPtr;
    }

    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::ComputeBackwardMatrix()
    //
    // Computes a set of backward probability matrices for aligning
    // seq1 and seq2.
    //
    // For efficiency reasons, a single-dimensional floating-point
    // array is used here, with the following indexing scheme:
    //
    //    backward[i + NumMatrixTypes * (j * (seq2Length+1) + k)]
    //    refers to the probability of starting in state i and
    //    aligning from character j+1 to the end of the first
    //    sequence and from character k+1 to the end of the second
    //    sequence.
    /////////////////////////////////////////////////////////////////

    VF *ComputeBackwardMatrix (Sequence *seq1, Sequence *seq2, bool flag=true) const
    {

        assert (seq1);
        assert (seq2);

        const int seq1Length = seq1->GetLength();
        const int seq2Length = seq2->GetLength();
        SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
        SafeVector<char>::iterator iter2 = seq2->GetDataPtr();

        // create matrix
        VF *backwardPtr;
        if(flag) backwardPtr = new VF (NumMatrixTypes * (seq1Length+1) * (seq2Length+1), LOG_ZERO);
        else backwardPtr = new VF (3 * (seq1Length+1) * (seq2Length+1), LOG_ZERO);
        assert (backwardPtr);
        VF &backward = *backwardPtr;

        // initialization condition
        if(flag)
        {
            for (int k = 0; k < NumMatrixTypes; k++)
                backward[NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1) + k] = initialDistribution[k];
        }
        // remember offset for each index combination
        int ij = (seq1Length+1) * (seq2Length+1) - 1;
        int i1j = ij + seq2Length + 1;
        int ij1 = ij + 1;
        int i1j1 = ij + seq2Length + 2;

        if(flag)
        {
            ij *= NumMatrixTypes;
            i1j *= NumMatrixTypes;
            ij1 *= NumMatrixTypes;
            i1j1 *= NumMatrixTypes;
        }
        else
        {
            ij *= 3;
            i1j *= 3;
            ij1 *= 3;
            i1j1 *= 3;
        }

        // compute backward scores
        for (int i = seq1Length; i >= 0; i--)
        {
            unsigned char c1 = (i == seq1Length) ? '~' : (unsigned char) iter1[i+1];
            for (int j = seq2Length; j >= 0; j--)
            {
                unsigned char c2 = (j == seq2Length) ? '~' : (unsigned char) iter2[j+1];

                if(!flag) backward[0 + ij] = LOG_ONE;//local
                if (i < seq1Length && j < seq2Length)
                {
                    if(flag)
                    {
                        const float ProbXY = backward[0 + i1j1] + matchProb[c1][c2];
                        for (int k = 0; k < NumMatrixTypes; k++)
                            LOG_PLUS_EQUALS (backward[k + ij], ProbXY + transProb[k][0]);
                    }
                    //local
                    else
                    {
                        const float ProbXY = backward[0 + i1j1] + matchProb[c1][c2] - insProb[c1][0] - insProb[c2][0];
                        for (int k = 0; k < 3; k++)
                            LOG_PLUS_EQUALS (backward[k + ij], ProbXY + local_transProb[k][0] - 2*random_transProb[1] );
                    }
                }
                if (i < seq1Length)
                {
                    if(flag)
                    {
                        for (int k = 0; k < NumInsertStates; k++)
                        {
                            LOG_PLUS_EQUALS (backward[0 + ij], backward[2*k+1 + i1j] + insProb[c1][k] + transProb[0][2*k+1]);
                            LOG_PLUS_EQUALS (backward[2*k+1 + ij], backward[2*k+1 + i1j] + insProb[c1][k] + transProb[2*k+1][2*k+1]);
                        }
                    }
                    //local
                    else
                    {
                        LOG_PLUS_EQUALS (backward[0 + ij], backward[1 + i1j] + local_transProb[0][1] - random_transProb[1]);
                        LOG_PLUS_EQUALS (backward[1 + ij], backward[1 + i1j] + local_transProb[1][1] - random_transProb[1]);
                    }
                }
                if (j < seq2Length)
                {
                    if(flag)
                    {
                        for (int k = 0; k < NumInsertStates; k++)
                        {
                            LOG_PLUS_EQUALS (backward[0 + ij], backward[2*k+2 + ij1] + insProb[c2][k] + transProb[0][2*k+2]);
                            LOG_PLUS_EQUALS (backward[2*k+2 + ij], backward[2*k+2 + ij1] + insProb[c2][k] + transProb[2*k+2][2*k+2]);
                        }
                    }
                    //local
                    else
                    {
                        LOG_PLUS_EQUALS (backward[0 + ij], backward[2 + ij1] + local_transProb[0][2] - random_transProb[1]);
                        LOG_PLUS_EQUALS (backward[2 + ij], backward[2 + ij1] + local_transProb[2][2] - random_transProb[1]);
                    }
                }
                if(flag)
                {
                    ij -= NumMatrixTypes;
                    i1j -= NumMatrixTypes;
                    ij1 -= NumMatrixTypes;
                    i1j1 -= NumMatrixTypes;
                }
                else
                {
                    ij -= 3;
                    i1j -= 3;
                    ij1 -= 3;
                    i1j1 -= 3;
                }
            }
        }

        return backwardPtr;
    }

    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::ComputeTotalProbability()
    //
    // Computes the total probability of an alignment given
    // the forward and backward matrices.
    // flag: 1 probcons, 0 local
    /////////////////////////////////////////////////////////////////

    float ComputeTotalProbability (Sequence *seq1, Sequence *seq2,
                                   const VF &forward, const VF &backward, bool flag=true) const
    {

        // compute total probability
        float totalForwardProb = LOG_ZERO;
        float totalBackwardProb = LOG_ZERO;
        const int seq1Length = seq1->GetLength();
        const int seq2Length = seq2->GetLength();

        if(flag)
        {
            for (int k = 0; k < NumMatrixTypes; k++)
            {
                LOG_PLUS_EQUALS (totalForwardProb,
                                 forward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] +
                                 backward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
            }

            totalBackwardProb =
                forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] +
                backward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)];

            for (int k = 0; k < NumInsertStates; k++)
            {
                LOG_PLUS_EQUALS (totalBackwardProb,
                                 forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] +
                                 backward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)]);
                LOG_PLUS_EQUALS (totalBackwardProb,
                                 forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] +
                                 backward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)]);
            }
        }
        else
        {
            SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
            SafeVector<char>::iterator iter2 = seq2->GetDataPtr();
            int ij = 0;
            for (int i = 0; i <= seq1Length; i++)
            {
                unsigned char c1 = (i == 0) ? '~' : (unsigned char) iter1[i];
                for (int j = 0; j <= seq2Length; j++)
                {
                    unsigned char c2 = (j == 0) ? '~' : (unsigned char) iter2[j];
                    if(i>0&&j>0)
                    {
                        LOG_PLUS_EQUALS (totalForwardProb,forward[ij]);
                        LOG_PLUS_EQUALS (totalBackwardProb,backward[ij] + matchProb[c1][c2]
                                         - insProb[c1][0] - insProb[c2][0] - 2*random_transProb[1]);
                    }
                    ij += 3;
                }
            }

        }

        return (totalForwardProb + totalBackwardProb) / 2;
    }

    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::ComputePosteriorMatrix()
    //
    // Computes the posterior probability matrix based on
    // the forward and backward matrices.
    // flag: 1 probcons, 0 local
    /////////////////////////////////////////////////////////////////

    VF *ComputePosteriorMatrix (Sequence *seq1, Sequence *seq2,
                                const VF &forward, const VF &backward, bool flag=true) const
    {

        assert (seq1);
        assert (seq2);

        const int seq1Length = seq1->GetLength();
        const int seq2Length = seq2->GetLength();

        float totalProb = ComputeTotalProbability (seq1, seq2,forward, backward, flag);

        // compute posterior matrices
        VF *posteriorPtr = new VF((seq1Length+1) * (seq2Length+1));
        assert (posteriorPtr);
        VF &posterior = *posteriorPtr;

        int ij = 0;
        VF::iterator ptr = posterior.begin();

        for (int i = 0; i <= seq1Length; i++)
        {
            for (int j = 0; j <= seq2Length; j++)
            {
                *(ptr++) = EXP (min (LOG_ONE, forward[ij] + backward[ij] - totalProb));
                if(flag) ij += NumMatrixTypes;
                else ij += 3;
            }
        }

        posterior[0] = 0;

        return posteriorPtr;
    }

    /*
    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::ComputeExpectedCounts()
    //
    // Computes the expected counts for the various transitions.
    /////////////////////////////////////////////////////////////////

    VVF *ComputeExpectedCounts () const {

      assert (seq1);
      assert (seq2);

      const int seq1Length = seq1->GetLength();
      const int seq2Length = seq2->GetLength();
      SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
      SafeVector<char>::iterator iter2 = seq2->GetDataPtr();

      // compute total probability
      float totalProb = ComputeTotalProbability (seq1Length, seq2Length,
                                                 forward, backward);

      // initialize expected counts
      VVF *countsPtr = new VVF(NumMatrixTypes + 1, VF(NumMatrixTypes, LOG_ZERO)); assert (countsPtr);
      VVF &counts = *countsPtr;

      // remember offset for each index combination
      int ij = 0;
      int i1j = -seq2Length - 1;
      int ij1 = -1;
      int i1j1 = -seq2Length - 2;

      ij *= NumMatrixTypes;
      i1j *= NumMatrixTypes;
      ij1 *= NumMatrixTypes;
      i1j1 *= NumMatrixTypes;

      // compute expected counts
      for (int i = 0; i <= seq1Length; i++){
        unsigned char c1 = (i == 0) ? '~' : (unsigned char) iter1[i];
        for (int j = 0; j <= seq2Length; j++){
          unsigned char c2 = (j == 0) ? '~' : (unsigned char) iter2[j];

          if (i > 0 && j > 0){
            for (int k = 0; k < NumMatrixTypes; k++)
              LOG_PLUS_EQUALS (counts[k][0],
                               forward[k + i1j1] + transProb[k][0] +
                               matchProb[c1][c2] + backward[0 + ij]);
          }
          if (i > 0){
            for (int k = 0; k < NumInsertStates; k++){
              LOG_PLUS_EQUALS (counts[0][2*k+1],
                               forward[0 + i1j] + transProb[0][2*k+1] +
                               insProb[c1][k] + backward[2*k+1 + ij]);
              LOG_PLUS_EQUALS (counts[2*k+1][2*k+1],
                               forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1] +
                               insProb[c1][k] + backward[2*k+1 + ij]);
            }
          }
          if (j > 0){
            for (int k = 0; k < NumInsertStates; k++){
              LOG_PLUS_EQUALS (counts[0][2*k+2],
                               forward[0 + ij1] + transProb[0][2*k+2] +
                               insProb[c2][k] + backward[2*k+2 + ij]);
              LOG_PLUS_EQUALS (counts[2*k+2][2*k+2],
                               forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2] +
                               insProb[c2][k] + backward[2*k+2 + ij]);
            }
          }

          ij += NumMatrixTypes;
          i1j += NumMatrixTypes;
          ij1 += NumMatrixTypes;
          i1j1 += NumMatrixTypes;
        }
      }

      // scale all expected counts appropriately
      for (int i = 0; i < NumMatrixTypes; i++)
        for (int j = 0; j < NumMatrixTypes; j++)
          counts[i][j] -= totalProb;

    }
    */

    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::ComputeNewParameters()
    //
    // Computes a new parameter set based on the expected counts
    // given.
    /////////////////////////////////////////////////////////////////

    void ComputeNewParameters (Sequence *seq1, Sequence *seq2,
                               const VF &forward, const VF &backward,
                               VF &initDistribMat, VF &gapOpen,
                               VF &gapExtend, VVF &emitPairs, VF &emitSingle, bool enableTrainEmissions) const
    {

        assert (seq1);
        assert (seq2);

        const int seq1Length = seq1->GetLength();
        const int seq2Length = seq2->GetLength();
        SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
        SafeVector<char>::iterator iter2 = seq2->GetDataPtr();

        // compute total probability
        float totalProb = ComputeTotalProbability (seq1, seq2,
                          forward, backward);

        // initialize expected counts
        VVF transCounts (NumMatrixTypes, VF (NumMatrixTypes, LOG_ZERO));
        VF initCounts (NumMatrixTypes, LOG_ZERO);
        VVF pairCounts (256, VF (256, LOG_ZERO));
        VF singleCounts (256, LOG_ZERO);

        // remember offset for each index combination
        int ij = 0;
        int i1j = -seq2Length - 1;
        int ij1 = -1;
        int i1j1 = -seq2Length - 2;

        ij *= NumMatrixTypes;
        i1j *= NumMatrixTypes;
        ij1 *= NumMatrixTypes;
        i1j1 *= NumMatrixTypes;

        // compute initial distribution posteriors
        initCounts[0] = LOG_ADD (forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] +
                                 backward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)],
                                 forward[0 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] +
                                 backward[0 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
        for (int k = 0; k < NumInsertStates; k++)
        {
            initCounts[2*k+1] = LOG_ADD (forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] +
                                         backward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)],
                                         forward[2*k+1 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] +
                                         backward[2*k+1 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
            initCounts[2*k+2] = LOG_ADD (forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] +
                                         backward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)],
                                         forward[2*k+2 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] +
                                         backward[2*k+2 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
        }

        // compute expected counts
        for (int i = 0; i <= seq1Length; i++)
        {
            unsigned char c1 = (i == 0) ? '~' : (unsigned char) toupper(iter1[i]);
            for (int j = 0; j <= seq2Length; j++)
            {
                unsigned char c2 = (j == 0) ? '~' : (unsigned char) toupper(iter2[j]);

                if (i > 0 && j > 0)
                {
                    if (enableTrainEmissions && i == 1 && j == 1)
                    {
                        LOG_PLUS_EQUALS (pairCounts[c1][c2],
                                         initialDistribution[0] + matchProb[c1][c2] + backward[0 + ij]);
                        LOG_PLUS_EQUALS (pairCounts[c2][c1],
                                         initialDistribution[0] + matchProb[c2][c1] + backward[0 + ij]);
                    }

                    for (int k = 0; k < NumMatrixTypes; k++)
                    {
                        LOG_PLUS_EQUALS (transCounts[k][0],
                                         forward[k + i1j1] + transProb[k][0] +
                                         matchProb[c1][c2] + backward[0 + ij]);
                        if (enableTrainEmissions && (i != 1 || j != 1))
                        {
                            LOG_PLUS_EQUALS (pairCounts[c1][c2],
                                             forward[k + i1j1] + transProb[k][0] +
                                             matchProb[c1][c2] + backward[0 + ij]);
                            LOG_PLUS_EQUALS (pairCounts[c2][c1],
                                             forward[k + i1j1] + transProb[k][0] +
                                             matchProb[c2][c1] + backward[0 + ij]);
                        }

                    }
                }
                if (i > 0)
                {
                    for (int k = 0; k < NumInsertStates; k++)
                    {
                        LOG_PLUS_EQUALS (transCounts[0][2*k+1],
                                         forward[0 + i1j] + transProb[0][2*k+1] +
                                         insProb[c1][k] + backward[2*k+1 + ij]);
                        LOG_PLUS_EQUALS (transCounts[2*k+1][2*k+1],
                                         forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1] +
                                         insProb[c1][k] + backward[2*k+1 + ij]);
                        if (enableTrainEmissions)
                        {
                            if (i == 1 && j == 0)
                            {
                                LOG_PLUS_EQUALS (singleCounts[c1],
                                                 initialDistribution[2*k+1] + insProb[c1][k] + backward[2*k+1 + ij]);
                            }
                            else
                            {
                                LOG_PLUS_EQUALS (singleCounts[c1],
                                                 forward[0 + i1j] + transProb[0][2*k+1] +
                                                 insProb[c1][k] + backward[2*k+1 + ij]);
                                LOG_PLUS_EQUALS (singleCounts[c1],
                                                 forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1] +
                                                 insProb[c1][k] + backward[2*k+1 + ij]);
                            }
                        }
                    }
                }
                if (j > 0)
                {
                    for (int k = 0; k < NumInsertStates; k++)
                    {
                        LOG_PLUS_EQUALS (transCounts[0][2*k+2],
                                         forward[0 + ij1] + transProb[0][2*k+2] +
                                         insProb[c2][k] + backward[2*k+2 + ij]);
                        LOG_PLUS_EQUALS (transCounts[2*k+2][2*k+2],
                                         forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2] +
                                         insProb[c2][k] + backward[2*k+2 + ij]);
                        if (enableTrainEmissions)
                        {
                            if (i == 0 && j == 1)
                            {
                                LOG_PLUS_EQUALS (singleCounts[c2],
                                                 initialDistribution[2*k+2] + insProb[c2][k] + backward[2*k+2 + ij]);
                            }
                            else
                            {
                                LOG_PLUS_EQUALS (singleCounts[c2],
                                                 forward[0 + ij1] + transProb[0][2*k+2] +
                                                 insProb[c2][k] + backward[2*k+2 + ij]);
                                LOG_PLUS_EQUALS (singleCounts[c2],
                                                 forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2] +
                                                 insProb[c2][k] + backward[2*k+2 + ij]);
                            }
                        }
                    }
                }

                ij += NumMatrixTypes;
                i1j += NumMatrixTypes;
                ij1 += NumMatrixTypes;
                i1j1 += NumMatrixTypes;
            }
        }

        // scale all expected counts appropriately
        for (int i = 0; i < NumMatrixTypes; i++)
        {
            initCounts[i] -= totalProb;
            for (int j = 0; j < NumMatrixTypes; j++)
                transCounts[i][j] -= totalProb;
        }
        if (enableTrainEmissions)
        {
            for (int i = 0; i < 256; i++)
            {
                for (int j = 0; j < 256; j++)
                    pairCounts[i][j] -= totalProb;
                singleCounts[i] -= totalProb;
            }
        }

        // compute new initial distribution
        float totalInitDistribCounts = 0;
        for (int i = 0; i < NumMatrixTypes; i++)
            totalInitDistribCounts += exp (initCounts[i]); // should be 2
        initDistribMat[0] = min (1.0f, max (0.0f, (float) exp (initCounts[0]) / totalInitDistribCounts));
        for (int k = 0; k < NumInsertStates; k++)
        {
            float val = (exp (initCounts[2*k+1]) + exp (initCounts[2*k+2])) / 2;
            initDistribMat[2*k+1] = initDistribMat[2*k+2] = min (1.0f, max (0.0f, val / totalInitDistribCounts));
        }

        // compute total counts for match state
        float inMatchStateCounts = 0;
        for (int i = 0; i < NumMatrixTypes; i++)
            inMatchStateCounts += exp (transCounts[0][i]);
        for (int i = 0; i < NumInsertStates; i++)
        {

            // compute total counts for gap state
            float inGapStateCounts =
                exp (transCounts[2*i+1][0]) +
                exp (transCounts[2*i+1][2*i+1]) +
                exp (transCounts[2*i+2][0]) +
                exp (transCounts[2*i+2][2*i+2]);

            gapOpen[2*i] = gapOpen[2*i+1] =
                               (exp (transCounts[0][2*i+1]) +
                                exp (transCounts[0][2*i+2])) /
                               (2 * inMatchStateCounts);

            gapExtend[2*i] = gapExtend[2*i+1] =
                                 (exp (transCounts[2*i+1][2*i+1]) +
                                  exp (transCounts[2*i+2][2*i+2])) /
                                 inGapStateCounts;
        }

        if (enableTrainEmissions)
        {
            float totalPairCounts = 0;
            float totalSingleCounts = 0;
            for (int i = 0; i < 256; i++)
            {
                for (int j = 0; j <= i; j++)
                    totalPairCounts += exp (pairCounts[j][i]);
                totalSingleCounts += exp (singleCounts[i]);
            }

            for (int i = 0; i < 256; i++) if (!islower ((char) i))
                {
                    int li = (int)((unsigned char) tolower ((char) i));
                    for (int j = 0; j <= i; j++) if (!islower ((char) j))
                        {
                            int lj = (int)((unsigned char) tolower ((char) j));
                            emitPairs[i][j] = emitPairs[i][lj] = emitPairs[li][j] = emitPairs[li][lj] =
                                    emitPairs[j][i] = emitPairs[j][li] = emitPairs[lj][i] = emitPairs[lj][li] = exp(pairCounts[j][i]) / totalPairCounts;
                        }
                    emitSingle[i] = emitSingle[li] = exp(singleCounts[i]) / totalSingleCounts;
                }
        }
    }

    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::ComputeAlignment()
    //
    // Computes an alignment based on the given posterior matrix.
    // This is done by finding the maximum summing path (or
    // maximum weight trace) through the posterior matrix.  The
    // final alignment is returned as a pair consisting of:
    //    (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and
    //        denote insertions in one of the two sequences and
    //        B's denote that both sequences are present (i.e.
    //        matches).
    //    (2) a float indicating the sum achieved
    /////////////////////////////////////////////////////////////////

    pair<SafeVector<char> *, float> ComputeAlignment (int seq1Length, int seq2Length,
            const VF &posterior) const
    {

        float *twoRows = new float[(seq2Length+1)*2];
        assert (twoRows);
        float *oldRow = twoRows;
        float *newRow = twoRows + seq2Length + 1;

        char *tracebackMatrix = new char[(seq1Length+1)*(seq2Length+1)];
        assert (tracebackMatrix);
        char *tracebackPtr = tracebackMatrix;

        VF::const_iterator posteriorPtr = posterior.begin() + seq2Length + 1;

        // initialization
        for (int i = 0; i <= seq2Length; i++)
        {
            oldRow[i] = 0;
            *(tracebackPtr++) = 'L';
        }

        // fill in matrix
        for (int i = 1; i <= seq1Length; i++)
        {

            // initialize left column
            newRow[0] = 0;
            posteriorPtr++;
            *(tracebackPtr++) = 'U';

            // fill in rest of row
            for (int j = 1; j <= seq2Length; j++)
            {
                ChooseBestOfThree (*(posteriorPtr++) + oldRow[j-1], newRow[j-1], oldRow[j],
                                   'D', 'L', 'U', &newRow[j], tracebackPtr++);
            }

            // swap rows
            float *temp = oldRow;
            oldRow = newRow;
            newRow = temp;
        }

        // store best score
        float total = oldRow[seq2Length];
        delete [] twoRows;

        // compute traceback
        SafeVector<char> *alignment = new SafeVector<char>;
        assert (alignment);
        int r = seq1Length, c = seq2Length;
        while (r != 0 || c != 0)
        {
            char ch = tracebackMatrix[r*(seq2Length+1) + c];
            switch (ch)
            {
            case 'L':
                c--;
                alignment->push_back ('Y');
                break;
            case 'U':
                r--;
                alignment->push_back ('X');
                break;
            case 'D':
                c--;
                r--;
                alignment->push_back ('B');
                break;
            default:
                assert (false);
            }
        }

        delete [] tracebackMatrix;

        reverse (alignment->begin(), alignment->end());

        return make_pair(alignment, total);
    }

    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::ComputeAlignmentWithGapPenalties()
    //
    // Similar to ComputeAlignment() except with gap penalties.
    /////////////////////////////////////////////////////////////////

    pair<SafeVector<char> *, float> ComputeAlignmentWithGapPenalties (MultiSequence *align1,
            MultiSequence *align2,
            const VF &posterior, int numSeqs1,
            int numSeqs2,
            float gapOpenPenalty,
            float gapContinuePenalty) const
    {
        int seq1Length = align1->GetSequence(0)->GetLength();
        int seq2Length = align2->GetSequence(0)->GetLength();
        SafeVector<SafeVector<char>::iterator > dataPtrs1 (align1->GetNumSequences());
        SafeVector<SafeVector<char>::iterator > dataPtrs2 (align2->GetNumSequences());

        // grab character data
        for (int i = 0; i < align1->GetNumSequences(); i++)
            dataPtrs1[i] = align1->GetSequence(i)->GetDataPtr();
        for (int i = 0; i < align2->GetNumSequences(); i++)
            dataPtrs2[i] = align2->GetSequence(i)->GetDataPtr();

        // the number of active sequences at any given column is defined to be the
        // number of non-gap characters in that column; the number of gap opens at
        // any given column is defined to be the number of gap characters in that
        // column where the previous character in the respective sequence was not
        // a gap
        SafeVector<int> numActive1 (seq1Length+1), numGapOpens1 (seq1Length+1);
        SafeVector<int> numActive2 (seq2Length+1), numGapOpens2 (seq2Length+1);

        // compute number of active sequences and gap opens for each group
        for (int i = 0; i < align1->GetNumSequences(); i++)
        {
            SafeVector<char>::iterator dataPtr = align1->GetSequence(i)->GetDataPtr();
            numActive1[0] = numGapOpens1[0] = 0;
            for (int j = 1; j <= seq1Length; j++)
            {
                if (dataPtr[j] != '-')
                {
                    numActive1[j]++;
                    numGapOpens1[j] += (j != 1 && dataPtr[j-1] != '-');
                }
            }
        }
        for (int i = 0; i < align2->GetNumSequences(); i++)
        {
            SafeVector<char>::iterator dataPtr = align2->GetSequence(i)->GetDataPtr();
            numActive2[0] = numGapOpens2[0] = 0;
            for (int j = 1; j <= seq2Length; j++)
            {
                if (dataPtr[j] != '-')
                {
                    numActive2[j]++;
                    numGapOpens2[j] += (j != 1 && dataPtr[j-1] != '-');
                }
            }
        }

        VVF openingPenalty1 (numSeqs1+1, VF (numSeqs2+1));
        VF continuingPenalty1 (numSeqs1+1);
        VVF openingPenalty2 (numSeqs1+1, VF (numSeqs2+1));
        VF continuingPenalty2 (numSeqs2+1);

        // precompute penalties
        for (int i = 0; i <= numSeqs1; i++)
            for (int j = 0; j <= numSeqs2; j++)
                openingPenalty1[i][j] = i * (gapOpenPenalty * j + gapContinuePenalty * (numSeqs2 - j));
        for (int i = 0; i <= numSeqs1; i++)
            continuingPenalty1[i] = i * gapContinuePenalty * numSeqs2;
        for (int i = 0; i <= numSeqs2; i++)
            for (int j = 0; j <= numSeqs1; j++)
                openingPenalty2[i][j] = i * (gapOpenPenalty * j + gapContinuePenalty * (numSeqs1 - j));
        for (int i = 0; i <= numSeqs2; i++)
            continuingPenalty2[i] = i * gapContinuePenalty * numSeqs1;

        float *twoRows = new float[6*(seq2Length+1)];
        assert (twoRows);
        float *oldRowMatch = twoRows;
        float *newRowMatch = twoRows + (seq2Length+1);
        float *oldRowInsertX = twoRows + 2*(seq2Length+1);
        float *newRowInsertX = twoRows + 3*(seq2Length+1);
        float *oldRowInsertY = twoRows + 4*(seq2Length+1);
        float *newRowInsertY = twoRows + 5*(seq2Length+1);

        char *tracebackMatrix = new char[3*(seq1Length+1)*(seq2Length+1)];
        assert (tracebackMatrix);
        char *tracebackPtr = tracebackMatrix;

        VF::const_iterator posteriorPtr = posterior.begin() + seq2Length + 1;

        // initialization
        for (int i = 0; i <= seq2Length; i++)
        {
            oldRowMatch[i] = oldRowInsertX[i] = (i == 0) ? 0 : LOG_ZERO;
            oldRowInsertY[i] = (i == 0) ? 0 : oldRowInsertY[i-1] + continuingPenalty2[numActive2[i]];
            *(tracebackPtr) = *(tracebackPtr+1) = *(tracebackPtr+2) = 'Y';
            tracebackPtr += 3;
        }

        // fill in matrix
        for (int i = 1; i <= seq1Length; i++)
        {

            // initialize left column
            newRowMatch[0] = newRowInsertY[0] = LOG_ZERO;
            newRowInsertX[0] = oldRowInsertX[0] + continuingPenalty1[numActive1[i]];
            posteriorPtr++;
            *(tracebackPtr) = *(tracebackPtr+1) = *(tracebackPtr+2) = 'X';
            tracebackPtr += 3;

            // fill in rest of row
            for (int j = 1; j <= seq2Length; j++)
            {

                // going to MATCH state
                ChooseBestOfThree (oldRowMatch[j-1],
                                   oldRowInsertX[j-1],
                                   oldRowInsertY[j-1],
                                   'M', 'X', 'Y', &newRowMatch[j], tracebackPtr++);
                newRowMatch[j] += *(posteriorPtr++);

                // going to INSERT X state
                ChooseBestOfThree (oldRowMatch[j] + openingPenalty1[numActive1[i]][numGapOpens2[j]],
                                   oldRowInsertX[j] + continuingPenalty1[numActive1[i]],
                                   oldRowInsertY[j] + openingPenalty1[numActive1[i]][numGapOpens2[j]],
                                   'M', 'X', 'Y', &newRowInsertX[j], tracebackPtr++);

                // going to INSERT Y state
                ChooseBestOfThree (newRowMatch[j-1] + openingPenalty2[numActive2[j]][numGapOpens1[i]],
                                   newRowInsertX[j-1] + openingPenalty2[numActive2[j]][numGapOpens1[i]],
                                   newRowInsertY[j-1] + continuingPenalty2[numActive2[j]],
                                   'M', 'X', 'Y', &newRowInsertY[j], tracebackPtr++);
            }

            // swap rows
            float *temp;
            temp = oldRowMatch;
            oldRowMatch = newRowMatch;
            newRowMatch = temp;
            temp = oldRowInsertX;
            oldRowInsertX = newRowInsertX;
            newRowInsertX = temp;
            temp = oldRowInsertY;
            oldRowInsertY = newRowInsertY;
            newRowInsertY = temp;
        }

        // store best score
        float total;
        char matrix;
        ChooseBestOfThree (oldRowMatch[seq2Length], oldRowInsertX[seq2Length], oldRowInsertY[seq2Length],
                           'M', 'X', 'Y', &total, &matrix);

        delete [] twoRows;

        // compute traceback
        SafeVector<char> *alignment = new SafeVector<char>;
        assert (alignment);
        int r = seq1Length, c = seq2Length;
        while (r != 0 || c != 0)
        {

            int offset = (matrix == 'M') ? 0 : (matrix == 'X') ? 1 : 2;
            char ch = tracebackMatrix[(r*(seq2Length+1) + c) * 3 + offset];
            switch (matrix)
            {
            case 'Y':
                c--;
                alignment->push_back ('Y');
                break;
            case 'X':
                r--;
                alignment->push_back ('X');
                break;
            case 'M':
                c--;
                r--;
                alignment->push_back ('B');
                break;
            default:
                assert (false);
            }
            matrix = ch;
        }

        delete [] tracebackMatrix;

        reverse (alignment->begin(), alignment->end());

        return make_pair(alignment, 1.0f);
    }

    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::ComputeViterbiAlignment()
    //
    // Computes the highest probability pairwise alignment using the
    // probabilistic model.  The final alignment is returned as a
    //  pair consisting of:
    //    (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and
    //        denote insertions in one of the two sequences and
    //        B's denote that both sequences are present (i.e.
    //        matches).
    //    (2) a float containing the log probability of the best
    //        alignment (not used)
    /////////////////////////////////////////////////////////////////


    pair<SafeVector<char> *, float> ComputeViterbiAlignment (Sequence *seq1, Sequence *seq2) const
    {

        assert (seq1);
        assert (seq2);

        const int seq1Length = seq1->GetLength();
        const int seq2Length = seq2->GetLength();

        // retrieve the points to the beginning of each sequence
        SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
        SafeVector<char>::iterator iter2 = seq2->GetDataPtr();

        // create viterbi matrix
        VF *viterbiPtr = new VF (3 * (seq1Length+1) * (seq2Length+1), LOG_ZERO);
        assert (viterbiPtr);
        VF &viterbi = *viterbiPtr;

        // create traceback matrix
        VI *tracebackPtr = new VI (3 * (seq1Length+1) * (seq2Length+1), -1);
        assert (tracebackPtr);
        VI &traceback = *tracebackPtr;

        // initialization condition
        /*
            for (int k = 0; k < NumMatrixTypes; k++)
              viterbi[k] = initialDistribution[k];
        */
        viterbi[0] = LOG(0.6080327034);
        viterbi[1] = LOG(0.1959836632);
        viterbi[2] = LOG(0.1959836632);

        // remember offset for each index combination
        int ij = 0;
        int i1j = -seq2Length - 1;
        int ij1 = -1;
        int i1j1 = -seq2Length - 2;

        ij *= 3;
        i1j *= 3;
        ij1 *= 3;
        i1j1 *= 3;

        // compute viterbi scores
        for (int i = 0; i <= seq1Length; i++)
        {
            unsigned char c1 = (i == 0) ? '~' : (unsigned char) iter1[i];
            for (int j = 0; j <= seq2Length; j++)
            {
                unsigned char c2 = (j == 0) ? '~' : (unsigned char) iter2[j];

                if (i > 0 && j > 0)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        float newVal = viterbi[k + i1j1] + local_transProb[k][0] + matchProb[c1][c2];
                        if (viterbi[0 + ij] < newVal)
                        {
                            viterbi[0 + ij] = newVal;
                            traceback[0 + ij] = k;
                        }
                    }
                }
                if (i > 0)
                {
                    for (int k = 0; k < 1; k++)
                    {
                        float valFromMatch = insProb[c1][k] + viterbi[0 + i1j] + local_transProb[0][2*k+1];
                        float valFromIns = insProb[c1][k] + viterbi[2*k+1 + i1j] + local_transProb[2*k+1][2*k+1];
                        if (valFromMatch >= valFromIns)
                        {
                            viterbi[2*k+1 + ij] = valFromMatch;
                            traceback[2*k+1 + ij] = 0;
                        }
                        else
                        {
                            viterbi[2*k+1 + ij] = valFromIns;
                            traceback[2*k+1 + ij] = 2*k+1;
                        }
                    }
                }
                if (j > 0)
                {
                    for (int k = 0; k < 1; k++)
                    {
                        float valFromMatch = insProb[c2][k] + viterbi[0 + ij1] + local_transProb[0][2*k+2];
                        float valFromIns = insProb[c2][k] + viterbi[2*k+2 + ij1] + local_transProb[2*k+2][2*k+2];
                        if (valFromMatch >= valFromIns)
                        {
                            viterbi[2*k+2 + ij] = valFromMatch;
                            traceback[2*k+2 + ij] = 0;
                        }
                        else
                        {
                            viterbi[2*k+2 + ij] = valFromIns;
                            traceback[2*k+2 + ij] = 2*k+2;
                        }
                    }
                }

                ij += 3;
                i1j += 3;
                ij1 += 3;
                i1j1 += 3;
            }
        }

        // figure out best terminating cell
        float bestProb = LOG_ZERO;
        int state = -1;
        viterbi[0] = LOG(0.6080327034);
        viterbi[1] = LOG(0.1959836632);
        viterbi[2] = LOG(0.1959836632);

        for (int k = 0; k < 3; k++)
        {
            float thisProb = viterbi[k + 3 * ((seq1Length+1)*(seq2Length+1) - 1)] + viterbi[k];
            if (bestProb < thisProb)
            {
                bestProb = thisProb;
                state = k;
            }
        }
        assert (state != -1);

        delete viterbiPtr;

        // compute traceback
        SafeVector<char> *alignment = new SafeVector<char>;
        assert (alignment);
        int r = seq1Length, c = seq2Length;
        while (r != 0 || c != 0)
        {
            int newState = traceback[state + 3 * (r * (seq2Length+1) + c)];
            if (state == 0)
            {
                c--;
                r--;
                alignment->push_back ('B');
            }
            else if (state % 2 == 1)
            {
                r--;
                alignment->push_back ('X');
            }
            else
            {
                c--;
                alignment->push_back ('Y');
            }
            state = newState;
        }

        delete tracebackPtr;

        reverse (alignment->begin(), alignment->end());

        return make_pair(alignment, bestProb);
    }

    /////////////////////////////////////////////////////////////////
    // ProbabilisticModel::BuildPosterior()
    //
    // Builds a posterior probability matrix needed to align a pair
    // of alignments.  Mathematically, the returned matrix M is
    // defined as follows:
    //    M[i,j] =     sum          sum      f(s,t,i,j)
    //             s in align1  t in align2
    // where
    //                  [  P(s[i'] <--> t[j'])
    //                  [       if s[i'] is a letter in the ith column of align1 and
    //                  [          t[j'] it a letter in the jth column of align2
    //    f(s,t,i,j) =  [
    //                  [  0    otherwise
    //
    /////////////////////////////////////////////////////////////////


#ifdef _OPENMP
    struct SeqsPair
    {
        int seq1;
        int seq2;
    };
#endif

    VF *BuildPosterior (MultiSequence *align1, MultiSequence *align2,
                        const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                        float cutoff = 0.0f) const
    {
        const int seq1Length = align1->GetSequence(0)->GetLength();
        const int seq2Length = align2->GetSequence(0)->GetLength();

        VF *posteriorPtr = new VF((seq1Length+1) * (seq2Length+1), 0);
        assert (posteriorPtr);
        VF &posterior = *posteriorPtr;

#ifdef _OPENMP
        //calculate sequence pairs for openmp model
        int pairIdx = 0;
        int numPairs = align1->GetNumSequences() * align2->GetNumSequences();
        SeqsPair* seqsPairs = new SeqsPair[numPairs];
        for(int a = 0; a < align1->GetNumSequences(); a++)
        {
            for(int b = 0; b < align2->GetNumSequences(); b++)
            {
                seqsPairs[pairIdx].seq1 = a;
                seqsPairs[pairIdx].seq2 = b;
                pairIdx++;
            }
        }
#endif

        // do all pairwise
#ifdef _OPENMP
        #pragma omp parallel for private(pairIdx) default(shared) schedule(dynamic)
        for(pairIdx = 0; pairIdx < numPairs; pairIdx++)
        {
            int i = seqsPairs[pairIdx].seq1;
            int j = seqsPairs[pairIdx].seq2;

#else
        for ( int i = 0; i < align1->GetNumSequences(); i++ )
        {
            for ( int j = 0; j < align2->GetNumSequences(); j++ )
            {
#endif

            // for each s in align1
            int first = align1->GetSequence(i)->GetLabel();
            SafeVector<int> *mapping1 = align1->GetSequence(i)->GetMapping();
            // for each t in align2
            int second = align2->GetSequence(j)->GetLabel();
            SafeVector<int> *mapping2 = align2->GetSequence(j)->GetMapping();

            if (first < second)
            {

                // get the associated sparse matrix
                SparseMatrix *matrix = sparseMatrices[first][second];

                for (int ii = 1; ii <= matrix->GetSeq1Length(); ii++)
                {
                    SafeVector<PIF>::iterator row = matrix->GetRowPtr(ii);
                    int base = (*mapping1)[ii] * (seq2Length+1);
                    int rowSize = matrix->GetRowSize(ii);

                    // add in all relevant values
                    for (int jj = 0; jj < rowSize; jj++)
                        posterior[base + (*mapping2)[row[jj].first]] += row[jj].second;

                    // subtract cutoff
                    for (int jj = 0; jj < matrix->GetSeq2Length(); jj++)
                        posterior[base + (*mapping2)[jj]] -= cutoff;
                }

            }
            else
            {

                // get the associated sparse matrix
                SparseMatrix *matrix = sparseMatrices[second][first];

                for (int jj = 1; jj <= matrix->GetSeq1Length(); jj++)
                {
                    SafeVector<PIF>::iterator row = matrix->GetRowPtr(jj);
                    int base = (*mapping2)[jj];
                    int rowSize = matrix->GetRowSize(jj);

                    // add in all relevant values
                    for (int ii = 0; ii < rowSize; ii++)
                        posterior[base + (*mapping1)[row[ii].first] * (seq2Length + 1)] += row[ii].second;

                    // subtract cutoff
                    for (int ii = 0; ii < matrix->GetSeq2Length(); ii++)
                        posterior[base + (*mapping1)[ii] * (seq2Length + 1)] -= cutoff;
                }
            }

            delete mapping2;
            delete mapping1;
#ifndef _OPENMP
        }
#endif
    }

    return posteriorPtr;
}

//added by Liu Yongchao.Feb 23, 2010
VF *BuildPosterior(int* seqsWeights, MultiSequence *align1,
                   MultiSequence *align2,
                   const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                   float cutoff = 0.0f) const
{
    const int seq1Length = align1->GetSequence(0)->GetLength();
    const int seq2Length = align2->GetSequence(0)->GetLength();

    VF *posteriorPtr = new VF((seq1Length + 1) * (seq2Length + 1), 0);
    assert(posteriorPtr);
    VF &posterior = *posteriorPtr;

    //compute the total sum of all weights
    float totalWeights = 0;
    for (int i = 0; i < align1->GetNumSequences(); i++)
    {
        int first = align1->GetSequence(i)->GetLabel();
        int w1 = seqsWeights[first];
        for (int j = 0; j < align2->GetNumSequences(); j++)
        {
            int second = align2->GetSequence(j)->GetLabel();
            int w2 = seqsWeights[second];

            totalWeights += w1 * w2;
        }
    }
    // for each s in align1
    for (int i = 0; i < align1->GetNumSequences(); i++)
    {
        int first = align1->GetSequence(i)->GetLabel();
        int w1 = seqsWeights[first];
        SafeVector<int> *mapping1 = align1->GetSequence(i)->GetMapping();
        // for each t in align2
        for (int j = 0; j < align2->GetNumSequences(); j++)
        {
            int second = align2->GetSequence(j)->GetLabel();
            int w2 = seqsWeights[second];
            SafeVector<int> *mapping2 =
                align2->GetSequence(j)->GetMapping();

            float w = (float) (w1 * w2) / totalWeights;
            if (first < second)
            {

                // get the associated sparse matrix
                SparseMatrix *matrix = sparseMatrices[first][second];

                for (int ii = 1; ii <= matrix->GetSeq1Length(); ii++)
                {
                    SafeVector<PIF>::iterator row = matrix->GetRowPtr(ii);
                    int base = (*mapping1)[ii] * (seq2Length + 1);
                    int rowSize = matrix->GetRowSize(ii);

                    // add in all relevant values
                    for (int jj = 0; jj < rowSize; jj++)
                        posterior[base + (*mapping2)[row[jj].first]] += w
                                * row[jj].second;

                    // subtract cutoff
                    for (int jj = 0; jj < matrix->GetSeq2Length(); jj++)
                        posterior[base + (*mapping2)[jj]] -= w * cutoff;
                }

            }
            else
            {

                // get the associated sparse matrix
                SparseMatrix *matrix = sparseMatrices[second][first];

                for (int jj = 1; jj <= matrix->GetSeq1Length(); jj++)
                {
                    SafeVector<PIF>::iterator row = matrix->GetRowPtr(jj);
                    int base = (*mapping2)[jj];
                    int rowSize = matrix->GetRowSize(jj);

                    // add in all relevant values
                    for (int ii = 0; ii < rowSize; ii++)
                        posterior[base
                                  + (*mapping1)[row[ii].first]
                                  * (seq2Length + 1)] += w
                                                         * row[ii].second;

                    // subtract cutoff
                    for (int ii = 0; ii < matrix->GetSeq2Length(); ii++)
                        posterior[base + (*mapping1)[ii] * (seq2Length + 1)] -=
                            w * cutoff;
                }

            }

            delete mapping2;
        }

        delete mapping1;
    }

    return posteriorPtr;
}
};

#endif
