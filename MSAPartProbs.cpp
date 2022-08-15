////////////////////////////////////////////////////////////////
// MSAPartProbs.cpp
//
// Utilities for computation on global pair-HMM based on partition function
/////////////////////////////////////////////////////////////////

#include "SafeVector.h"
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include "MultiSequence.h"
#include "ScoreType.h"

#define  TRACE 0		// 0: NOTRACE 1: TRACE
//proba like settings
#define  endgaps 1		// 1: engap penaties enabled 0: disabled
#define  PART_FULL_MEMORY 0	//0: LOW MEM OPTION
#define  REVPART_FULL_MEMORY 0	//0: LOW MEM OPTION
using namespace std;

#ifdef _WIN32
#define OS_HUGE_VALL	HUGE_VAL
#else
#define OS_HUGE_VALL	HUGE_VALL
#endif

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

typedef struct sequence
{
    char *title;
    char *text;
    int length;
} fasta;

typedef struct alignment
{
    char *title;
    char *text;
    int length;
} align;

////////////////////////////////////////////////////////
//externs related to scoring matrix and input arguments
///////////////////////////////////////////////////////////
extern float g_gap_open1, g_gap_open2, g_gap_ext1, g_gap_ext2;
extern char aminos[26], matrixtype[20], bases[26];

extern double sub_matrix[26][26];
extern double normalized_matrix[26][26]; // be used to compute sequences' similarity
extern int subst_index[26];

extern float TEMPERATURE;
extern int MATRIXTYPE;

extern float GAPOPEN;
extern float GAPEXT;
extern argument_decl argument;

//////////////////////////////////////////////////////////////////////////////
//calculates reverse partition function values based on z matrices
//and also simulaneously calculates the propability of each basepair
//or aminoacid residue pair i,j
//////////////////////////////////////////////////////////////////////////////

VF *revers_partf(fasta sequences[2], const double termgapopen,
                 const double termgapextend, long double **Zfm, const double d,
                 const double e)
{
    // printf("revpart\n");
    //rest of the declarations
    int i, j;
    long double **Zm = NULL;
    long double **Ze = NULL;
    long double **Zf = NULL;
    int len0, len1;
    float probability;
    long double tempvar;
    int Si, Tj;
    double endgapopen, endgapextend;
    FILE *fo;

    //Init lengths of sequences
    len0 = strlen(sequences[0].text);
    len1 = strlen(sequences[1].text);

    //Safe vector declared
    VF *posteriorPtr = new VF((len0 + 1) * (len1 + 1));
    VF & posterior = *posteriorPtr;
    VF::iterator ptr = posterior.begin();

    if (TRACE)			//open the trace file
        fo = fopen("revpartdump", "a");

    //default:
    endgapopen = termgapopen;
    endgapextend = termgapextend;

    //instantiate the z matrix
    if (REVPART_FULL_MEMORY)
    {

        Ze = new long double *[sequences[1].length + 1];
        Zf = new long double *[sequences[1].length + 1];
        Zm = new long double *[sequences[1].length + 1];

        if (TRACE)
            printf("\n\n %e %e\n", d, e);

        //DYNAMICALLY GROW 2D Zm Zf Ze MARICES (long double)
        for (i = 0; i <= sequences[1].length; i++)
        {
            Ze[i] = new long double[sequences[0].length + 1];
            Zf[i] = new long double[sequences[0].length + 1];
            Zm[i] = new long double[sequences[0].length + 1];
        }
    }
    else
    {
        Zm = new long double *[2];
        Ze = new long double *[2];
        Zf = new long double *[2];
        for (i = 0; i <= 1; i++)
        {
            Zm[i] = new long double[sequences[0].length + 1];
            Ze[i] = new long double[sequences[0].length + 1];
            Zf[i] = new long double[sequences[0].length + 1];
        }

    }

    if (TRACE)
    {
        printf("in rev partf---");
        printf("\n\n");
    }

    if (REVPART_FULL_MEMORY)
    {
        for (i = 0; i <= len1; i++)
            for (j = 0; j <= len0; j++)
            {
                Zm[i][j] = 0.0;
                Zf[i][j] = 0.0;
                Ze[i][j] = 0.0;
            }
    }
    else
    {

        for (j = 0; j <= len0; j++)
        {
            Zm[0][j] = 0;
            Zf[0][j] = 0;
            Ze[0][j] = 0;
            Zf[1][j] = 0;
            Ze[1][j] = 0;
            Zm[1][j] = 0;
        }
    }

    //fill the probability matrix with 0s
    for (i = 0; i <= len1; i++)
        for (j = 0; j <= len0; j++)
            ptr[j * (len1 + 1) + i] = 0;

    if (endgaps == 0)
    {
        Zm[len1][len0] = 1;
        Ze[len1][len0] = Zf[len1][len0] = 0;
        Zf[len1 - 1][len0] = Zm[len1][len0] * d;
        Ze[len1][len0 - 1] = Zm[len1][len0] * d;

        //>=2ND ROW INIT
        if (REVPART_FULL_MEMORY)
        {
            for (i = len1 - 2; i >= 0; i--)
            {
                Zf[i][len0] = Zf[i + 1][len0] * e;
            }
        }

        //>=2ND COL INIT
        if (REVPART_FULL_MEMORY)
        {
            for (j = len0 - 2; j >= 0; j--)
            {
                Ze[len1][j] = Ze[len1][j + 1] * e;
            }
        }
        else
        {
            for (j = len0 - 2; j >= 0; j--)
            {
                Ze[0][j] = Ze[0][j + 1] * e;
            }
        }
    }
    else
    {

        if (REVPART_FULL_MEMORY)
        {

            Zm[len1][len0] = 1;
            Ze[len1][len0] = Zf[len1][len0] = 0;
            Zf[len1 - 1][len0] = Zm[len1][len0] * endgapopen;
            Ze[len1][len0 - 1] = Zm[len1][len0] * endgapopen;

            //>=2ND ROW INIT
            for (i = len1 - 2; i >= 0; i--)
            {
                Zf[i][len0] = Zf[i + 1][len0] * endgapextend;
            }

            //M Iy= d+j*e

            //>=2ND COL INIT
            for (j = len0 - 2; j >= 0; j--)
            {
                Ze[len1][j] = Ze[len1][j + 1] * endgapextend;
            }

        }
        else
        {
            //in Zm
            //let:
            //  Zm(0) be the current row being filled/computed
            //  Zm(1) be the previous row

            Zm[1][len0] = 1;
            Ze[0][len0] = Zf[0][len0] = 0;
            Zf[1][len0] = Zm[1][len0] * endgapopen;
            Ze[0][len0 - 1] = Zm[1][len0] * endgapopen;

            //>=2ND COL INIT
            for (j = len0 - 2; j >= 0; j--)
            {
                Ze[0][j] = Ze[0][j + 1] * endgapextend;
            }

        }			//END ELSE

    }				//END FULL MEMORY and GAP enablement IF STATEMENT

    double scorez, zz = 0;

    for (i = len1 - 1; i >= 0; i--)
    {

        for (j = len0 - 1; j >= 0; j--)
        {
            Si = subst_index[sequences[1].text[i] - 'A'];
            Tj = subst_index[sequences[0].text[j] - 'A'];
            scorez = sub_matrix[Si][Tj];

            //endgaps modification aug 10
            double open0, extend0, open1, extend1;

            open0 = open1 = d;
            extend0 = extend1 = e;

            if (endgaps == 1)
            {

                //check to see if one of the 2 sequences or both reach the end

                if (i == 0)
                {
                    open0 = endgapopen;
                    extend0 = endgapextend;

                }

                if (j == 0)
                {
                    open1 = endgapopen;
                    extend1 = endgapextend;
                }

            }

            if (REVPART_FULL_MEMORY)
            {
                //z computation

                Ze[i][j] = Zm[i][j + 1] * open0 + Ze[i][j + 1] * extend0;
                Zf[i][j] = Zm[i + 1][j] * open1 + Zf[i + 1][j] * extend1;
                Zm[i][j] = (Zm[i + 1][j + 1] + Zf[i + 1][j + 1]
                            + Ze[i + 1][j + 1]) * scorez;
                zz = Zm[i][j] + Zf[i][j] + Ze[i][j];

            }
            else
            {

                //2 ROW zE zF ALGORITHM GOES...:
                //Ze[1][j] =Zm[i][j + 1] * exp(beta * open0) + Ze[1][j + 1] *exp(beta * extend0);
                //Zf[1][j] = Zm[i + 1][j] * exp(beta * open1) + Zf[0][j] * exp(beta * extend1);
                //Zm[i][j] = (Zm[i + 1][j + 1] + Zf[0][j + 1] + Ze[0][j + 1]) * exp(beta * scorez);
                //zz = Zm[0][j] + Zf[1][j] + Ze[1][j];

                //lowmem code for merging probability calculating module
                //Here we make use of Zm as a 2 row matrix

                Zf[1][j] = Zm[1][j] * open1 + Zf[0][j] * extend1;
                Ze[1][j] = Zm[0][j + 1] * open0 + Ze[1][j + 1] * extend0;
                Zm[0][j] = (Zm[1][j + 1] + Zf[0][j + 1] + Ze[0][j + 1])
                           * scorez;

                tempvar = Zfm[i + 1][j + 1] * Zm[0][j];
                //divide P(i,j) i.e. pairwise probability by denominator
                tempvar /= (scorez * Zfm[0][0]);
                probability = (float) tempvar;

                //store only noticable probabilities
                //if (probability <= 1 && probability >= 0.001) {
                //algorithm goes...
                //validprob[i + 1][j + 1] = probability;
                ptr[(j + 1) * (len1 + 1) + (i + 1)] = probability;
                //}
                //lowmem code ends here

            }

        }			//end of for

        if (REVPART_FULL_MEMORY == 0)
        {
            for (int t = 0; t <= sequences[0].length; t++)
            {
                Ze[0][t] = Ze[1][t];
                Ze[1][t] = 0;

                Zf[0][t] = Zf[1][t];
                Zf[1][t] = 0;

                Zm[1][t] = Zm[0][t];
                Zm[0][t] = 0;

            }
            Zf[0][len0] = 1;

        }

    }				//end of for

    if (TRACE)
    {
        printf("\n\nrM:....\n\n");
        if (REVPART_FULL_MEMORY)
        {
            for (i = 0; i <= len1; i++)
            {
                for (j = 0; j <= len0; j++)
                    printf("%.2Le ", Zm[i][j]);
                printf("\n");
            }

            printf("\n\nrE:....\n\n");
            for (i = 0; i <= len1; i++)
            {
                for (j = 0; j <= len0; j++)
                    printf("%.2Le ", Ze[i][j]);
                printf("\n");

            }

            printf("\n\nrF:....\n\n");
            for (i = 0; i <= len1; i++)
            {
                for (j = 0; j <= len0; j++)
                    printf("%.2Le ", Zf[i][j]);
                printf("\n");

            }

        }

    }

    if (TRACE)
    {
        fprintf(fo, "\n");
        fclose(fo);
    }

    //delete unused memory

    if (REVPART_FULL_MEMORY)
    {
        for (i = 0; i <= len1; i++)
        {
            delete (Zm[i]);
            delete (Zf[i]);
            delete (Ze[i]);
        }
    }
    else
    {
        delete (Zf[0]);
        delete (Ze[0]);
        delete (Zm[0]);

        delete (Zm[1]);
        delete (Zf[1]);
        delete (Ze[1]);
    }

    for (i = 0; i <= len1; i++)
    {
        delete (Zfm[i]);
    }

    if (Zf != NULL)
        delete (Zf);

    if (Ze != NULL)
        delete (Ze);

    if (Zm != NULL)
        delete (Zm);

    if (Zfm != NULL)
        delete (Zfm);

    posterior[0] = 0;
    return (posteriorPtr);

}

//////////////////////////////////////////////////////////////
//forward partition function
/////////////////////////////////////////////////////////////

long double **partf(fasta sequences[2], const double termgapopen,
                    const double termgapextend, const double d, const double e)
{
    //printf("partf\n");
    int i, j, len1, len0;
    long double **Zm = NULL, **Zf = NULL, **Ze = NULL, zz = 0;
    double endgapopen, endgapextend;

    //default:
    endgapopen = termgapopen;
    endgapextend = termgapextend;

    //the flag endgaps is set at the #define section
    if (PART_FULL_MEMORY)
    {

        Zf = new long double *[sequences[1].length + 1];
        Ze = new long double *[sequences[1].length + 1];
        Zm = new long double *[sequences[1].length + 1];

        //comment
        if (TRACE)
            printf("\nPARTF:====\n");

        //DYNAMICALLY GROW 2D M,IX,IY,PIX,PIY MARICES
        for (i = 0; i <= sequences[1].length; i++)
        {
            Zf[i] = new long double[sequences[0].length + 1];
            Ze[i] = new long double[sequences[0].length + 1];
            Zm[i] = new long double[sequences[0].length + 1];
        }
    }
    else
    {
        Zm = new long double *[sequences[1].length + 1];
        Ze = new long double *[2];
        Zf = new long double *[2];
        for (i = 0; i <= sequences[1].length; i++)
        {
            Zm[i] = new long double[sequences[0].length + 1];
        }
        Ze[0] = new long double[sequences[0].length + 1];
        Zf[0] = new long double[sequences[0].length + 1];
        Ze[1] = new long double[sequences[0].length + 1];
        Zf[1] = new long double[sequences[0].length + 1];
    }

    len0 = strlen(sequences[0].text);
    len1 = strlen(sequences[1].text);

    if (PART_FULL_MEMORY)
    {
        for (i = 0; i <= sequences[1].length; i++)
            for (j = 0; j <= sequences[0].length; j++)
            {
                Zm[i][j] = 0.00;
                Zf[i][j] = 0.00;
                Ze[i][j] = 0.00;
            }
    }
    else
    {
        for (i = 0; i <= len1; i++)
        {
            for (j = 0; j <= len0; j++)
            {
                Zm[i][j] = 0;
            }
        }
        for (j = 0; j <= len0; j++)
        {
            Zf[0][j] = 0;
            Ze[0][j] = 0;
            Zf[1][j] = 0;
            Ze[1][j] = 0;
        }
    }

    //INTITIALIZE THE DP

    if (endgaps == 0)
    {
        Zm[0][0] = 1.00;

        Zf[0][0] = Ze[0][0] = 0;
        Zf[1][0] = Zm[0][0] * d;
        Ze[0][1] = Zm[0][0] * d;

        //>=2ND ROW INIT
        if (PART_FULL_MEMORY)
        {
            for (i = 2; i <= sequences[1].length; i++)
            {
                Zf[i][0] = Zf[i - 1][0] * e;
            }
        }

        //>=2ND COL INIT
        for (j = 2; j <= sequences[0].length; j++)
        {
            Ze[0][j] = Ze[0][j - 1] * e;
        }
    }
    else
    {
        //init z
        Zm[0][0] = 1.00;
        Zf[0][0] = Ze[0][0] = 0;
        Zf[1][0] = Zm[0][0] * endgapopen;
        Ze[0][1] = Zm[0][0] * endgapopen;

        //>=2ND ROW INIT
        if (PART_FULL_MEMORY)
        {
            for (i = 2; i <= sequences[1].length; i++)
            {
                Zf[i][0] = Zf[i - 1][0] * endgapextend;
            }
        }

        //>=2ND COL INIT
        for (j = 2; j <= sequences[0].length; j++)
        {
            Ze[0][j] = Ze[0][j - 1] * endgapextend;
        }
    }

    //1ST ROW/COL INIT

    int Si, Tj;
    double score;

    for (i = 1; i <= sequences[1].length; i++)
    {

        for (j = 1; j <= sequences[0].length; j++)
        {

            Si = subst_index[sequences[1].text[i - 1] - 'A'];
            Tj = subst_index[sequences[0].text[j - 1] - 'A'];

            score = sub_matrix[Si][Tj];

            double open0, extend0, open1, extend1;

            open0 = open1 = d;
            extend0 = extend1 = e;

            if (endgaps == 1)
            {
                //check to see if one of the 2 sequences or both reach the end

                if (i == sequences[1].length)
                {
                    open0 = endgapopen;
                    extend0 = endgapextend;

                }

                if (j == sequences[0].length)
                {
                    open1 = endgapopen;
                    extend1 = endgapextend;
                }
            }

            //
            //z computation using open and extend temp vars
            //open0 is gap open in seq0 and open1 is gap open in seq1
            //entend0 is gap extend in seq0 and extend1 is gap extend in seq1

            if (PART_FULL_MEMORY)
            {
                Ze[i][j] = Zm[i][j - 1] * open0 + Ze[i][j - 1] * extend0;

                if (Ze[i][j] >= OS_HUGE_VALL)
                {
                    printf("ERROR: huge val error for Ze\n");
                    exit(1);
                }

                Zf[i][j] = Zm[i - 1][j] * open1 + Zf[i - 1][j] * extend1;

                if (Zf[i][j] >= OS_HUGE_VALL)
                {
                    printf("ERROR: huge val error for Zf\n");
                    exit(1);
                }

                Zm[i][j] = (Zm[i - 1][j - 1] + Ze[i - 1][j - 1]
                            + Zf[i - 1][j - 1]) * score;

                if (Zm[i][j] >= OS_HUGE_VALL)
                {
                    printf("ERROR: huge val error for Zm\n");
                    exit(1);
                }

                zz = Zm[i][j] + Ze[i][j] + Zf[i][j];
            }
            else
            {
                Ze[1][j] = Zm[i][j - 1] * open0 + Ze[1][j - 1] * extend0;

                if (Ze[1][j] >= OS_HUGE_VALL)
                {
                    printf("ERROR: huge val error for zE\n");
                    exit(1);
                }

                Zf[1][j] = Zm[i - 1][j] * open1 + Zf[0][j] * extend1;

                if (Zf[1][j] >= OS_HUGE_VALL)
                {
                    printf("ERROR: huge val error for zF\n");
                    exit(1);
                }

                Zm[i][j] = (Zm[i - 1][j - 1] + Ze[0][j - 1] + Zf[0][j - 1])
                           * score;

                if (Zm[i][j] >= OS_HUGE_VALL)
                {
                    printf("ERROR: huge val error for zM\n");
                    exit(1);
                }

                zz = Zm[i][j] + Ze[1][j] + Zf[1][j];
            }

        }			//end for

        if (!PART_FULL_MEMORY)
        {
            for (int t = 0; t <= sequences[0].length; t++)
            {
                Ze[0][t] = Ze[1][t];
                Ze[1][t] = 0;

                Zf[0][t] = Zf[1][t];
                Zf[1][t] = 0;
            }

            Zf[1][0] = 1;

        }

    }				//end for

    //store the sum of zm zf ze (m,n)s in zm's 0,0 th position
    Zm[0][0] = zz;

    if (TRACE)
    {
        //debug code aug 3
        //print the 3 Z matrices namely Zm Zf and Ze

        printf("\n\nFINAL Zm:\n");
        for (i = 0; i <= sequences[1].length; i++)
        {
            for (j = 0; j <= sequences[0].length; j++)
                printf("%.2Le ", Zm[i][j]);
            printf("\n");
        }

        printf("FINAL Zf \n");
        for (i = 0; i <= sequences[1].length; i++)
        {
            for (j = 0; j <= sequences[0].length; j++)
                printf("%.2Le ", Zf[i][j]);
            printf("\n");
        }

        printf("FINAL Ze \n");
        for (i = 0; i <= sequences[1].length; i++)
        {
            for (j = 0; j <= sequences[0].length; j++)
                printf("%.2Le ", Ze[i][j]);
            printf("\n");
        }

        //end debug dump code

    }

    if (PART_FULL_MEMORY)
    {
        for (i = 0; i <= sequences[1].length; i++)
        {
            delete (Zf[i]);
            delete (Ze[i]);
        }
    }
    else
    {
        delete (Zf[0]);
        delete (Ze[0]);
        delete (Zf[1]);
        delete (Ze[1]);
    }

    delete (Zf);
    delete (Ze);

    return Zm;

}				//end of forward partition function

/////////////////////////////////////////////////////////////////////////////////////////
//entry point (was the main function) , returns the posterior probability safe vector
////////////////////////////////////////////////////////////////////////////////////////
VF *ComputePostProbs(int a, int b, string seq1, string seq2)
{
    //printf("probamod\n");
    double gap_open = -22, gap_ext = -1, beta = 0.2;//T = 5, beta = 1/T = 0.2, by default
    int stock_loop = 1;
    int le = 160;
    double termgapopen = 1.0f;	//exp(0)
    double termgapextend = 1.0f;	//exp(0)

    //initialize the sequence structure
    fasta sequences[2];

    sequences[0].length = strlen((char *) seq1.c_str());
    sequences[0].text = (char *) seq1.c_str();
    sequences[0].title = new char[10];
    strcpy(sequences[0].title, "seq0");
    sequences[1].length = strlen((char *) seq2.c_str());
    sequences[1].text = (char *) seq2.c_str();
    sequences[1].title = new char[10];
    strcpy(sequences[1].title, "seq1");

    if (TRACE)

    {
        printf("%d %d %s\n%d %d %s\n--\n", a, sequences[0].length,
               sequences[0].text, b, sequences[1].length, sequences[1].text);
        printf("after init\n");

        FILE *dump1 = fopen("dump1", "a");
        fprintf(dump1, "%d %d %s\n%d %d %s\n--\n", a, sequences[0].length,
                sequences[0].text, b, sequences[1].length, sequences[1].text);
        fclose(dump1);
    }

    gap_open = argument.gapopen;
    gap_ext = argument.gapext;
    beta = argument.beta;

    stock_loop = argument.N;
    le = argument.matrix;

    //compute the values of exp(beta * ?)
    termgapopen = exp(beta * 0.0);
    termgapextend = exp(beta * 0.0);
    gap_open = exp(beta * gap_open);
    gap_ext = exp(beta * gap_ext);

    if (TRACE)
        printf("%f %f %f %d\n", gap_open, gap_ext, beta, le);

    //call for calculating the posterior probabilities
    // 1. call partition function partf
    // 2. calculate revpartition using revers_parf
    // 3. calculate probabilities
    /// MODIFICATION... POPULATE SAFE VECTOR

    long double **MAT1;

    MAT1 = partf(sequences, termgapopen, termgapextend, gap_open, gap_ext);

    return revers_partf(sequences, termgapopen, termgapextend, MAT1, gap_open,
                        gap_ext);

}

//////////////////////////////////////////////////////////////
//Compute Viterbi Alignment
/////////////////////////////////////////////////////////////

pair<SafeVector<char> *, float> partViterbi(string seq1, string seq2)
{


    double gap_open = -12, gap_ext = -1, beta = 0.2;//T = 5, beta = 1/T = 0.2, by default
    int stock_loop = 1;
    int le = 160;
    //double termgapopen = 1.0f;	//exp(0)
    //double termgapextend = 1.0f;	//exp(0)

    //initialize the sequence structure
    fasta sequences[2];
    sequences[0].length = strlen((char *) seq1.c_str());
    sequences[0].text = (char *) seq1.c_str();
    sequences[0].title = new char[10];
    strcpy(sequences[0].title, "seq0");
    sequences[1].length = strlen((char *) seq2.c_str());
    sequences[1].text = (char *) seq2.c_str();
    sequences[1].title = new char[10];
    strcpy(sequences[1].title, "seq1");

    gap_open = argument.gapopen;
    gap_ext = argument.gapext;
    beta = argument.beta;

    stock_loop = argument.N;
    le = argument.matrix;

    //compute the values of exp(beta * ?)
    double endgapopen = exp(beta * 0.0);
    double endgapextend = exp(beta * 0.0);
    double d = exp(beta * gap_open);
    double e = exp(beta * gap_ext);

    int i, j, len1, len0;
    long double **Zm = NULL, **Zf = NULL, **Ze = NULL;
    int **traceZm = NULL, **traceZf = NULL, **traceZe = NULL;

    //the flag endgaps is set at the #define section
    Zf = new long double *[sequences[1].length + 1];
    Ze = new long double *[sequences[1].length + 1];
    Zm = new long double *[sequences[1].length + 1];

    traceZf = new int *[sequences[1].length + 1];
    traceZe = new int *[sequences[1].length + 1];
    traceZm = new int *[sequences[1].length + 1];

    //DYNAMICALLY GROW 2D M,IX,IY,PIX,PIY MARICES
    for (i = 0; i <= sequences[1].length; i++)
    {
        Zf[i] = new long double[sequences[0].length + 1];
        Ze[i] = new long double[sequences[0].length + 1];
        Zm[i] = new long double[sequences[0].length + 1];

        traceZf[i] = new int[sequences[0].length + 1];
        traceZe[i] = new int[sequences[0].length + 1];
        traceZm[i] = new int[sequences[0].length + 1];
    }

    len0 = strlen(sequences[0].text);
    len1 = strlen(sequences[1].text);


    for (i = 0; i <= sequences[1].length; i++)
        for (j = 0; j <= sequences[0].length; j++)
        {
            Zm[i][j] = 0.00;
            Zf[i][j] = 0.00;
            Ze[i][j] = 0.00;

            traceZm[i][j] = -1;
            traceZf[i][j] = -1;
            traceZe[i][j] = -1;
        }


    //INTITIALIZE THE DP
    if (endgaps == 0)
    {
        Zm[0][0] = 1.00;

        Zf[0][0] = Ze[0][0] = 0;
        Zf[1][0] = Zm[0][0] * d;
        Ze[0][1] = Zm[0][0] * d;

        //>=2ND ROW INIT

        for (i = 2; i <= sequences[1].length; i++)
        {
            Zf[i][0] = Zf[i - 1][0] * e;
            traceZf[i][0] = 2;
        }


        //>=2ND COL INIT
        for (j = 2; j <= sequences[0].length; j++)
        {
            Ze[0][j] = Ze[0][j - 1] * e;
            traceZe[0][j] = 1;
        }
    }
    else
    {
        //init z
        Zm[0][0] = 1.00;
        Zf[0][0] = Ze[0][0] = 0;
        Zf[1][0] = Zm[0][0] * endgapopen;
        Ze[0][1] = Zm[0][0] * endgapopen;

        //>=2ND ROW INIT

        for (i = 2; i <= sequences[1].length; i++)
        {
            Zf[i][0] = Zf[i - 1][0] * endgapextend;
            traceZf[i][0] = 2;
        }
        //>=2ND COL INIT
        for (j = 2; j <= sequences[0].length; j++)
        {
            Ze[0][j] = Ze[0][j - 1] * endgapextend;
            traceZe[0][j] = 1;
        }
    }

    //1ST ROW/COL INIT

    int Si, Tj;
    double score;

    for (i = 1; i <= sequences[1].length; i++)
    {

        for (j = 1; j <= sequences[0].length; j++)
        {

            Si = subst_index[sequences[1].text[i - 1] - 'A'];
            Tj = subst_index[sequences[0].text[j - 1] - 'A'];

            score = sub_matrix[Si][Tj];

            double open0, extend0, open1, extend1;

            open0 = open1 = d;
            extend0 = extend1 = e;

            if (endgaps == 1)
            {
                //check to see if one of the 2 sequences or both reach the end

                if (i == sequences[1].length)
                {
                    open0 = endgapopen;
                    extend0 = endgapextend;

                }

                if (j == sequences[0].length)
                {
                    open1 = endgapopen;
                    extend1 = endgapextend;
                }
            }

            //
            //z computation using open and extend temp vars
            //open0 is gap open in seq0 and open1 is gap open in seq1
            //entend0 is gap extend in seq0 and extend1 is gap extend in seq1
            Zf[i][j] = Zf[i - 1][j] * extend1;
            traceZf[i][j] = 2;

            if(Zm[i - 1][j] * open1 > Zf[i][j])
            {
                Zf[i][j] = Zm[i - 1][j] * open1;
                traceZf[i][j] = 0;
            }
            if (Zf[i][j] >= OS_HUGE_VALL)
            {
                printf("ERROR: huge val error for Zf\n");
                exit(1);
            }
            Ze[i][j] = Ze[i][j - 1] * extend0;
            traceZe[i][j] = 1;
            if(Zm[i][j - 1] * open0 > Ze[i][j])
            {
                Ze[i][j] = Zm[i][j - 1] * open0;
                traceZe[i][j] = 0;
            }

            if (Ze[i][j] >= OS_HUGE_VALL)
            {
                printf("ERROR: huge val error for Ze\n");
                exit(1);
            }

            Zm[i][j] = Zm[i - 1][j - 1] * score;
            traceZm[i][j] = 0;
            if(Zf[i - 1][j - 1] * score > Zm[i][j])
            {
                Zm[i][j] = Zf[i - 1][j - 1] * score;
                traceZm[i][j] = 2;
            }
            if(Ze[i - 1][j - 1] * score > Zm[i][j])
            {
                Zm[i][j] = Ze[i - 1][j - 1] * score;
                traceZm[i][j] = 1;
            }
            if (Zm[i][j] >= OS_HUGE_VALL)
            {
                printf("ERROR: huge val error for Zm\n");
                exit(1);
            }

        }//end for
    }//end for
    // figure out best terminating cell

    float bestProb = Zm[sequences[1].length][sequences[0].length];
    int state = 0;
    if( bestProb < Zf[sequences[1].length][sequences[0].length])
    {
        bestProb = Zf[sequences[1].length][sequences[0].length];
        state = 2;
    }
    if( bestProb < Ze[sequences[1].length][sequences[0].length])
    {
        bestProb = Ze[sequences[1].length][sequences[0].length];
        state = 1;
    }
    assert (state != -1);

    // compute traceback
    SafeVector<char> *alignment = new SafeVector<char>;
    assert (alignment);
    int c = sequences[1].length, r = sequences[0].length;
    while (r != 0 || c != 0)
    {
        int newState;
        if(state == 0)
        {
            newState = traceZm[c][r];
            c--;
            r--;
            alignment->push_back ('B');
        }
        else if(state == 1)
        {
            newState = traceZe[c][r];
            r--;
            alignment->push_back ('X');
        }
        else
        {
            newState = traceZf[c][r];
            c--;
            alignment->push_back ('Y');
        }
        state = newState;
    }

    reverse (alignment->begin(), alignment->end());

    for (i = 0; i <= sequences[1].length; i++)
    {
        delete (Zf[i]);
        delete (Ze[i]);
        delete (Zm[i]);
        delete (traceZf[i]);
        delete (traceZe[i]);
        delete (traceZm[i]);
    }

    delete (Zf);
    delete (Ze);
    delete (Zm);
    delete (traceZf);
    delete (traceZe);
    delete (traceZm);

    return make_pair(alignment, bestProb);
}

//////////////////////////////////////////////////////////////
// Compute two sequences' similarity defined as the normalized alignment score without gap penalties
/////////////////////////////////////////////////////////////

float computeSimilarity(string seq1, string seq2, SafeVector<char> * alignment)
{

    //initialize the sequence structure
    fasta sequences[2];
    sequences[0].length = strlen((char *) seq1.c_str());
    sequences[0].text = (char *) seq1.c_str();
    sequences[0].title = new char[10];
    strcpy(sequences[0].title, "seq0");
    sequences[1].length = strlen((char *) seq2.c_str());
    sequences[1].text = (char *) seq2.c_str();
    sequences[1].title = new char[10];
    strcpy(sequences[1].title, "seq1");

    float bestProb = 0;
    int Si, Tj;
    double score;
    int i = 1;
    int j = 1;
    for (SafeVector<char>::iterator iter = alignment->begin();
            iter != alignment->end(); ++iter)
    {
        if (*iter == 'B')
        {
            Si = subst_index[sequences[1].text[j - 1] - 'A'];
            Tj = subst_index[sequences[0].text[i - 1] - 'A'];
            score = normalized_matrix[Si][Tj];
            bestProb += score;
            i++;
            j++;
        }
        else if(*iter == 'X') i++;
        else if(*iter == 'Y') j++;
    }
    if(i!= sequences[0].length + 1 || j!= sequences[1].length + 1 ) cerr << "similarity error"<< endl;
    bestProb /= alignment->size();
    //bestProb /= min(sequences[0].length, sequences[1].length);
    return bestProb;
}
//end of posterior probability  module
