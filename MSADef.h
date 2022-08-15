#ifndef _MSA_DEF_H
#define _MSA_DEF_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

//maximum number
#define MAX_INT_NUM				0x7FFFFFFF
#define MAX_FLOAT_NUM			FLT_MAX
#define INT_MULTIPLY			1000

#define SUBMATRIX_INT_SCALE		100

//a tree node is a leaf or a node
enum
{
    NONE, NODE, LEAF
};

#endif

