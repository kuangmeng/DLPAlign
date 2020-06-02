/////////////////////////////////////////////////////////////////
// MSAClusterTree.cpp
//
// Routines for UPGMA guide tree
//
/////////////////////////////////////////////////////////////////

#include "MSAClusterTree.h"

MSAClusterTree::MSAClusterTree(MSA* msa, VVF& distMatrix, int numSeqs) :
    MSAGuideTree(msa, distMatrix, numSeqs)
{
}

MSAClusterTree::~MSAClusterTree()
{
}

void MSAClusterTree::create()
{
    //generate the neighbor-joining tree
    this->generateClusterTree();

    //calculate sequence weights
    this->getSeqsWeights();

    //construct the alignment orders
    this->createAlignmentOrders();
}

void MSAClusterTree::create(int second_step)
{
    //generate the neighbor-joining tree
    this->generateClusterTree(second_step);

    //calculate sequence weights
    this->getSeqsWeights();

    //construct the alignment orders
    this->createAlignmentOrders();
}

void MSAClusterTree::generateClusterTree()
{
    int i;
    ValidNode* validNodes, *headValidNodes;
    ValidNode* miniPtr, *minjPtr, *ivalid, *jvalid;
    int mini, minj;
    float* joins;
    unsigned int* clusterLeafs;

    //initialize the valid nodes link list
    validNodes = new ValidNode[leafsNum + 1];
    joins = new float[leafsNum + 1];
    clusterLeafs = new unsigned int[nodesNum + 1];
    if (!validNodes || !joins || !clusterLeafs)
    {
        cerr << "Out of memory of the reconstruction of cluster tree" << endl;
    }
    //initialize cluster size
    for (i = 0; i < this->leafsNum; i++)
    {
        clusterLeafs[i] = 1;
    }

    headValidNodes = &validNodes[0];
    headValidNodes->next = &validNodes[1];
    headValidNodes->n = -1;
    headValidNodes->node = -1;
    headValidNodes->prev = NULL;

    //build an initial link list
    ValidNode* curr = &validNodes[1];
    ValidNode* prev = headValidNodes;
    ValidNode* next = &validNodes[2];
    for (i = 0; i < leafsNum; i++)
    {
        curr->n = i;
        curr->node = i;
        curr->prev = prev;
        curr->next = next;
        prev = curr;
        curr = next;
        next++;
    }
    prev->next = NULL;

    //to generate the cluster tree
    int nodeIdx;	//the index of an internal node
    int firstNode = leafsNum;	//the index of the first internal node
    int lastNode = firstNode + leafsNum - 1;//the index of the last internal node

    for (nodeIdx = firstNode; nodeIdx < lastNode; nodeIdx++)
    {
        //find closest pair of clusters
        float minDist = 1.1f;
        miniPtr = headValidNodes;
        minjPtr = headValidNodes;

        for (ivalid = headValidNodes->next; ivalid != NULL;
                ivalid = ivalid->next)
        {
            mini = ivalid->n;

            for (jvalid = headValidNodes->next;
                    jvalid != NULL && jvalid->n < mini; jvalid = jvalid->next)
            {
                minj = jvalid->n;
                float dist = (*distMatrix)[mini][minj];
                if (dist < 0)
                {
                    cerr
                            << "ERROR: It is impossible to have distance value less than zero"
                            << endl;
                    dist = 0;
                }
                if (dist < minDist)
                {
                    minDist = dist;
                    miniPtr = ivalid;
                    minjPtr = jvalid;
                }
                //printf("dist %g mini %d minj %d\n", dist, ivalid->node, jvalid->node);
            }
        }
        //printf("**** mini %d minj %d minDist %g *****\n", miniPtr->node, minjPtr->node, minDist);
        //check the validity of miniPtr and minjPtr;
        if (miniPtr == headValidNodes || minjPtr == headValidNodes)
        {
            cerr << "OOPS: Error occurred while constructing the cluster tree\n"
                 << endl;
            exit(-1);
        }
        //computing branch length and join the two nodes
        float branchLength = minDist * 0.5f;
        this->connectNodes(&nodes[nodeIdx], nodeIdx, &nodes[miniPtr->node],
                           branchLength, &nodes[minjPtr->node], branchLength);
        clusterLeafs[nodeIdx] = clusterLeafs[miniPtr->node]
                                + clusterLeafs[minjPtr->node];

        //remove the valid node minjPtr from the list
        minjPtr->prev->next = minjPtr->next;
        if (minjPtr->next != NULL)
        {
            minjPtr->next->prev = minjPtr->prev;
        }
        minjPtr->prev = minjPtr->next = NULL;

        //compute the distance of each remaining valid node to the new node
        for (ivalid = headValidNodes->next; ivalid != NULL;
                ivalid = ivalid->next)
        {
            int idx = ivalid->n;

            float idist = (*distMatrix)[miniPtr->n][idx];
            float jdist = (*distMatrix)[minjPtr->n][idx];

            unsigned int isize = clusterLeafs[miniPtr->node];
            unsigned int jsize = clusterLeafs[minjPtr->node];
            joins[idx] = (idist * isize + jdist * jsize) / (isize + jsize);
            //joins[idx] = (idist + jdist )/ 2;
        }
        //update the distance to the new node
        miniPtr->node = nodeIdx;
        mini = miniPtr->n;
        for (jvalid = headValidNodes->next; jvalid != NULL;
                jvalid = jvalid->next)
        {
            minj = jvalid->n;

            float dist = joins[minj];
            (*distMatrix)[mini][minj] = dist;
            (*distMatrix)[minj][mini] = dist;
        }
    }
    //add a pseudo root to this unrooted NJ tree
    this->root = &nodes[lastNode - 1];

    delete[] validNodes;
    delete[] joins;
    delete[] clusterLeafs;
}

//new add
void MSAClusterTree::generateClusterTree(int second_step)
{
    int i;
    ValidNode* validNodes, *headValidNodes;
    ValidNode* miniPtr, *minjPtr, *ivalid, *jvalid;
    int mini, minj;
    float* joins;
    unsigned int* clusterLeafs;

    //initialize the valid nodes link list
    validNodes = new ValidNode[leafsNum + 1];
    joins = new float[leafsNum + 1];
    clusterLeafs = new unsigned int[nodesNum + 1];
    if (!validNodes || !joins || !clusterLeafs)
    {
        cerr << "Out of memory of the reconstruction of cluster tree" << endl;
    }
    //initialize cluster size
    for (i = 0; i < this->leafsNum; i++)
    {
        clusterLeafs[i] = 1;
    }

    headValidNodes = &validNodes[0];
    headValidNodes->next = &validNodes[1];
    headValidNodes->n = -1;
    headValidNodes->node = -1;
    headValidNodes->prev = NULL;

    //build an initial link list
    ValidNode* curr = &validNodes[1];
    ValidNode* prev = headValidNodes;
    ValidNode* next = &validNodes[2];
    for (i = 0; i < leafsNum; i++)
    {
        curr->n = i;
        curr->node = i;
        curr->prev = prev;
        curr->next = next;
        prev = curr;
        curr = next;
        next++;
    }
    prev->next = NULL;

    //to generate the cluster tree
    int nodeIdx;	//the index of an internal node
    int firstNode = leafsNum;	//the index of the first internal node
    int lastNode = firstNode + leafsNum - 1;//the index of the last internal node

    for (nodeIdx = firstNode; nodeIdx < lastNode; nodeIdx++)
    {
        //find closest pair of clusters
        float minDist = 1.1f;
        miniPtr = headValidNodes;
        minjPtr = headValidNodes;

        for (ivalid = headValidNodes->next; ivalid != NULL;
                ivalid = ivalid->next)
        {
            mini = ivalid->n;

            for (jvalid = headValidNodes->next;
                    jvalid != NULL && jvalid->n < mini; jvalid = jvalid->next)
            {
                minj = jvalid->n;
                float dist = (*distMatrix)[mini][minj];
                if (dist < 0)
                {
                    cerr
                            << "ERROR: It is impossible to have distance value less than zero"
                            << endl;
                    dist = 0;
                }
                if (dist < minDist)
                {
                    minDist = dist;
                    miniPtr = ivalid;
                    minjPtr = jvalid;
                }
                //printf("dist %g mini %d minj %d\n", dist, ivalid->node, jvalid->node);
            }
        }
        //printf("**** mini %d minj %d minDist %g *****\n", miniPtr->node, minjPtr->node, minDist);
        //check the validity of miniPtr and minjPtr;
        if (miniPtr == headValidNodes || minjPtr == headValidNodes)
        {
            cerr << "OOPS: Error occurred while constructing the cluster tree\n"
                 << endl;
            exit(-1);
        }
        //computing branch length and join the two nodes
        float branchLength = minDist * 0.5f;
        this->connectNodes(&nodes[nodeIdx], nodeIdx, &nodes[miniPtr->node],
                           branchLength, &nodes[minjPtr->node], branchLength);
        clusterLeafs[nodeIdx] = clusterLeafs[miniPtr->node]
                                + clusterLeafs[minjPtr->node];

        //remove the valid node minjPtr from the list
        minjPtr->prev->next = minjPtr->next;
        if (minjPtr->next != NULL)
        {
            minjPtr->next->prev = minjPtr->prev;
        }
        minjPtr->prev = minjPtr->next = NULL;

        //compute the distance of each remaining valid node to the new node
        for (ivalid = headValidNodes->next; ivalid != NULL;
                ivalid = ivalid->next)
        {
            int idx = ivalid->n;

            float idist = (*distMatrix)[miniPtr->n][idx];
            float jdist = (*distMatrix)[minjPtr->n][idx];

            unsigned int isize = clusterLeafs[miniPtr->node];
            unsigned int jsize = clusterLeafs[minjPtr->node];
            if(second_step == 1){
                //cerr << "Step 2: using WPGMA." << endl;
                joins[idx] = (idist + jdist )/ 2;
            }
            else{
                //cerr << "Step 2: using UPGMA." << endl;
                joins[idx] = (idist * isize + jdist * jsize) / (isize + jsize);
            }
        }
        //update the distance to the new node
        miniPtr->node = nodeIdx;
        mini = miniPtr->n;
        for (jvalid = headValidNodes->next; jvalid != NULL;
                jvalid = jvalid->next)
        {
            minj = jvalid->n;

            float dist = joins[minj];
            (*distMatrix)[mini][minj] = dist;
            (*distMatrix)[minj][mini] = dist;
        }
    }
    //add a pseudo root to this unrooted NJ tree
    this->root = &nodes[lastNode - 1];

    delete[] validNodes;
    delete[] joins;
    delete[] clusterLeafs;
}
