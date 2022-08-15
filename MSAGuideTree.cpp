////////////////////////////////////////////////////////////////
// MSAGuideTree.cpp
//
// Utilities for tree data structure
/////////////////////////////////////////////////////////////////

#include "MSAGuideTree.h"
#include "MSA.h"
MSAGuideTree::MSAGuideTree(MSA* msa, VVF& distances, int numSeqs)
{
    int i;
    TreeNode* node;
    //system configuration
    this->msa = msa;
    this->distMatrix = &distances;
    this->numSeqs = numSeqs;
    this->seqsWeights = msa->getSeqsWeights();

    //tree structure
    this->nodesSize = this->numSeqs * 2 + 1;
    this->nodes = new TreeNode[this->nodesSize];
    if (!this->nodes)
    {
        cerr << "TreeNodes memory allocation failed" << endl;
        exit(-1);
    }
    //initialize all the tree nodes
    this->leafs = this->nodes;
    this->leafsNum = this->numSeqs;
    this->nodesNum = 2 * this->leafsNum - 1;
    for (i = 0; i < this->nodesSize; i++)
    {
        node = &nodes[i];
        node->left = 0;
        node->right = 0;
        node->parent = 0;
        node->leftIdx = -1;
        node->rightIdx = -1;
        node->parentIdx = -1;
        node->idx = -1;
        node->dist = 0;
        node->leaf = NODE;		//setted to be NODE, by default
        node->order = 0;
        node->depth = 0;
    }
    //initialize the leaf nodes
    for (i = 0; i < this->leafsNum; i++)
    {
        node = &this->leafs[i];
        node->idx = i;
        node->leaf = LEAF;
    }
}
MSAGuideTree::~MSAGuideTree()
{
    //release tree nodes
    delete[] this->nodes;

    //release alignment orders
    releaseAlignmentOrders();

}
//get the tree nodes
TreeNode* MSAGuideTree::getNodes()
{
    return nodes;
}
//get the leaf nodes
TreeNode* MSAGuideTree::getLeafs()
{
    return leafs;
}
//get the number of nodes;
int MSAGuideTree::getNodesNum()
{
    return nodesNum;
}
//get the number of leaf nodes
int MSAGuideTree::getLeafsNum()
{
    return leafsNum;
}
//get the alignment orders
AlignmentOrder* MSAGuideTree::getAlignOrders()
{
    return alignOrders;
}
int MSAGuideTree::getAlignOrdersNum()
{
    return alignOrdersNum;
}
/****************************************************
 create the evolutionary relationship
 ****************************************************/
void MSAGuideTree::connectNodes(TreeNode* parent, int parentIdx,
                                TreeNode* leftChild, float leftDist, TreeNode* rightChild,
                                float rightDist)
{
    //save the parents index for each child
    leftChild->parent = parent;
    leftChild->parentIdx = parentIdx;
    rightChild->parent = parent;
    rightChild->parentIdx = parentIdx;

    //save the branch lengths (i.e. distance) from each child to its parent
    leftChild->dist = leftDist;
    rightChild->dist = rightDist;

    //save the indices of itself and its children for this new tree node
    parent->idx = parentIdx;
    parent->left = leftChild;
    parent->leftIdx = leftChild->idx;
    parent->right = rightChild;
    parent->rightIdx = rightChild->idx;
}
/*****************************************
 compute the alignment order of the phylogentic tree
 *****************************************/
void MSAGuideTree::createAlignmentOrders()
{
    int i;

    AlignmentOrder* order;
    //allocate memory space for alignment orders vector
    this->alignOrdersNum = 0;//for alignment orders, it starts from 1 instead of 0
    this->alignOrdersSize = numSeqs;//the number of internal nodes of the phylogentic tree + 1
    this->alignOrders = new AlignmentOrder[this->alignOrdersSize];
    if (!this->alignOrders)
    {
        cerr << "OOPS: Alignment orders memory allocation failed" << endl;
        exit(-1);
    }
    //initialize the alignment orders vector
    for (i = 0; i < this->alignOrdersSize; i++)
    {
        order = &this->alignOrders[i];
        order->leftOrder = 0;
        order->rightOrder = 0;
        order->leftLeafs = 0;
        order->leftNum = 0;
        order->rightLeafs = 0;
        order->rightNum = 0;
    }
    //starting out constructing the alignment orders
    int subLeafsNum;
    int nodeDepth = 1;
    int subOrder = recursiveCreateAlignmentOrders(this->root, 0, subLeafsNum,
                   nodeDepth);

    //check whether the function works well
    if (subLeafsNum != numSeqs || this->alignOrdersNum != subOrder)
    {
        fprintf(stderr,
                "The alignment orders constructed were wrong (subLeafsNum %d, alignOrdersNum %d, subOrder %d)\n",
                subLeafsNum, alignOrdersNum, subOrder);
    }

}
int MSAGuideTree::recursiveCreateAlignmentOrders(TreeNode* subRoot,
        int* subLeafs, int& subLeafsNum, int nodeDepth)
{
    int leftNum, rightNum;
    int leftOrder, rightOrder;
    int* leftLeafs, *rightLeafs;

    if (subRoot->leaf == LEAF)
    {
        subLeafs[0] = subRoot->idx;
        subLeafsNum = 1;

        return 0;			//if it is a leaf, return the index 0
    }
    leftOrder = rightOrder = 0;
    leftNum = rightNum = 0;
    leftLeafs = new int[numSeqs];
    rightLeafs = new int[numSeqs];

    //check the left subtree
    if (subRoot->left)
    {
        //recursively tranverse the left subtree
        leftOrder = recursiveCreateAlignmentOrders(subRoot->left, leftLeafs,
                    leftNum, nodeDepth + 1);
    }
    //check the right subtree
    if (subRoot->right)
    {
        rightOrder = recursiveCreateAlignmentOrders(subRoot->right, rightLeafs,
                     rightNum, nodeDepth + 1);
    }
    //save the leafs in the left and right subtrees of the current subtree
    if (this->alignOrdersNum > this->alignOrdersSize)
    {
        fprintf(stderr, "the alignment order function works bad\n");
        \
        exit(-1);
    }

    AlignmentOrder* order = &this->alignOrders[++this->alignOrdersNum];
    order->nodeDepth = nodeDepth;
    order->leftOrder = leftOrder;
    order->rightOrder = rightOrder;
    order->leftNum = leftNum;
    order->rightNum = rightNum;
    order->leftLeafs = new int[order->leftNum];
    order->rightLeafs = new int[order->rightNum];
    if (!order->leftLeafs || !order->rightLeafs)
    {
        fprintf(stderr,
                "memory allocation failed while recursively constructing alignment orders\n");
        exit(-1);
    }
    memcpy(order->leftLeafs, leftLeafs, order->leftNum * sizeof(int));
    memcpy(order->rightLeafs, rightLeafs, order->rightNum * sizeof(int));

    delete[] leftLeafs;
    delete[] rightLeafs;

    //for the root of the tree, subLeafs buffer is set to 0
    if (subLeafs)
    {
        //copy the results to the parent tree node
        memcpy(subLeafs, order->leftLeafs, order->leftNum * sizeof(int));
        memcpy(subLeafs + order->leftNum, order->rightLeafs,
               order->rightNum * sizeof(int));
    }
    //compute the total number of leafs in this subtree
    subLeafsNum = order->leftNum + order->rightNum;

    return this->alignOrdersNum;//return the index of itself, starting from 1, instead of 0
}
void MSAGuideTree::releaseAlignmentOrders()
{
    if (!this->alignOrders)
    {
        return;
    }
    for (int i = 0; i < this->alignOrdersNum; i++)
    {
        AlignmentOrder* order = &this->alignOrders[i];
        if (order->leftLeafs)
        {
            delete[] order->leftLeafs;
        }
        if (order->rightLeafs)
        {
            delete[] order->rightLeafs;
        }
    }
    delete[] alignOrders;
}
/********************************
 display the alignment orders
 ********************************/
void MSAGuideTree::displayAlignmentOrders()
{
    int i, j;
    AlignmentOrder* order;
    fprintf(stderr, "************DISPLAY ALIGNMENT ORDER***************\n");
    for (i = 1; i <= this->alignOrdersNum; i++)
    {
        order = &this->alignOrders[i];

        fprintf(stderr, "GROUP (%d depth %d):\n---LEFT ORDER: %d\n", i,
                order->nodeDepth, order->leftOrder);
        fprintf(stderr, "---LEFT: ");
        for (j = 0; j < order->leftNum; j++)
        {
            fprintf(stderr, "%d ", order->leftLeafs[j]);
        }

        fprintf(stderr, "\n---RIGHT ORDER: %d\n", order->rightOrder);
        fprintf(stderr, "\n---RIGHT: ");
        for (j = 0; j < order->rightNum; j++)
        {
            fprintf(stderr, "%d ", order->rightLeafs[j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "*******************************************\n");
}
/*********************************
 display the tree
 *********************************/
void MSAGuideTree::displayTree()
{
    fprintf(stderr, "**************DISPLAY TREE*********************\n");
    for (int i = 0; i < nodesNum; i++)
    {
        TreeNode* node = &nodes[i];

        fprintf(stderr,
                "%d(%p): left(%p) %d, right(%p) %d, parent(%p) %d, dist %f\n",
                (node == &nodes[node->idx]) ? node->idx : -2, node, node->left,
                (!node->left || node->left == &nodes[node->leftIdx]) ?
                node->leftIdx : -2, node->right,
                (!node->right || node->right == &nodes[node->rightIdx]) ?
                node->rightIdx : -2, node->parent,
                (!node->parent || node->parent == &nodes[node->parentIdx]) ?
                node->parentIdx : -2, node->dist);
    }
    fprintf(stderr, "*******************************************\n");
}
/*********************************
 compute the sequence weights
 *********************************/
void MSAGuideTree::getSeqsWeights()
{
    int i;
    TreeNode* curr;

    //compute the order of each node, which represents the number of leaf nodes in the substree rooting from it.
    for (i = 0; i < leafsNum; i++)
    {
        //for each leaf nodes
        curr = &this->leafs[i];
        while (curr != 0)
        {
            curr->order++;

            curr = curr->parent;
        }
    }
    //compute the weight of each sequence, which corresponds to a leaf node
    for (i = 0; i < numSeqs; i++)
    {
        //compute the weight of each sequence
        float weights = 0;
        curr = &this->leafs[i];
        while (curr->parent != 0)
        {
            weights += curr->dist / curr->order;
            curr = curr->parent;
            //printf("order:%d weights: %f\n", curr->order, weights);
        }
        //save the weight of this sequence
        seqsWeights[i] = (int) (100 * weights);
        //printf("%d\n", seqsWeights[i]);
    }
    //normalize the weights
    int wsum = 0;
    for (i = 0; i < numSeqs; i++)
    {
        wsum += seqsWeights[i];
    }
    if (wsum == 0)
    {
        //in this case, every sequence is assumed to have an identical weight
        for (i = 0; i < numSeqs; i++)
        {
            seqsWeights[i] = 1;
        }
        wsum = numSeqs;
    }
    //printf("wsum:%d \n", wsum);
    for (i = 0; i < numSeqs; i++)
    {
        seqsWeights[i] = (seqsWeights[i] * INT_MULTIPLY) / wsum;
        if (seqsWeights[i] < 1)
        {
            seqsWeights[i] = 1;
        }
        //printf("%d \n", seqsWeights[i]);
    }
}
void MSAGuideTree::create()
{
    //do nothing
}

void MSAGuideTree::create(int variance)
{
    //do nothing
}
