#ifndef _MSA_GUIDE_TREE_H
#define _MSA_GUIDE_TREE_H
#include "MSADef.h"
#include "MSA.h"

#include "SafeVector.h"
#include "MultiSequence.h"
#include "ScoreType.h"
#include "ProbabilisticModel.h"
#include "SparseMatrix.h"

class MSA;
struct ValidNode
{
    ValidNode* prev;
    ValidNode* next;
    int n;				//the index in the distance matrix
    int node;			//the index in the tree node entries
};

struct TreeNode
{
    struct TreeNode *left;			//the pointer to its left child
    struct TreeNode *right;			//the pointer to its right child
    struct TreeNode *parent;		//the pointer to its parent
    int leftIdx;					//the index of the left child
    int rightIdx;					//the index of the right child
    int parentIdx;					//the index of its parent
    int idx;						//the index of itself
    float dist;						//the distance to its parent
    int leaf;						//whether it is a leaf node or not
    int order;			//the number of generations dating back to its ancestor
    int depth;						//the depth of the node
};
struct AlignmentOrder
{
    int nodeDepth;			//the depth of the internal node
    int leftOrder;			//the order number of the right child
    int rightOrder;			//the order number of the left child
    int* leftLeafs;			//the indices of leafs in the left subtree
    int leftNum;			//the number of leafs in the left subtree
    int* rightLeafs;			//the indices of leafs in the right subtree
    int rightNum;			//the number of leafs in the right substree
};

class MSAGuideTree
{
public:
    MSAGuideTree(MSA* msa, VVF& distMatrix, int numSeqs);
    virtual ~MSAGuideTree() = 0;	//abstract class

    //get the tree nodes
    TreeNode* getNodes();
    //get the leaf nodes
    TreeNode* getLeafs();
    //get the number of nodes;
    int getNodesNum();
    //get the number of leaf nodes
    int getLeafsNum();
    //get the root of the tree
    TreeNode* getRoot()
    {
        return this->root;
    }
    //get the alignment orders
    AlignmentOrder* getAlignOrders();
    int getAlignOrdersNum();
    //construct the alignment orders
    void createAlignmentOrders();

    //construct the guide tree
    virtual void create();

    virtual void create(int); //new add
    //calculate the sequence weights
    virtual void getSeqsWeights();

    /**********DEBUGING****************/
    //display the tree
    void displayTree();
    //display the alignment orders
    void displayAlignmentOrders();

protected:
    //join two nodes
    void connectNodes(TreeNode* parent, int parentIdx, TreeNode* leftChild,
                      float leftDist, TreeNode* rightChild, float rightDist);
    //release the alignment orders vector
    void releaseAlignmentOrders();
    //recursive implemenation of constructing the alignment orders
    int recursiveCreateAlignmentOrders(TreeNode* subRoot, int* subLeafs,
                                       int& subLeafsNum, int nodeDepth);

    //system configurations
    MSA* msa;
    VVF* distMatrix;
    int numSeqs;
    int* seqsWeights;

    //all the tree nodes
    TreeNode* nodes;
    int nodesNum;
    int nodesSize;
    //the root tree node
    TreeNode* root;
    //leaf node
    TreeNode* leafs;
    int leafsNum;

    //alignment order
    AlignmentOrder* alignOrders;
    int alignOrdersNum;
    int alignOrdersSize;
};
#endif

