#ifndef _MSA_CLUSTER_TREE_H
#define _MSA_CLUSTER_TREE_H

#include "MSAGuideTree.h"

class MSAClusterTree: public MSAGuideTree
{
public:
    MSAClusterTree(MSA* msa, VVF& distMatrix, int numSeqs);
    ~MSAClusterTree();

    //construct the cluster tree
    void create();
    void create(int second_step);
private:
    //generate the cluster tree
    void generateClusterTree();
    void generateClusterTree(int second_step);
};
#endif
