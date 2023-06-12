#pragma once

#include "util.hpp"
#include "abstracttrigraph.hpp"
#include "bitsettrigraph.hpp"
#include "trigraph.hpp"
//#include "ModDecomp.h"

class Solver{
public:
    int UB;
    int LB;
    bool VERBOSE = true;
    bool USE_BITSET = false;
    bool TIMEOUT = false;
    vector<pair<int,int>> ansseq;
    vector<pair<int,int>> seq;
    // sometimes this is a bitsettrigraph as soon as we have <= 64 vertices
    ITriGraph* WorkingGraph;
    // This will always be a non-bitset trigraph
    TriGraph* TGTriGraph;
    vector<unordered_set<int>>& ignore;
    bool WITH_TIMEOUT = false;
    milliseconds t_until;
    Solver(TriGraph &TG, int UB, int LB, vector<unordered_set<int>>& ignore);
    ~Solver() {
        if (USE_BITSET) delete this->WorkingGraph;
    }
    vector<Reversecontract> contractTwins(const int justcontracted);
    bool cost(const int u, const int v);
    vector<pair<int, int>> ignoreRemOnConstract(const int u);
    void computeUpperBound(int LB);
    void bruteForce(int LB);
    void upperBoundDegeneracy();
    void upperBoundFast();
    void upperBoundFast2(); // not fast enough for heuristic track
    void upperBoundRandom(milliseconds t_budget);
    //vector<Reversecontract> contractModularDecomposition();
};