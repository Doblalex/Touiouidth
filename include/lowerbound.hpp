#pragma once

#include "util.hpp"
#include "abstracttrigraph.hpp"
#include "bitsettrigraph.hpp"
#include "trigraph.hpp"
#include "solver.hpp"

class LowerBound {
public:
    TriGraph TG;
    milliseconds timeout;
    int lastLB;
    vector<unordered_set<int>> ignore;
    milliseconds start_t;
    LowerBound(TriGraph& TG, milliseconds timeout);
    int getLB(TriGraph* currTG, milliseconds t_left, vector<pair<int,int>>* lastseq);
    int iterativeLB();
    int iterativeLB2();
    void insertV(TriGraph* currTG, int v);
};