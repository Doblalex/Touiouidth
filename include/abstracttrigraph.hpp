#pragma once

#include "util.hpp"

class ITriGraph {
public:
    vector<int> backpointer;
    unordered_map<int, int> frontpointer;
    virtual ~ITriGraph() = default;
    virtual int redDegree() = 0;
    virtual vector<int> getV() = 0;
    virtual int numV() = 0;
    virtual Reversecontract contract(const int a, const int b) = 0;
    virtual void revContract(const Reversecontract& rc) = 0;
    virtual void revContractSeq(const vector<Reversecontract>& rc) = 0;
    virtual void addEdge(const int u, const int v, const bool red) = 0;

    // this is non-symmetric!
    virtual bool isTwins(const int u, const int v) = 0;
    virtual vector<int> adjRedBlack(const int u) = 0;
    virtual vector<int> adjRed(const int u) = 0;
    virtual vector<int> adjBlack(const int u) = 0;
    virtual bool isActive(const int u) = 0;
    virtual unordered_set<int> symdiffblack(const int u, const int v) = 0;
    virtual bool contractPossible(const int u, const int v, const int sup) = 0;
    virtual vector<int> rneigh(const int v, const int r) = 0;
    virtual int symdiffsize(const int u, const int v) = 0;
    virtual int degree(const int u) = 0;
    virtual int redblackDegree(const int u) = 0;
    virtual int redDegree(const int u) = 0;
    virtual int getNodeID(const int u) = 0;
};