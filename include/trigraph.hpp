#pragma once
#include "util.hpp"
#include "bitsettrigraph.hpp"
#include "abstracttrigraph.hpp"

// TODO, maybe frontpointer can be removed
class TriGraph : public ITriGraph{
private:
    bool makeRed(const int u, const int v);
public:
    vector<unordered_set<int>> edges;
    vector<unordered_set<int>> red_edges;
    unordered_set<int> active;
    vector<bool> vis;

    TriGraph(int n);
    TriGraph(vector<int> V);
    TriGraph() {}
    ~TriGraph() = default;


    int redDegree();
    vector<int> getV();
    int numV();
    Reversecontract contract(const int a, const int b);
    void revContract(const Reversecontract& rc);
    void revContractSeq(const vector<Reversecontract>& rc);
    void addEdge(const int u, const int v, const bool red);

    // this is non-symmetric!
    bool isTwins(const int u, const int v);
    vector<int> ullToV(unsigned long long x);
    vector<int> adjRedBlack(const int u);
    vector<int> adjRed(const int u);
    vector<int> adjBlack(const int u);
    bool isActive(const int u);
    unordered_set<int> symdiffblack(const int u, const int v);
    bool contractPossible(const int u, const int v, const int sup);
    vector<int> rneigh(const int v, const int r);
    int symdiffsize(const int u, const int v);
    int redblackDegree(const int u);
    int degree(const int u);
    int redDegree(const int u);
    int getNodeID(const int u);
    BitsetTriGraph* toBitsetTrigraph();
    TriGraph extractComponent(const vector<int>& vertices);
    vector<TriGraph> connectedComponents();
};