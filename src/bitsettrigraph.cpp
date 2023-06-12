#include "bitsettrigraph.hpp"

BitsetTriGraph::BitsetTriGraph(int n)
{
    assert(n <= 64);
    edges = vector<unsigned long long>(n, 0);
    red_edges = vector<unsigned long long>(n, 0);
    backpointer = vector<int>(n);
    bit = vector<unsigned long long>(n);
    unsigned long long b = 1ULL;
    active = 0;
    for (char i = 0; i < n; i++)
    {
        backpointer[i] = i;
        frontpointer[i] = i;
        bit[i] = b;
        active |= bit[i];
        b <<= 1ULL;
    }
}

BitsetTriGraph::BitsetTriGraph(vector<int> V)
{
    int n = V.size();
    assert(n <= 64);
    edges = vector<unsigned long long>(n, 0);
    red_edges = vector<unsigned long long>(n, 0);
    backpointer = vector<int>(n);
    bit = vector<unsigned long long>(n);
    unsigned long long b = 1ULL;
    active = 0;
    for (char i = 0; i < n; i++)
    {
        backpointer[i] = V[i];
        frontpointer[V[i]] = i;
        bit[i] = b;
        active |= bit[i];
        b <<= 1ULL;
    }
}

int BitsetTriGraph::redDegree()
{
    int mx = 0;
    for (auto v: getV()) {
        mx = max(mx, __builtin_popcountll(red_edges[v]));
    }
    return mx;
}

vector<int> BitsetTriGraph::getV()
{
    unsigned long long temp = active;
    int w;
    vector<int> activevec;
    while (temp)
    {
        w = __builtin_ctzll(temp);
        activevec.push_back(w);
        temp &= (~bit[w]);
    }
    return activevec;
}

Reversecontract BitsetTriGraph::contract(const int u, const int v)
{
    assert(u != v);
    vector<pair<int,int>> resetaddred;
    vector<pair<int,int>> resetaddblack;
    vector<pair<int,int>> resetremovered;

    // black intersection
    unsigned long long nonred = this->edges[u] & this->edges[v];

    // everything that will be red neighbors of a now
    unsigned long long red1 = (this->edges[u] | this->edges[v] | this->red_edges[u] | this->red_edges[v]) & ~(nonred | this->bit[u] | this->bit[v]);
    unsigned long long red2 = red1;

    int w;
    while (red1) {
        w = __builtin_ctzll(red1);
        if (!(bit[w] & this->red_edges[u])) {
            // w is not in red edges of u but will be now
            resetremovered.push_back({u,w});
            if (bit[w] & this->edges[u]) {
                // but w was a bleck neighbor of u
                resetaddblack.push_back({u,w});
            }
        }
        // remove black u-v (if exists) and add red u-v
        this->edges[u] &= ~(bit[w]);
        this->edges[w] &= ~(bit[u]);
        this->red_edges[u] |= bit[w];
        this->red_edges[w] |= bit[u];

        red1 &= (~bit[w]);
    }
    unsigned long long blackv = this->edges[v];
    unsigned long long redv = this->red_edges[v];


    // remove black v-w
    while (blackv) {
        w = __builtin_ctzll(blackv);
        resetaddblack.push_back({v,w});
        this->edges[w] &= ~(bit[v]);
        blackv &= (~bit[w]);
    }
    // remove red v-w
    while (redv) {
        w = __builtin_ctzll(redv);
        resetaddred.push_back({v,w});
        this->red_edges[w] &= ~(bit[v]);
        redv &= (~bit[w]);
    }
    active &= (~bit[v]);
    // compute maximum red degree among all affected vertices
    int maxred = __builtin_popcountll(red_edges[u]);
    while(red2) {
        w = __builtin_ctzll(red2);
        maxred = max(maxred, __builtin_popcountll(red_edges[w]));
        red2 &= (~bit[w]);
    }    
    return Reversecontract(v, resetaddblack, resetaddred, resetremovered, maxred);
}

void BitsetTriGraph::revContract(const Reversecontract& rc)
{
    for (auto p: rc.eadd) {
        this->edges[p.first] |= this->bit[p.second];
        this->edges[p.second] |= this->bit[p.first];
    }
    for (auto p: rc.readd) {
        this->red_edges[p.first] |= this->bit[p.second];
        this->red_edges[p.second] |= this->bit[p.first];
    }
    for (auto p: rc.redel) {
        this->red_edges[p.first] &= ~this->bit[p.second];
        this->red_edges[p.second] &= ~this->bit[p.first];
    }
    this->active |= bit[rc.v];
}

void BitsetTriGraph::revContractSeq(const vector<Reversecontract>& seq)
{
    for (int i = seq.size() - 1; i >= 0; i--)
    {
        revContract(seq[i]);
    }
}

void BitsetTriGraph::addEdge(const int u, const int v, const bool red) {
    if (red) {
        this->red_edges[u] |= bit[v];
        this->red_edges[v] |= bit[u];
    }
    else {
        this->edges[u] |= bit[v];
        this->edges[v] |= bit[u];
    }
}

bool BitsetTriGraph::isTwins(const int u, const int v)
{
    if (u == v)
        return false;
    const unsigned long long mask = ~(bit[u] | bit[v]);
    
    const unsigned long long blacku = edges[u] & mask;
    const unsigned long long redu = red_edges[u] & mask;
    const unsigned long long redblacku = redu | blacku;
    const unsigned long long blackv = edges[v] & mask;
    const unsigned long long redv = red_edges[v] & mask;
    // (1) if black neighbor of v then black or red neighbor of u
    // (2) if red neighbor of v then red neighbor of u
    // (3) if black neighbor of u then black neighbor of v
    return ((blackv | redblacku) == redblacku) && ((redv | redu) == redu) && ((blacku | blackv) == blackv);
}

int BitsetTriGraph::numV() {
    return __builtin_popcountll(active);
}

vector<int> BitsetTriGraph::ullToV(unsigned long long x){
    int w;
    vector<int> vec;
    while (x)
    {
        w = __builtin_ctzll(x);
        vec.push_back(w);
        x &= (~bit[w]);
    }
    return vec;
}

vector<int> BitsetTriGraph::adjRedBlack(const int u) {
    return ullToV(this->edges[u] | this->red_edges[u]);
}

vector<int> BitsetTriGraph::adjRed(const int u) {
    return ullToV(this->red_edges[u]);
}

vector<int> BitsetTriGraph::adjBlack(const int u) {
    return ullToV(this->edges[u]);
}

bool BitsetTriGraph::isActive(const int u) {
    return bit[u] & active;
}

unordered_set<int> BitsetTriGraph::symdiffblack(const int u, const int v) {
    auto Nu = adjBlack(u);
    auto Nv = adjBlack(v);

    unordered_set<int> symDifference;
    for (auto w: Nu) {
        if (w != v && (edges[v] & bit[w]) == 0)symDifference.insert(w);
    }
    for (auto w: Nv) {
        if (w != u && (edges[u] & bit[w]) == 0)symDifference.insert(w);
    }
    return symDifference;
}

bool BitsetTriGraph::contractPossible(const int u, const int v, const int sup) {
    auto symdiff = ((edges[u]&(~edges[v]))|(edges[v]&(~edges[u])))&(~(bit[u]|bit[v]));
    if (__builtin_popcountll((red_edges[u]|red_edges[v]|symdiff) & (~(bit[u]|bit[v]))) >= sup) {
        return false;
    }
    for (auto w : symdiffblack(u,v)) {
        if ((red_edges[w] & bit[u]) || (red_edges[w] & bit[v])) continue;
        if (__builtin_popcountll(red_edges[w]) >= sup-1) {
            return false;
        }
    }
    return true;
}

vector<int> BitsetTriGraph::rneigh(const int v, const int r)
{
    queue<pair<int, int>> q;
    q.push({v, 0});
    unsigned long long vis = 0;
    vis |= bit[v];
    vector<int> ans;
    while (q.size() > 0)
    {
        auto p = q.front();
        q.pop();
        ans.push_back(p.first);

        if (p.second < r)
        {
            auto ws = this->edges[p.first] | this->red_edges[p.first];
            int w;
            while (ws)
            {
                w = __builtin_ctzll(ws);
                ws &= (~bit[w]);
                if (!(vis & bit[w]))
                {
                    vis |= bit[w];
                    q.push({w, p.second + 1});
                }
            }
        }
    }
    return ans;
}

int BitsetTriGraph::symdiffsize(const int u, const int v)
{
    unsigned long long Nu = edges[u];
    unsigned long long Nv = edges[v];

    auto symdiff = Nu ^ Nv;
    symdiff |= (red_edges[u] | red_edges[v]);
    symdiff &= (~bit[u]);
    symdiff &= (~bit[v]);

    return __builtin_popcountll(symdiff);
}

int BitsetTriGraph::redblackDegree(const int u) {
    return __builtin_popcountll(edges[u]|red_edges[u]);
}

int BitsetTriGraph::redDegree(const int u) {
    return __builtin_popcountll(red_edges[u]);
}

int BitsetTriGraph::degree(const int u) {
    return __builtin_popcountll(edges[u]);
}

int BitsetTriGraph::getNodeID(const int u) {
    return backpointer[u];
}