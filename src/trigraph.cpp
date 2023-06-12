#include "trigraph.hpp"

TriGraph::TriGraph(int n)
{
    edges = vector<unordered_set<int>>(n);
    red_edges = vector<unordered_set<int>>(n);
    this->backpointer = vector<int>(n);
    for (int i = 0; i < n; i++)
    {
        this->backpointer[i] = i;
        this->frontpointer[i] = i;
        active.insert(i);
    }
    vis = vector<bool>(n, false);
}

TriGraph::TriGraph(vector<int> V)
{
    edges = vector<unordered_set<int>>(V.size());
    red_edges = vector<unordered_set<int>>(V.size());
    backpointer = vector<int>(V.size());
    int n = V.size();
    for (size_t i = 0; i < V.size(); i++)
    {
        frontpointer[V[i]] = i;
        backpointer[i] = V[i];
        active.insert(i);
    }
    vis = vector<bool>(n, false);
}

int TriGraph::redDegree()
{
    int mx = 0;
    for (auto v: active) {
        mx = max(mx, (int)red_edges[v].size());
    }
    return mx;
}

bool TriGraph::makeRed(int u, int v)
{
    if (red_edges[u].find(v) != red_edges[u].end())
        return false;
    edges[u].erase(v);
    edges[v].erase(u);
    red_edges[u].insert(v);
    red_edges[v].insert(u);
    return true;
}

void TriGraph::revContract(const Reversecontract& rc)
{
    edges[rc.v] = unordered_set<int>();
    red_edges[rc.v] = unordered_set<int>();
    unordered_set<int> update_v;
    for (auto uv : rc.eadd)
    {
        edges[uv.first].insert(uv.second);
        edges[uv.second].insert(uv.first);
    }
    for (auto uv : rc.readd)
    {
        red_edges[uv.first].insert(uv.second);
        red_edges[uv.second].insert(uv.first);
    }
    for (auto uv : rc.redel)
    {
        red_edges[uv.first].erase(uv.second);
        red_edges[uv.second].erase(uv.first);
    }
    active.insert(rc.v);
}

vector<int> TriGraph::getV()
{
    return vector<int>(active.begin(), active.end());
}

Reversecontract TriGraph::contract(int u, int v)
{
    // TODO: this can probably be optimized
    vector<pair<int, int>> eadd;
    vector<pair<int, int>> readd;
    vector<pair<int, int>> redel;

    if (edges[u].find(v) != edges[u].end())
    {
        edges[u].erase(v);
        eadd.push_back({u, v});
    }
    if (red_edges[u].find(v) != red_edges[u].end())
    {
        red_edges[u].erase(v);
        readd.push_back({u, v});
    }

    unordered_set<int> neighborset;
    neighborset.insert(edges[u].begin(), edges[u].end());
    neighborset.insert(red_edges[u].begin(), red_edges[u].end());
    neighborset.insert(edges[v].begin(), edges[v].end());
    neighborset.insert(red_edges[v].begin(), red_edges[v].end());

    for (auto w : neighborset)
    {
        if (w == u || w == v)
            continue;
        if (!(edges[u].find(w) != edges[u].end() && edges[v].find(w) != edges[v].end()))
        {
            if (edges[u].find(w) != edges[u].end())
            {
                eadd.push_back({u, w});
            }
            if (makeRed(u, w))
            {
                redel.push_back({u, w});
            }
        }
    }
    
    // delete v
    for (auto w : edges[v])
    {
        edges[w].erase(v);
        eadd.push_back({v, w});
    }
    for (auto w : red_edges[v])
    {
        red_edges[w].erase(v);
        readd.push_back({v, w});
    }
    active.erase(v);
    int mxred = red_edges[u].size();
    for (auto w: neighborset) {
        if (w == v) continue;
        mxred = max(mxred, (int)red_edges[w].size());
    }
    return Reversecontract(v, eadd, readd, redel, mxred);
}

int TriGraph::numV()
{
    return active.size();
}

void TriGraph::revContractSeq(const vector<Reversecontract>& seq)
{
    for (int i = seq.size() - 1; i >= 0; i--)
    {
        revContract(seq[i]);
    }
}

vector<int> TriGraph::rneigh(const int v, const int r)
{
    queue<pair<int, int>> q;
    q.push({v, 0});
    vis[v] = true;
    vector<int> ans;
    while (q.size() > 0)
    {
        auto p = q.front();
        q.pop();
        ans.push_back(p.first);

        if (p.second < r)
        {
            for (auto w : this->edges[p.first])
            {
                if (!vis[w])
                {
                    vis[w] = true;
                    q.push({w, p.second + 1});
                }
            }
            for (auto w : this->red_edges[p.first])
            {
                if (!vis[w])
                {
                    vis[w] = true;
                    q.push({w, p.second + 1});
                }
            }
        }
    }
    for (auto v : ans)
        vis[v] = false;
    return ans;
}

bool TriGraph::isTwins(int u, int v)
{
    // TODO: can this be faster?
    if (u == v)
        return false;
    for (auto w: edges[v]) {
        if (w == u)continue;
        if (edges[u].find(w)==edges[u].end() && red_edges[u].find(w) == red_edges[u].end())
            return false;
    }
    for (auto w: red_edges[v]) {
        if (w == u)continue;
        if (red_edges[u].find(w) == red_edges[u].end())
            return false;
    }
    for (auto w: edges[u]) {
        if (w == v) continue;
        if (edges[v].find(w) == edges[v].end())
            return false;
    }
    return true;
}

vector<int> TriGraph::adjBlack(const int u)
{
    return vector<int>(edges[u].begin(), edges[u].end());
}

vector<int> TriGraph::adjRed(const int u) {
    return vector<int>(red_edges[u].begin(), red_edges[u].end());
}

vector<int> TriGraph::adjRedBlack(const int u) {
    vector<int> ans(edges[u].begin(), edges[u].end());
    ans.insert(ans.begin(), red_edges[u].begin(), red_edges[u].end());
    return ans;
}

unordered_set<int> TriGraph::symdiffblack(const int u, const int v) {
    auto Nu = adjBlack(u);
    auto Nv = adjBlack(v);
    unordered_set<int> Nus(Nu.begin(), Nu.end());
    unordered_set<int> Nvs(Nv.begin(), Nv.end());
    unordered_set<int> ans;
    for (auto w: Nu) {
        if (w != v)
            if (Nvs.find(w) == Nvs.end())ans.insert(w);
    }
    for (auto w: Nv) {
        if (w != u)
            if (Nus.find(w) == Nus.end())ans.insert(w);
    }
    return ans;
}

int TriGraph::symdiffsize(const int u, const int v)
{
    auto symdiff = symdiffblack(u,v);

    for (const auto w: red_edges[u]) {
        if (w != v)symdiff.insert(w);
    }
    for (const auto w: red_edges[v]) {
        if (w != u)symdiff.insert(w);
    }
    return symdiff.size();
}


int TriGraph::redblackDegree(const int u) {
    return edges[u].size() + red_edges[u].size();
}
int TriGraph::redDegree(const int u) {
    return red_edges[u].size();
}
int TriGraph::degree(const int u) {
    return edges[u].size();
}

int TriGraph::getNodeID(const int u) {
    return backpointer[u];
}

BitsetTriGraph* TriGraph::toBitsetTrigraph()
{
    auto V = this->getV();
    vector<int> Vid;
    for (auto v: V)Vid.push_back(backpointer[v]);
    BitsetTriGraph* g = new BitsetTriGraph(Vid);

    for (auto u: V) {
        for (auto v: edges[u]) {
            g->addEdge(g->frontpointer[backpointer[u]], g->frontpointer[backpointer[v]], false);
        }
        for (auto v: red_edges[u]) {
            g->addEdge(g->frontpointer[backpointer[u]], g->frontpointer[backpointer[v]], true);
        }
    }
    return g;
}


bool TriGraph::isActive(const int u)
{
    return active.find(u) != active.end();
}

void TriGraph::addEdge(const int u, const int v, const bool red) {
    if (red) {
        red_edges[u].insert(v);
        red_edges[v].insert(u);
    }
    else {
        edges[u].insert(v);
        edges[v].insert(u);
    }
}

bool TriGraph::contractPossible(const int u, const int v, const int sup) {
    auto redneighu = symdiffblack(u,v);
    for (auto w : red_edges[u]) {
        if (w != v)
            redneighu.insert(w);
    }
    for (auto w : red_edges[v]) {
        if (w != u)
            redneighu.insert(w);
    }
    if ((int)redneighu.size() >= sup) return false;
    for (auto w : symdiffblack(u,v)) {
        if ((red_edges[w].find(u) != red_edges[w].end()) || (red_edges[w].find(u) != red_edges[w].end())) continue;
        if ((int)red_edges[w].size() >= sup-1) return false;;
    }
    return true;
}

TriGraph TriGraph::extractComponent(const vector<int>& vertices) {
    
    vector<int> Vs;
    for (auto v: vertices) {
        Vs.push_back(backpointer[v]);
    }
    TriGraph TG(Vs);

    for (const auto v: vertices) {
        for (const auto w: edges[v]) {
            TG.edges[TG.frontpointer[backpointer[v]]].insert(TG.frontpointer[backpointer[w]]);
        }
        for (const auto w: red_edges[v]) {
            TG.red_edges[TG.frontpointer[backpointer[v]]].insert(TG.frontpointer[backpointer[w]]);
        }
    }
    return TG;
}

vector<TriGraph> TriGraph::connectedComponents() {
    
    unordered_set<int> vis;
    vector<TriGraph> ans;
    for (auto v: getV()) {
        
        if (vis.find(v) == vis.end()) {
            queue<int> q;
            q.push(v);
            vis.insert(v);
            vector<int> V;
            while(q.size() > 0) {                
                auto w = q.front();
                V.push_back(w);
                q.pop();
                for (auto x: edges[w]) {
                    if (vis.find(x) == vis.end()) {
                        vis.insert(x);
                        q.push(x);
                    }
                }
                for (auto x: red_edges[w]) {
                    if (vis.find(x) == vis.end()) {
                        vis.insert(x);
                        q.push(x);
                    }
                }
            }            
            ans.push_back(extractComponent(V));
        }
        
    }
    return ans;
}