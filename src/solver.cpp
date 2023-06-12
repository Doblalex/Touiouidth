#include "solver.hpp"

Solver::Solver(TriGraph &TG, int UB, int LB, vector<unordered_set<int>> &ignore) : ignore(ignore)
{
    this->TGTriGraph = &TG;
    if (TGTriGraph->numV() > 64)
    {
        this->WorkingGraph = &TG;
        USE_BITSET = false;
    }
    else
    {
        this->WorkingGraph = TGTriGraph->toBitsetTrigraph();
        USE_BITSET = true;
    }
    this->LB = LB;
    this->UB = UB;
    this->seq = vector<pair<int, int>>();
    this->ansseq = vector<pair<int, int>>();
}

vector<Reversecontract> Solver::contractTwins(const int justcontracted)
{
    // TODO If no red edges then this is possible in linear time!
    vector<int> trytwins;
    if (justcontracted != -1)
    {
        trytwins = this->WorkingGraph->adjRedBlack(justcontracted);
        trytwins.push_back(justcontracted);
    }
    bool contracted = true;
    vector<Reversecontract> seq;
    while (contracted)
    {
        contracted = false;
        vector<int> V;
        if (justcontracted == -1)
        {
            V = this->WorkingGraph->getV();
        }
        else
        {
            V = trytwins;
        }
        for (auto u : V)
        {
            if (!this->WorkingGraph->isActive(u))
                continue;
            ;
            for (auto v : V)
            {
                if (!this->WorkingGraph->isActive(v) || !this->WorkingGraph->isActive(u))
                    continue;
                if (this->WorkingGraph->isTwins(u, v))
                {
                    this->seq.push_back({this->WorkingGraph->backpointer[u], this->WorkingGraph->backpointer[v]});
                    seq.push_back(this->WorkingGraph->contract(u, v));
                    contracted = true;
                }
            }
        }
    }
    return seq;
}

vector<pair<int, int>> Solver::ignoreRemOnConstract(const int u)
{
    vector<pair<int, int>> ignoresremoved;
    auto update = this->WorkingGraph->rneigh(u, 2);

    for (auto x : update)
    {
        for (auto y : vector<int>(this->ignore[this->WorkingGraph->backpointer[x]].begin(), this->ignore[this->WorkingGraph->backpointer[x]].end()))
        {
            this->ignore[this->WorkingGraph->backpointer[x]].erase(y);
            this->ignore[y].erase(this->WorkingGraph->backpointer[x]);
            ignoresremoved.push_back({this->WorkingGraph->backpointer[x], y});
        }
    }
    return ignoresremoved;
}

void Solver::computeUpperBound(int LB)
{
    vector<Reversecontract> seq;

    while (this->WorkingGraph->numV() > 1)
    {
        auto V = this->WorkingGraph->getV();
        int uc = -1;
        int vc = -1;
        pair<int, int> mindiff = {INT_MAX, INT_MAX};

        for (auto u : V)
        {
            auto N2u = this->WorkingGraph->rneigh(u, 2);
            for (auto v : N2u)
            {
                if (u >= v)
                    continue;
                int compval = this->WorkingGraph->symdiffsize(u, v);
                pair<int, int> comp = {compval, this->WorkingGraph->redblackDegree(u) + this->WorkingGraph->redblackDegree(v)};
                if (comp < mindiff)
                {
                    mindiff = comp;
                    uc = u;
                    vc = v;
                }
            }
        }
        if (vc != -1)
        {
            auto rc = this->WorkingGraph->contract(uc, vc);
            LB = max(LB, rc.maxred);
            this->seq.push_back({this->WorkingGraph->backpointer[uc], this->WorkingGraph->backpointer[vc]});
            seq.push_back(rc);
            auto seqq = this->contractTwins(uc);
            seq.insert(seq.end(), seqq.begin(), seqq.end());
        }
        else
        {
            auto rc = this->WorkingGraph->contract(V[0], V[1]);
            this->seq.push_back({this->WorkingGraph->backpointer[V[0]], this->WorkingGraph->backpointer[V[1]]});
            LB = max(LB, rc.maxred);
            seq.push_back(rc);
            auto seqq = this->contractTwins(V[0]);
            seq.insert(seq.end(), seqq.begin(), seqq.end());
        }
    }
    if (LB < this->UB)
    {
        if (this->VERBOSE)
            cout << "c solution computed by upper bound " << LB << endl;
        this->UB = LB;
        this->ansseq = vector<pair<int, int>>(this->seq);
    }
    for (size_t i = 0; i < seq.size(); i++)
        this->seq.pop_back();
    this->WorkingGraph->revContractSeq(seq);
}

void Solver::bruteForce(int LB)
{
    if (max(LB, this->LB) >= this->UB)
        return;

    if (WITH_TIMEOUT && duration_cast<milliseconds>(
                            system_clock::now().time_since_epoch()) > t_until)
    {
        TIMEOUT = true;
        return;
    }
    bool switched_modes = false;
    if (!USE_BITSET && this->WorkingGraph->numV() <= 64)
    {
        this->WorkingGraph = this->TGTriGraph->toBitsetTrigraph();
        USE_BITSET = true;
        switched_modes = true;
    }

    auto V = this->WorkingGraph->getV();
    if (V.size() == 1)
    {
        if (LB < this->UB)
        {
            if (this->VERBOSE)
                cout << "c New UB: " << LB << endl;
            this->ansseq = vector<pair<int, int>>(this->seq);
            this->UB = LB;
        }
    }

    vector<vector<pair<int, int>>> trypairs(this->UB + 1);
    vector<pair<int, int>> ignoresadded;
    for (size_t i = 0; i < V.size(); i++)
    {
        int u = V[i];
        for (size_t j = i + 1; j < V.size(); j++)
        {
            int v = V[j];
            if (this->WorkingGraph->contractPossible(u, v, this->UB) && this->ignore[this->WorkingGraph->backpointer[u]].find(this->WorkingGraph->backpointer[v]) == this->ignore[this->WorkingGraph->backpointer[u]].end())
            {
                trypairs[this->WorkingGraph->symdiffsize(u, v)].push_back({u, v});
            }
        }
    }
    /*
    shuffle(V.begin(), V.end(), rng);
    vector<pair<int, int>> ignoresadded;
    for (size_t i = 0; i < V.size(); i++)
    {
        int u = V[i];
        // for (auto v: this->WorkingGraph->rneigh(u, 1)) { if (u >= v) continue;
        for (size_t j = i + 1; j < V.size(); j++)
        {
            int v = V[j]; */
    for (auto vec : trypairs)
    {
        for (auto pp : vec)
        {
            int u = pp.first;
            int v = pp.second;
            if (TIMEOUT)
                return;
            auto rc = this->WorkingGraph->contract(u, v);
            auto ignoresremoved = this->ignoreRemOnConstract(u);
            this->seq.push_back({this->WorkingGraph->backpointer[u], this->WorkingGraph->backpointer[v]});
            auto rcseqtwins = this->contractTwins(u);
            this->bruteForce(max(LB, rc.maxred));
            auto nContracted = 1 + rcseqtwins.size();
            while (nContracted)
            {
                this->seq.pop_back();
                nContracted--;
            }
            this->WorkingGraph->revContractSeq(rcseqtwins);
            this->WorkingGraph->revContract(rc);
            for (const auto &p : ignoresremoved)
            {
                this->ignore[p.first].insert(p.second);
                this->ignore[p.second].insert(p.first);
            }
            this->ignore[this->WorkingGraph->backpointer[u]].insert(this->WorkingGraph->backpointer[v]);
            this->ignore[this->WorkingGraph->backpointer[v]].insert(this->WorkingGraph->backpointer[u]);
            ignoresadded.push_back({this->WorkingGraph->backpointer[u], this->WorkingGraph->backpointer[v]});

            if (this->UB <= max(LB, this->LB))
            {
                // the branch found the same as the lower bound
                goto GETMEOUT;
            }
            if (TIMEOUT)
                goto GETMEOUT;
        }
    }
GETMEOUT:
    for (auto p : ignoresadded)
    {
        this->ignore[p.first].erase(p.second);
        this->ignore[p.second].erase(p.first);
    }
    if (switched_modes)
    {
        // switch back to not use bitset graph
        this->USE_BITSET = false;
        delete this->WorkingGraph;
        this->WorkingGraph = this->TGTriGraph;
    }
}

void Solver::upperBoundDegeneracy()
{
    vector<Reversecontract> rcseq;

    set<pair<int, int>> backdegrees;
    unordered_map<int, int> valinbackdegrees;
    vector<int> ordering;
    unordered_set<int> inordering;

    for (auto v : this->WorkingGraph->getV())
    {
        int deg = this->WorkingGraph->redblackDegree(v);
        backdegrees.insert({deg, v});
        valinbackdegrees[v] = deg;
    }

    while ((int)ordering.size() < this->WorkingGraph->numV())
    {
        pair<int, int> p = *backdegrees.begin();
        backdegrees.erase(p);

        ordering.push_back(p.second);
        inordering.insert(p.second);
        for (auto v : this->WorkingGraph->adjBlack(p.second))
        {
            if (inordering.find(v) == inordering.end())
            {
                backdegrees.erase({valinbackdegrees[v], v});
                valinbackdegrees[v]--;
                backdegrees.insert({valinbackdegrees[v], v});
            }
        }
    }
    reverse(ordering.begin(), ordering.end());
    int tww = 0;
    for (size_t i = 1; i < ordering.size(); i++)
    {
        auto rc = this->WorkingGraph->contract(ordering[0], ordering[i]);
        rcseq.push_back(rc);
        this->seq.push_back({ordering[0], ordering[i]});
        tww = max(tww, rc.maxred);
    }

    if (tww < this->UB)
    {
        if (this->VERBOSE)
            cout << "c solution computed by upper bound " << tww << endl;
        this->UB = tww;
        this->ansseq = vector<pair<int, int>>(this->seq);
    }
    for (size_t i = 0; i < rcseq.size(); i++)
        this->seq.pop_back();
    this->WorkingGraph->revContractSeq(rcseq);
}

void Solver::upperBoundFast()
{
    set<pair<pair<int, int>, int>> priority;
    int tww = 0;
    for (auto v : this->WorkingGraph->getV())
    {
        priority.insert({{this->WorkingGraph->redblackDegree(v), 0}, v});
    }

    while (this->WorkingGraph->numV() > 1)
    {
        auto p1 = *priority.begin();
        priority.erase(p1);
        auto p2 = *priority.begin();
        priority.erase(p2);

        auto rc = this->WorkingGraph->contract(p1.second, p2.second);
        auto v = p1.second;
        tww = max(tww, rc.maxred);
        priority.insert({{this->WorkingGraph->redblackDegree(v), this->WorkingGraph->redDegree(v)}, v});
    }
    this->UB = tww;
}

void Solver::upperBoundFast2()
{
    set<pair<pair<int, int>, int>> priority;
    unordered_map<int, pair<int, int>> inpriority;
    int tww = 0;
    for (auto v : this->WorkingGraph->getV())
    {
        this->WorkingGraph->redblackDegree(v);
        priority.insert({{this->WorkingGraph->redblackDegree(v), this->WorkingGraph->degree(v)}, v});
        inpriority[v] = {this->WorkingGraph->redblackDegree(v), this->WorkingGraph->degree(v)};
    }

    while (priority.size() > 1)
    {
        cout << tww << " " << priority.size() << endl;
        auto p1 = *priority.begin();
        priority.erase(p1);

        int v = p1.second;
        auto Nv = this->WorkingGraph->adjRedBlack(v);

        if (Nv.size())
        {
            int minsymdiff = INT_MAX;
            int wm;
            for (auto w : Nv)
            {
                int symdiff = this->WorkingGraph->symdiffsize(v, w);
                if (symdiff < minsymdiff)
                {
                    minsymdiff = symdiff;
                    wm = w;
                }
            }
            priority.erase({inpriority[wm], wm});
            auto rc = this->WorkingGraph->contract(v, wm);
            tww = max(rc.maxred, tww);
            pair<int, int> val = {this->WorkingGraph->redblackDegree(v), this->WorkingGraph->degree(v)};
            inpriority[v] = val;
            priority.insert({val, v});
        }
    }
    this->UB = tww;
}

void Solver::upperBoundRandom(milliseconds t_budget)
{
    auto t_start = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch());
    while (duration_cast<milliseconds>(system_clock::now().time_since_epoch())-t_start<=t_budget && this->UB > this->LB) {
        vector<Reversecontract> rcseq;

        int tww = 0;
        while (this->WorkingGraph->numV() > 1) {
            if (duration_cast<milliseconds>(system_clock::now().time_since_epoch())-t_start>t_budget)break;
            // vector<pair<int,int>> contractions_possible;
            unordered_map<int,vector<pair<int,int>>> contractions_possible;
            auto V = this->WorkingGraph->getV();
            for (size_t i = 0; i < V.size(); i++) {
                for (size_t j = i+1; j < V.size(); j++) {
                    int u = V[i];
                    int v = V[j];
                    if (this->WorkingGraph->contractPossible(u, v, this->UB)) {
                        contractions_possible[this->WorkingGraph->symdiffsize(u,v)].push_back({u,v});
                        // contractions_possible.push_back({u,v});
                    }
                }
            }
            if (contractions_possible.size() == 0)break;
            vector<int> symdiffsizes;
            for (const auto& [k, v]: contractions_possible)symdiffsizes.push_back(k);
            sort(symdiffsizes.begin(), symdiffsizes.end());
            // select something small here
            std::uniform_int_distribution<int> uni(0,symdiffsizes.size()-1);
            int indexselected = min(uni(rng), min(uni(rng), uni(rng)));
            auto p = *select_randomly(contractions_possible[symdiffsizes[indexselected]].begin(), contractions_possible[symdiffsizes[indexselected]].end(), rng);
            // auto p = *select_randomly(contractions_possible.begin(), contractions_possible.end(), rng)
            auto rc = this->WorkingGraph->contract(p.first, p.second);
            this->seq.push_back({this->WorkingGraph->backpointer[p.first], this->WorkingGraph->backpointer[p.second]});
            tww = max(tww, rc.maxred);
            rcseq.push_back(rc);
            auto rcseqtwins = this->contractTwins(p.first);
            rcseq.insert(rcseq.end(), rcseqtwins.begin(), rcseqtwins.end());
        }

        if (this->WorkingGraph->numV() == 1) {
            this->UB = tww;
            cout<<"c New UB: "<<this->UB<<endl;
            this->ansseq = vector<pair<int, int>>(this->seq);
        }
        this->WorkingGraph->revContractSeq(rcseq);
        for (auto rc: rcseq)this->seq.pop_back();
    }
}

/*vector<Reversecontract> Solver::contractModularDecomposition() {
    Graph g;
    int cnt = 0;
    for (auto v: this->WorkingGraph->getV()) {
        for (auto w: this->WorkingGraph->adjRedBlack(v)) {
            if (v < w) {
                boost::add_edge(v,w,cnt++,g);
            }
        }
    }
    ModDecomp md;
    Tree md_tree = md.decompose(g);
    printTree(md_tree);
    vector<Reversecontract> seq;
    return seq;
}*/