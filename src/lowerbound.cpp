#include "lowerbound.hpp"

LowerBound::LowerBound(TriGraph &TG, milliseconds timeout) : TG(TG), timeout(timeout)
{
    lastLB = 0;
    start_t = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch());
    auto V = TG.getV();
    auto maxv = *max_element(V.begin(), V.end());

    ignore = vector<unordered_set<int>>(maxv + 1);
}

int LowerBound::getLB(TriGraph* currTG, milliseconds t_left, vector<pair<int,int>>* lastseq)
{

    Solver solver(*currTG, TG.numV(), lastLB, this->ignore);
    solver.WITH_TIMEOUT = true;
    solver.t_until = duration_cast<milliseconds>(
                         system_clock::now().time_since_epoch()) +
                     t_left;
    solver.VERBOSE = false;
    vector<Reversecontract> rcsec;
    int twwub = TG.numV();
    if ((int)lastseq->size() == currTG->numV()-2) {
        int tww = 0;
        for (auto p : *lastseq)
        {
            rcsec.push_back(solver.WorkingGraph->contract(solver.WorkingGraph->frontpointer[p.first], solver.WorkingGraph->frontpointer[p.second]));
            tww = max(tww, rcsec[rcsec.size() - 1].maxred);
        }
        auto V = solver.WorkingGraph->getV();
        rcsec.push_back(solver.WorkingGraph->contract(V[0], V[1]));
        tww = max(tww, rcsec[rcsec.size() - 1].maxred);
        solver.UB = tww;
        solver.ansseq = vector<pair<int, int>>(solver.seq.begin(), solver.seq.end());
        solver.seq.clear();
        solver.WorkingGraph->revContractSeq(rcsec);
        twwub = min(twwub, tww);
    }
    

    if (twwub > lastLB)
    {
        auto rcseq = solver.contractTwins(-1);
        solver.bruteForce(0);
        if (solver.TIMEOUT)
            return -1;
        solver.WorkingGraph->revContractSeq(rcseq);
    }
    lastseq->clear();
    for (auto p: solver.ansseq) {
        lastseq->push_back({p.first,p.second});
    }
    return solver.UB;
}

void LowerBound::insertV(TriGraph* currTG, int v)
{
    int vindex = currTG->edges.size();
    currTG->edges.push_back(unordered_set<int>());
    currTG->red_edges.push_back(unordered_set<int>());
    currTG->active.insert(vindex);
    currTG->backpointer.push_back(v);
    currTG->frontpointer[v] = vindex;
    currTG->vis.push_back(false);

    for (auto w : TG.edges[v])
    {
        if (currTG->frontpointer.find(w) != currTG->frontpointer.end())
        {
            currTG->edges[vindex].insert(currTG->frontpointer[w]);
            currTG->edges[currTG->frontpointer[w]].insert(vindex);
        }
    }
}

int LowerBound::iterativeLB2()
{
    if (TG.numV() <= 3)
    {
        return 0;
    }
    set<set<int>> setsseen;
    vector<unordered_set<int>*> vertsfree(TG.numV());
    vector<unordered_set<int>*> vertshave(TG.numV());
    vector<TriGraph*> currTGs(TG.numV());
    vector<vector<pair<int, int>>*> lastseqs(TG.numV());

    auto V = TG.getV();
    for (size_t i = 0; i < V.size(); i++)
    {
        vertsfree[i] = new unordered_set<int>(V.begin(), V.end());
        vertsfree[i]->erase(V[i]);
        vertshave[i] = new unordered_set<int>({V[i]}); // unordered_set<int>({V[i]});
        currTGs[i] = new TriGraph();
        lastseqs[i] = new vector<pair<int,int>>();
        insertV(currTGs[i], V[i]);
    }

    while (currTGs[0]->numV() < (int)V.size())
    {
        unordered_set<int> keepindices;
        for (size_t i = 0; i < vertsfree.size(); i++)
        {
            int nextv = -1;
            int maxval = -1;

            for (auto v : *vertsfree[i])
            {
                unordered_set<int> inverts;
                for (auto w : this->TG.edges[v])
                {
                    if (vertshave[i]->find(w) != vertshave[i]->end())
                        inverts.insert(w);
                }

                if (inverts.size())
                {
                    vector<int> symdiffs;
                    for (auto u : *vertshave[i])
                    {
                        vector<int> symdiff;
                        for (auto w : currTGs[i]->edges[currTGs[i]->frontpointer[u]])
                        {
                            if (inverts.find(currTGs[i]->backpointer[w]) == inverts.end())
                                symdiff.push_back(currTGs[i]->backpointer[w]);
                        }
                        for (auto w : inverts)
                        {
                            if (currTGs[i]->edges[currTGs[i]->frontpointer[u]].find(currTGs[i]->frontpointer[w]) == currTGs[i]->edges[currTGs[i]->frontpointer[u]].end())
                                symdiff.push_back(w);
                        }
                        symdiffs.push_back(symdiff.size());
                    }
                    if (*min_element(symdiffs.begin(), symdiffs.end()) > maxval)
                    {
                        maxval = *min_element(symdiffs.begin(), symdiffs.end());
                        nextv = v;
                    }
                }
            }

            vertshave[i]->insert(nextv);
            vertsfree[i]->erase(nextv);
            this->insertV(currTGs[i], nextv);
            milliseconds t_now = duration_cast<milliseconds>(
                system_clock::now().time_since_epoch());
            milliseconds t_left = timeout - (t_now - start_t);
            set<int> vertshaveset(vertshave[i]->begin(), vertshave[i]->end());
            if (setsseen.find(vertshaveset) == setsseen.end()) {
                // cout << "c Trying to compute lower bound for number of vertices: " << vertshave[i]->size() << endl;
            
                setsseen.insert(vertshaveset);
                if (t_left <= chrono::milliseconds(0))
                    goto GETMEOUTHERE;
                int ans = getLB(currTGs[i], t_left, lastseqs[i]);
                if (ans == -1)
                    goto GETMEOUTHERE;
                keepindices.insert(i);
                if (ans > lastLB)
                {
                    lastLB = ans;
                    cout << "c New lower bound: " << lastLB << endl;
                }
            }
            
        }
        vector<unordered_set<int>*> tempvertsfree;
        vector<unordered_set<int>*> tempvertshave;
        vector<TriGraph*> tempcurrTGs;
        vector<vector<pair<int,int>>*> templastseqs;
        for (size_t i = 0; i < vertsfree.size(); i++) {
            if (keepindices.find(i) != keepindices.end()) {
                tempvertsfree.push_back(vertsfree[i]);
                tempvertshave.push_back(vertshave[i]);
                tempcurrTGs.push_back(currTGs[i]);
                templastseqs.push_back(lastseqs[i]);
            }
            else {
                delete vertsfree[i];
                delete vertshave[i];
                delete currTGs[i];
                delete lastseqs[i];
            }
        }
        vertsfree = tempvertsfree;
        vertshave = tempvertshave;
        currTGs = tempcurrTGs;
        lastseqs = templastseqs;        
    }
    GETMEOUTHERE:
    for (auto p: vertsfree)delete p;
    for (auto p: vertshave)delete p;
    for (auto p: currTGs)delete p;
    for (auto p: lastseqs)delete p;
    return lastLB;
}

int LowerBound::iterativeLB()
{
    if (TG.numV() <= 3)
        return 0;
    TriGraph* currTG = new TriGraph();
    vector<int> V = TG.getV();

    unordered_set<int> vertsfree(V.begin(), V.end());
    unordered_set<int> vertshave;
    vector<pair<int,int>>* lastseq = new vector<pair<int,int>>;

    // select some element with maximum degree
    int maxdeg = -1;
    int nextv = -1;
    for (auto v : vertsfree)
    {
        if (TG.redblackDegree(v) > maxdeg)
        {
            nextv = v;
            maxdeg = TG.redblackDegree(v);
        }
    }
    vertshave.insert(nextv);
    vertsfree.erase(nextv);
    this->insertV(currTG, nextv);

    while (currTG->numV() < (int)V.size())
    {
        int nextv = -1;
        int maxval = -1;

        for (auto v : vertsfree)
        {
            unordered_set<int> inverts;
            for (auto w : this->TG.edges[v])
            {
                if (vertshave.find(w) != vertshave.end())
                    inverts.insert(w);
            }

            if (inverts.size())
            {
                vector<int> symdiffs;
                for (auto u : vertshave)
                {
                    vector<int> symdiff;
                    for (auto w : currTG->edges[currTG->frontpointer[u]])
                    {
                        if (inverts.find(currTG->backpointer[w]) == inverts.end())
                            symdiff.push_back(currTG->backpointer[w]);
                    }
                    for (auto w : inverts)
                    {
                        if (currTG->edges[currTG->frontpointer[u]].find(currTG->frontpointer[w]) == currTG->edges[currTG->frontpointer[u]].end())
                            symdiff.push_back(w);
                    }
                    symdiffs.push_back(symdiff.size());
                }
                if (*min_element(symdiffs.begin(), symdiffs.end()) > maxval)
                {
                    maxval = *min_element(symdiffs.begin(), symdiffs.end());
                    nextv = v;
                }
            }
        }

        vertshave.insert(nextv);
        vertsfree.erase(nextv);
        this->insertV(currTG, nextv);
        milliseconds t_now = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch());
        milliseconds t_left = timeout - (t_now - start_t);
        cout << "c Trying to compute lower bound for number of vertices: " << vertshave.size() << endl;
        if (t_left <= chrono::milliseconds(0))
            break;
        int ans = getLB(currTG, t_left, lastseq);
        if (ans == -1)
            break;
        if (ans > lastLB)
        {
            lastLB = ans;
            cout << "c New lower bound: " << lastLB << endl;
        }
    }
    delete currTG;
    delete lastseq;
    return lastLB;
}