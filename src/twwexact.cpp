#include "util.hpp"
#include "bitsettrigraph.hpp"
#include "solver.hpp"
#include "lowerbound.hpp"



TriGraph readGraph()
{
    string line;
    TriGraph g;

    while (!cin.eof())
    {
        getline(cin, line);
        if (line.size() == 0)
        {
            continue;
        }
        if (line[0] == 'c')
        {
            continue;
        }
        else if (line[0] == 'p')
        {
            stringstream stream(line);
            string a, b;
            int n, m;
            stream >> a >> b >> n >> m;
            g = TriGraph(n);
        }
        else
        {
            stringstream stream(line);
            int u, v;
            stream >> u >> v;
            u--;
            v--;
            g.addEdge(u,v,false);
        }
    }
    cout << "c Read graph" << endl;
    return g;
}


milliseconds BUDGET = milliseconds(15*60*1000);

int main()
{

    milliseconds ms_start = duration_cast< milliseconds >(
        system_clock::now().time_since_epoch()
    );
    TriGraph graph = readGraph();
    // initialize here and use for all solvers, saves runtime to use vector over unordered_set
    vector<unordered_set<int>> ignore(graph.numV(), unordered_set<int>());

    vector<int> Vsleft;
    int LB = 0;
    vector<pair<int,int>> ansseq;
    
    
    for(auto TG: graph.connectedComponents()) {        
        Solver solver(TG, INT_MAX, LB, ignore); 

        
        solver.contractTwins(-1);

        if (solver.USE_BITSET) {
            for (auto p: solver.seq) {
                TG.contract(TG.frontpointer[p.first], TG.frontpointer[p.second]);
            }
        }
        LowerBound lowerbound(TG, BUDGET);
        int LBsubgraph = lowerbound.iterativeLB2();
        LB = max(LB, LBsubgraph);
        solver.LB = LB;
        
        solver.computeUpperBound(0);
        solver.upperBoundRandom(milliseconds(5*60*1000));
        //solver.UB = 1000;
        solver.bruteForce(0);
        LB = max(LB, solver.UB);
        
        ansseq.insert(ansseq.end(), solver.ansseq.begin(), solver.ansseq.end());
        // the remaining vertex is the first vertex in the last contraction
        if (ansseq.size())
            Vsleft.push_back(solver.ansseq[solver.ansseq.size()-1].first);
        else
            Vsleft.push_back(TG.getV()[0]);
    }

    for(size_t i = 1; i < Vsleft.size(); i++) {
        ansseq.push_back({Vsleft[0], Vsleft[i]});
    }

    cout<<"c computed solution with tww "<<LB<<endl;

    for (auto p : ansseq)
    {
        cout << p.first + 1 << " " << p.second + 1 << endl;
    }
}