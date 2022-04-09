#include <fstream>
#include <iostream>

// #define HEURISTIC
#define EXACT

#include "local_search.hpp"
#include "solver.hpp"

using namespace std;

#ifdef EXACT

int main(int, char **) {
    graph g;
    cin >> g;

    auto fvs = solve(g);

    cout << fvs.popcount() << endl;

    // for (size_t v : fvs) {
    //     cout << v + 1 << endl;
    // }

    return 0;
}

#endif

#ifdef HEURISTIC

int main(int, char **) {
    graph g;
    cin >> g;

    graph_search gs(g.size());
    bitvector dfvs(g.size());
    vector<bitvector> SCC;
    reduce_graph(g, dfvs, g.active_vertices(), gs, SCC);

    size_t largest = 0;

    for (auto &&c : SCC) {
        largest = max(largest, c.popcount());
    }

    cout << g.active_vertices().popcount() << " " << SCC.size() << " " << largest << endl;

    ofstream fs("scripts/plot_data");
    g.print_edgelist(fs, g.active_vertices());

    // cout << g.active_vertices().popcount() << " " << dfvs.popcount() << endl;

    // reducing_peeling(g, dfvs, gs, SCC);

    // cout << g.active_vertices().popcount() << " " << dfvs.popcount() << endl;

    // remove_redundant(g, dfvs);

    // cout << g.active_vertices().popcount() << " " << dfvs.popcount() << endl;

    // two_one_swap(g, dfvs);

    // cout << g.active_vertices().popcount() << " " << dfvs.popcount() << endl;
}

#endif