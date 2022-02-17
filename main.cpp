#include <iostream>

#include "reductions.hpp"

using namespace std;

int main(int, char**) {
    graph g;
    cin >> g;

    graph_search gs;
    vector<size_t> fvs;
    reduce_graph(g, fvs, gs);
    

    cout << popcount(g.active_vertices()) << endl;
    g.print_edgelist(cout);

    return 0;
}
