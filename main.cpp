#include <fstream>
#include <iostream>

#include "reductions.hpp"

using namespace std;

int main(int, char **) {
    graph g;
    cin >> g;

    graph_search gs;
    vector<size_t> fvs;
    reduce_graph(g, fvs, gs);

    size_t n = popcount(g.active_vertices());
    cout << n << endl;

    if (n == 0) {
        ofstream fs("data/solution/tmp");
        for (size_t u : fvs)
            fs << u + 1 << endl;
    } else if (n < 100) {
        ofstream fs("scripts/visualize.py");
        g.print_edgelist(fs);
    }

    return 0;
}
