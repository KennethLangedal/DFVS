#include <fstream>
#include <iostream>

#include "reductions.hpp"

using namespace std;

int main(int, char **) {
    graph g;
    cin >> g;

    graph_search gs;
    bitvector<N> fvs{};
    reduce_graph(g, fvs, gs, true);

    size_t n = popcount(g.active_vertices());
    cout << n << endl;

    if (n == 0) {
        g.unfold_graph(0, fvs);
        ofstream fs("data/solution/tmp");
        visit(fvs, [&](size_t v) {
            fs << v + 1 << endl;
        });
    } else if (n < 160) {
        ofstream fs("scripts/plot_data");
        g.print_edgelist(fs);
    }

    return 0;
}
