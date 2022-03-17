#include <fstream>
#include <iostream>

#include "solver.hpp"

using namespace std;

// int main(int, char **) {
//     graph g;
//     cin >> g;

//     graph_search gs;

//     auto res = analysis(g, gs);
//     cout << res.best_inc << "," << res.avg_inc << "," << res.worst_inc << endl;
//     cout << res.best_exc << "," << res.avg_exc << "," << res.worst_exc << endl;

//     return 0;
// }

// int main(int, char **) {
//     graph g;
//     cin >> g;

//     size_t lb_counter = 0;
//     graph_search gs(g.size());
//     bitvector fvs = solve(g, g.active_vertices(), gs, 0, g.size(), 0, lb_counter);

//     cout << 0 << endl;

//     ofstream fs("data/solution/tmp");
//     for (size_t v : fvs) {
//         fs << v + 1 << endl;
//     }

//     return 0;
// }

int main(int, char **) {
    graph g;
    cin >> g;

    graph_search gs(g.size());
    bitvector fvs(g.size());
    vector<bitvector> SCC;
    reduce_graph(g, fvs, g.active_vertices(), gs, SCC);

    size_t n = g.active_vertices().popcount();
    cout << n << endl;

    if (n == 0) {
        g.unfold_graph(0, fvs);
        ofstream fs("data/solution/tmp");
        for (size_t v : fvs) {
            fs << v + 1 << endl;
        }
    } else {
        ofstream fs("scripts/plot_data");
        g.print_edgelist(fs, g.active_vertices());
    }

    return 0;
}
