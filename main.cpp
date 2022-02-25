#include <fstream>
#include <iostream>

#include "solver.hpp"

using namespace std;

int main(int, char **) {
    graph g;
    cin >> g;

    graph_search gs;
    bitvector<N> fvs = solve(g, g.active_vertices(), gs, 0);

    cout << 0 << endl;

    ofstream fs("data/solution/tmp");
    visit(fvs, [&](size_t v) {
        fs << v + 1 << endl;
    });

    return 0;
}

// int main(int, char **) {
//     graph g;
//     cin >> g;

//     graph_search gs;
//     bitvector<N> fvs{};
//     vector<bitvector<N>> SCC;
//     reduce_graph(g, fvs, g.active_vertices(), gs, SCC);

//     size_t n = popcount(g.active_vertices());
//     cout << n << endl;

//     if (n == 0) {
//         g.unfold_graph(0, fvs);
//         ofstream fs("data/solution/tmp");
//         visit(fvs, [&](size_t v) {
//             fs << v + 1 << endl;
//         });
//     } else if (n < 220) {
//         ofstream fs("scripts/plot_data");
//         g.print_edgelist(fs);
//     }

//     return 0;
// }
