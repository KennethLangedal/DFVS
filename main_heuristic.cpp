#include "dense_graph.hpp"
#include "sparse_graph.hpp"
#include "reductions.hpp"
#include "local_search.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <signal.h>
#include <sstream>
#include <unistd.h>

using namespace std;

volatile sig_atomic_t tle = 0;

void term(int signum) {
    tle = 1;
}

template <class graph>
void solve_heuristic(graph &g) {
    graph_search gs(g.size());
    bitvector res(g.size());
    reduction_engine<graph> re;
    reduce_graph(g, res, gs, re);

    if (g.active_vertices().popcount() == 0) {
        re.unfold_graph(g, 0, res);
#ifdef VERBOSE
        cout << res.popcount() << endl;
#else
        for (auto u : res) {
            cout << u + 1 << endl;
        }
#endif
        return;
    }

#ifdef VERBOSE
    cout << g.size() << " " << g.active_vertices().popcount() << " " << res.popcount() << endl;

    size_t total_degree = 0, total_pi_degree = 0;
    for (auto u : g.active_vertices()) {
        total_degree += g.out_degree(u);
        total_pi_degree += g.pi_degree(u);
    }
    cout << (double)total_degree / (double)g.active_vertices().popcount() << " " << ((double)total_pi_degree / 2.0) / (double)g.active_vertices().popcount() << endl;
#endif

    if (g.active_vertices().popcount() <= (1 << 16)) {
        dense_graph _g(g);
    } else {
        sparse_graph _g(g);
    }
}

int main(int, char **) {

    struct sigaction action;
    memset(&action, 0, sizeof(struct sigaction));
    action.sa_handler = term;
    sigaction(SIGTERM, &action, NULL);

    size_t N, E, t;
    cin >> N >> E >> t;

    if (N <= (1 << 16) || true) {
        dense_graph g;
        g.parse_graph(cin, N);
        solve_heuristic(g);
    } else {
        sparse_graph g;
        g.parse_graph(cin, N);
        solve_heuristic(g);
    }

    //     sparse_local_search ls(_g, 0.4, 0);

    //     double T = 0.25;

    //     uint32_t imp = _g.size();

    //     ls.greedy_one_zero_swaps(_g);

    //     while (!tle) {
    //         ls.set_temperature(T);
    //         ls.search(_g, std::max(_g.size(), 2500ul));

    //         if (imp > ls.get_best().popcount()) {
    //             imp = ls.get_best().popcount();
    //         } else {
    //             T *= 0.999;
    //             if (T < 0.1) {
    //                 ls.greedy_one_zero_swaps(_g);
    //                 ls.shuffle_solution(_g);
    //                 T = 0.25;
    //             }
    //         }

    // #ifdef VERBOSE
    //         std::cout << "\x1b[2K" << ls.get_best().popcount() + res.popcount() << " " << ls.get_current() + res.popcount() << " " << T << '\r' << std::flush;
    // #endif
    //     }

    //     for (auto u : ls.get_best()) {
    //         res.set(org_label[u]);
    //     }

    // #ifdef VERBOSE
    //     cout << "\x1b[2K" << res.popcount() << endl;
    // #else
    //     for (auto u : res) {
    //         cout << u + 1 << endl;
    //     }
    // #endif
    return 0;
}