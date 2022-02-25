#include "solver.hpp"
#include <optional>

size_t find_dense_branch_node(const graph &g, const bitvector<N> &nodes) {
    size_t best = N * 64, best_cnt;
    visit(nodes, [&](size_t u) {
        size_t u_cnt = g.degree(u) + g.in_degree(u);
        if (u_cnt > 0 && (best == N * 64 || u_cnt > best_cnt)) {
            best = u;
            best_cnt = u_cnt;
        }
    });

    return best;
}

#include <fstream>
bitvector<N> solve(graph &g, const bitvector<N> &nodes, graph_search &gs, size_t d) {
    size_t t0 = g.get_timestamp();
    bitvector<N> fvs{};
    std::vector<bitvector<N>> SCC;
    reduce_graph(g, fvs, nodes, gs, SCC);
    bitvector<N> remaining = nodes & g.active_vertices();
    if (popcount(remaining) == 0) {
        g.unfold_graph(t0, fvs);
        return fvs;
    }
    if (d > 20) {
        std::cout << popcount(g.active_vertices()) << std::endl;
        if (popcount(g.active_vertices()) < 200) {
            std::ofstream fs("scripts/plot_data");
            g.print_edgelist(fs);
        }
        exit(0);
    }
    // for each CC
    for (auto &&CC : SCC) {
        if (popcount(CC) == 0) {
            continue;
        }
        size_t u = find_dense_branch_node(g, CC);
        size_t t1 = g.get_timestamp();
        // include u
        bitvector<N> fvs_inc_u{};
        add_to_fvs(g, fvs_inc_u, gs, u);
        fvs_inc_u |= solve(g, CC, gs, d + 1);
        g.unfold_graph(t1, fvs_inc_u);

        // exclude u
        bitvector<N> fvs_ex_u{};
        exclude_from_fvs(g, fvs_ex_u, gs, u);
        fvs_ex_u |= solve(g, CC, gs, d + 1);
        g.unfold_graph(t1, fvs_ex_u);

        if (popcount(fvs_inc_u) < popcount(fvs_ex_u)) {
            fvs |= fvs_inc_u;
        } else {
            fvs |= fvs_ex_u;
        }
    }

    g.unfold_graph(t0, fvs);

    return fvs;
}