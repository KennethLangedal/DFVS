#include "solver.hpp"
#include <optional>

std::optional<size_t> find_dense_branch_node(const graph &g, const bitvector<N> &nodes) {
    size_t best = N * 64, best_cnt;
    visit(nodes, [&](size_t u) {
        size_t u_cnt = popcount(g.out(u) & g.in(u));
        if (u_cnt > 0 && (best == N * 64 || u_cnt > best_cnt)) {
            best = u;
            best_cnt = u_cnt;
        }
    });
    if (best == N * 64)
        return std::nullopt;

    return best;
}

bitvector<N> solve(graph &g, const bitvector<N> &nodes, graph_search &gs, size_t d) {
    if (d > 15) {
        std::cout << "Max search depht exceeded" << std::endl;
        exit(0);
    }
    size_t t0 = g.get_timestamp();
    bitvector<N> fvs{};
    std::vector<bitvector<N>> SCC;
    reduce_graph(g, fvs, nodes, gs, SCC);
    bitvector<N> remaining = nodes & g.active_vertices();
    if (popcount(remaining) == 0) {
        g.unfold_graph(t0, fvs);
        return fvs;
    }
    // for each CC
    for (auto &&CC : SCC) {
        auto u = find_dense_branch_node(g, CC);
        if (!u.has_value()) {
            std::cout << "Failed to find branch" << std::endl;
            exit(0);
        }
        size_t t1 = g.get_timestamp();
        // include u
        bitvector<N> fvs_inc_u{};
        add_to_fvs(g, fvs_inc_u, gs, *u);
        fvs_inc_u |= solve(g, CC, gs, d + 1);
        g.unfold_graph(t1, fvs_inc_u);

        // exclude u
        bitvector<N> fvs_ex_u{};
        visit(g.in(*u) & g.out(*u), [&](size_t v){
            add_to_fvs(g, fvs_ex_u, gs, v);
        });
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