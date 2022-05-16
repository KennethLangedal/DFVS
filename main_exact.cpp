#include "solver.hpp"

using namespace std;

int main(int, char **) {
    sparse_graph g;
    cin >> g;

    auto fvs = solve(g);

#ifdef VERBOSE
    cout << "Solution size " << fvs.popcount() << endl;
#else
    for (size_t v : fvs) {
        cout << v + 1 << endl;
    }
#endif

    return 0;
}