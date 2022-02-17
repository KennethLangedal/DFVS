#pragma once
#include <array>
#include <functional>
#include <cstdint>
#include <cassert>

template <size_t N>
using bitvector = std::array<uint64_t, N>;

template <size_t N>
void visit(const bitvector<N>& v, std::function<void(size_t)>&& f) {
    uint64_t t, r;
    for (size_t i = 0; i < N; ++i) {
        t = v[i];
        while (t != 0) {
            f(i * 64 + __builtin_ctzll(t));
            t ^= t & -t;
        }
    }
}

template <size_t N>
void set(bitvector<N>& v, size_t i) {
    assert(i < 64 * N);
    v[i / 64] |= 1ull << (i & 63);
}

template <size_t N>
void reset(bitvector<N>& v, size_t i) {
    assert(i < 64 * N);
    v[i / 64] &= ~(1ull << (i & 63));
}

template <size_t N>
bool test(const bitvector<N>& v, size_t i) {
    assert(i < 64 * N);
    return (v[i / 64] & (1ull << (i & 63))) > 0;
}

template <size_t N>
bitvector<N>& operator|=(bitvector<N>& lhs, const bitvector<N>& rhs) {
    for (size_t i = 0; i < N; ++i) {
        lhs[i] |= rhs[i];
    }
    return lhs;
}

template <size_t N>
bitvector<N>& operator&=(bitvector<N>& lhs, const bitvector<N>& rhs) {
    for (size_t i = 0; i < N; ++i) {
        lhs[i] &= rhs[i];
    }
    return lhs;
}

template <size_t N>
bitvector<N> operator|(const bitvector<N>& lhs, const bitvector<N>& rhs) {
    bitvector<N> res;
    for (size_t i = 0; i < N; ++i) {
        res[i] = lhs[i] | rhs[i];
    }
    return res;
}

template <size_t N>
bitvector<N> operator&(const bitvector<N>& lhs, const bitvector<N>& rhs) {
    bitvector<N> res;
    for (size_t i = 0; i < N; ++i) {
        res[i] = lhs[i] & rhs[i];
    }
    return res;
}

template <size_t N>
bitvector<N> operator~(const bitvector<N>& rhs) {
    bitvector<N> res;
    for (size_t i = 0; i < N; ++i) {
        res[i] = ~rhs[i];
    }
    return res;
}

template <size_t N>
size_t popcount(const bitvector<N>& v) {
    size_t count = 0;
    for (size_t i = 0; i < N; ++i) {
        count += __builtin_popcountll(v[i]);
    }
    return count;
}