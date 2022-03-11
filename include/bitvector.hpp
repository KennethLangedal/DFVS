#pragma once
#include <cstdint>
#include <vector>

class bitvector_iterator {
private:
    std::vector<uint64_t>::const_iterator _it;
    uint64_t _t, _i, _N;

    void _advance();

public:
    bitvector_iterator(std::vector<uint64_t>::const_iterator it, size_t N, bool start);

    size_t operator*() const;

    bitvector_iterator &operator++();
    bitvector_iterator operator++(int);

    bool operator==(const bitvector_iterator &v);
    bool operator!=(const bitvector_iterator &v);
};

class bitvector {
private:
    std::vector<uint64_t> _data;
    size_t _N;

    mutable size_t _count;
    mutable bool _changed;

public:
    bitvector() = default;
    bitvector(size_t N);

    size_t popcount() const;
    size_t intersection_size(const bitvector &v) const;

    bool get(size_t i) const;

    void set(size_t i);
    void reset(size_t i);
    void clear();
    void fill();

    bool subset_eq(const bitvector &v) const;

    bool operator==(const bitvector &v) const;
    bool operator!=(const bitvector &v) const;

    bitvector &operator&=(const bitvector &v);
    bitvector &operator|=(const bitvector &v);
    bitvector &operator^=(const bitvector &v);

    bitvector &set_and(const bitvector &a, const bitvector &b);
    bitvector &set_or(const bitvector &a, const bitvector &b);
    bitvector &set_xor(const bitvector &a, const bitvector &b);
    bitvector &set_and_not(const bitvector &a, const bitvector &b);

    bitvector_iterator begin() const;
    bitvector_iterator end() const;
};
