#include "bitvector.hpp"
#include <cassert>
#include <numeric>

void bitvector_iterator::_advance() {
    while (_i < _N && _t == 0) {
        ++_i;
        ++_it;
        if (_i < _N)
            _t = *_it;
    }
}

bitvector_iterator::bitvector_iterator(std::vector<uint64_t>::const_iterator it, size_t N, bool start)
    : _it(it), _t(0), _i(0), _N(N) {
    if (start) {
        _t = *_it;
        _advance();
    } else {
        _i = N;
    }
}

uint32_t bitvector_iterator::operator*() const {
    return _i * 64 + __builtin_ctzll(_t);
}

bitvector_iterator &bitvector_iterator::operator++() {
    _t ^= _t & -_t;
    _advance();
    return *this;
}

bitvector_iterator bitvector_iterator::operator++(int) {
    bitvector_iterator tmp = *this;
    ++(*this);
    return tmp;
}

bool bitvector_iterator::operator==(const bitvector_iterator &v) {
    return _it == v._it && _t == v._t && _i == v._i && _N == v._N;
}

bool bitvector_iterator::operator!=(const bitvector_iterator &v) {
    return _it != v._it || _t != v._t || _i != v._i || _N != v._N;
}

bitvector_iterator &bitvector_iterator::advance(size_t dist) {
    size_t pop;
    while (_i < _N && dist > 0) {
        pop = __builtin_popcountll(_t);
        if (pop <= dist) {
            _t = 0;
            ++_i;
            ++_it;
            dist -= pop;
            if (_i < _N)
                _t = *_it;
        } else {
            while (dist > 0) {
                _t ^= _t & -_t;
                dist--;
            }
        }
    }
    return *this;
}

bitvector::bitvector(size_t N)
    : _data((N + 63) / 64, 0), _N(N), _count(0), _changed(false) {
}

size_t bitvector::popcount() const {
    if (_changed) {
        _count = std::accumulate(std::begin(_data), std::end(_data), 0ull, [&](size_t v, uint64_t d) { return v + __builtin_popcountll(d); });
        _changed = false;
    }
    return _count;
}

size_t bitvector::size() const {
    return _N;
}

size_t bitvector::intersection_size(const bitvector &v) const {
    assert(_N == v._N);
    size_t res = 0;
    for (size_t i = 0; i < _data.size(); ++i) {
        res += __builtin_popcountll(_data[i] & v._data[i]);
    }
    return res;
}

bool bitvector::get(size_t i) const {
    assert(i < _N);
    return _data[i / 64] & (1ull << (i & 63));
}

void bitvector::set(size_t i) {
    assert(i < _N);
    _data[i / 64] |= 1ull << (i & 63);
    _changed = true;
}

void bitvector::reset(size_t i) {
    assert(i < _N);
    _data[i / 64] &= ~(1ull << (i & 63));
    _changed = true;
}

void bitvector::clear() {
    for (size_t i = 0; i < _data.size(); ++i) {
        _data[i] = 0;
    }
    _changed = true;
}

void bitvector::fill() {
    for (size_t i = 0; i < _data.size(); ++i) {
        _data[i] = ~0ull;
    }
    if (_N & 63)
        _data.back() = (1ull << (_N & 63)) - 1;
    _changed = true;
}

bool bitvector::subset_eq(const bitvector &v) const {
    assert(_N == v._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        if (_data[i] & ~v._data[i])
            return false;
    }
    return true;
}

bool bitvector::operator==(const bitvector &v) const {
    assert(_N == v._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        if (_data[i] != v._data[i])
            return false;
    }
    return true;
}

bool bitvector::operator!=(const bitvector &v) const {
    assert(_N == v._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        if (_data[i] != v._data[i])
            return true;
    }
    return false;
}

bitvector &bitvector::operator&=(const bitvector &v) {
    assert(_N == v._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        _data[i] &= v._data[i];
    }
    _changed = true;
    return *this;
}

bitvector &bitvector::operator|=(const bitvector &v) {
    assert(_N == v._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        _data[i] |= v._data[i];
    }
    _changed = true;
    return *this;
}

bitvector &bitvector::operator^=(const bitvector &v) {
    assert(_N == v._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        _data[i] ^= v._data[i];
    }
    _changed = true;
    return *this;
}

bitvector &bitvector::set_and(const bitvector &a, const bitvector &b) {
    assert(_N == a._N && _N == b._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        _data[i] = a._data[i] & b._data[i];
    }
    _changed = true;
    return *this;
}

bitvector &bitvector::set_or(const bitvector &a, const bitvector &b) {
    assert(_N == a._N && _N == b._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        _data[i] = a._data[i] | b._data[i];
    }
    _changed = true;
    return *this;
}

bitvector &bitvector::set_xor(const bitvector &a, const bitvector &b) {
    assert(_N == a._N && _N == b._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        _data[i] = a._data[i] ^ b._data[i];
    }
    _changed = true;
    return *this;
}

bitvector &bitvector::set_not(const bitvector &a) {
    assert(_N == a._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        _data[i] = ~a._data[i];
    }
    if (_N & 63)
        _data.back() &= (1ull << (_N & 63)) - 1;
    _changed = true;
    return *this;
}

bitvector &bitvector::set_and_not(const bitvector &a, const bitvector &b) {
    assert(_N == a._N && _N == b._N);
    for (size_t i = 0; i < _data.size(); ++i) {
        _data[i] = a._data[i] & ~b._data[i];
    }
    _changed = true;
    return *this;
}

bitvector_iterator bitvector::begin() const {
    return bitvector_iterator(std::begin(_data), _data.size(), true);
}

bitvector_iterator bitvector::end() const {
    return bitvector_iterator(std::end(_data), _data.size(), false);
}