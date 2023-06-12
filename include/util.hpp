#pragma once

#include <bits/stdc++.h>
#include <random>
#include <tuple>
#include <string>
#define MP make_pair
#define PB push_back
#define st first
#define nd second
#define NOT_PLACED INT_MAX

using namespace std;
using namespace std::chrono;

static std::random_device rd; // random device engine, usually based on /dev/random on UNIX-like systems
// initialize Mersennes' twister using rd to generate the seed
static std::mt19937 rng{rd()};

template <typename TH>
void _dbg(const char *sdbg, TH h) { cerr << sdbg << "=" << h << "\n"; }
template <typename TH, typename... TA>
void _dbg(const char *sdbg, TH h, TA... t)
{
    while (*sdbg != ',')
        cerr << *sdbg++;
    cerr << "=" << h << ",";
    _dbg(sdbg + 1, t...);
}
#ifndef DEBUG
#define debug(...)
#else
#ifdef LOCAL
#define debug(...) _dbg(#__VA_ARGS__, __VA_ARGS__)
#define debugv(x)                 \
    {                             \
        {                         \
            cerr << #x << " = ";  \
            FORE(itt, (x))        \
            cerr << *itt << ", "; \
            cerr << "\n";         \
        }                         \
    }
#else
#define debug(...) (__VA_ARGS__)
#define debugv(x)
#define cerr \
    if (0)   \
    cout
#endif
#endif

#if defined(__GNUC__)

#define MY_LIB_IGNORE_DEPRECATED_BEGIN \
    _Pragma("GCC diagnostic push")     \
        _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")

#define MY_LIB_IGNORE_DEPRECATED_END \
    _Pragma("GCC diagnostic pop")

#elif defined(_MSC_VER)

#define MY_LIB_IGNORE_DEPRECATED_BEGIN \
    _Pragma("warning(push)")           \
        _Pragma("warning(disable : 4996)")

#define MY_LIB_IGNORE_DEPRECATED_END \
    _Pragma("warning(pop)")

#else

#define MY_LIB_IGNORE_DEPRECATED_BEGIN
#define MY_LIB_IGNORE_DEPRECATED_END

#endif

template <class T>
ostream &operator<<(ostream &out, vector<T> vec)
{
    out << "(";
    for (auto &v : vec)
        out << v << ", ";
    return out << ")";
}
template <class T>
ostream &operator<<(ostream &out, set<T> vec)
{
    out << "(";
    for (auto &v : vec)
        out << v << ", ";
    return out << ")";
}

class NotImplemented : public std::logic_error
{
public:
    NotImplemented() : std::logic_error("Function not yet implemented") { };
};

struct Reversecontract
{
    vector<pair<int, int>> eadd;
    vector<pair<int, int>> readd;
    vector<pair<int, int>> redel;
    int v;
    int maxred;

    Reversecontract(int v, vector<pair<int, int>> eadd,
    vector<pair<int, int>> readd,
    vector<pair<int, int>> redel, int maxred) {
        this->v = v;
        this->eadd = eadd;
        this->readd = readd;
        this->redel = redel;
        this->maxred = maxred;
    }
};

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}