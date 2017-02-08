#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <fstream>
#include <limits>
#include <chrono>
#include <vector>
#include <queue>
#include <array>
#include <cmath>
#include <map>

// Uncovering some classes from std
using std::vector;
using std::string;
using std::array;
using std::pair;
using std::cerr;
using std::cout;
using std::map;
using std::set;

#define val(x) #x "=", x

#define fprint(os, ...) ;
#define fprintln(os, ...) ;

// #define fprint(os, ...) debug_print(os, __VA_ARGS__)
// #define fprintln(os, ...)             \
// 	{                                 \
// 		debug_print(os, __VA_ARGS__); \
// 		os << '\n';                   \
// 	}

#define print(...) debug_print(std::cout, __VA_ARGS__)
#define println(...)                         \
	{                                        \
		debug_print(std::cout, __VA_ARGS__); \
		std::cout << '\n';                   \
	}

void debug_print(std::ostream &os) {
}
template <typename T, typename... Args>
void debug_print(std::ostream &os, T t, Args... args) {
	os << t << " ";
	debug_print(os, args...);
}

template <class T>
size_t getIndex(const std::vector<T> &v, const T &elem) {
	auto ptr = lower_bound(v.begin(), v.end(), elem);
	if (ptr == v.end())
		throw std::runtime_error("World is coming to an end with elem");
	if (elem != *ptr)
		throw std::runtime_error("World is coming to an end with no matching element");
	return ptr - v.begin();
}

template <class T>
void sortAndRemoveDuplicate(std::vector<T> &v) {
	std::sort(v.begin(), v.end());			// vector may have repeated elements like 1 1 2 2 3 3 3 4 4 5 5 6 7
	auto last = unique(v.begin(), v.end()); // vector now holds {1 2 3 4 5 6 7 x x x x x x}, where 'x' is indeterminate
	v.erase(last, v.end());
}

template <typename T, typename T2>
size_t addToMap(std::map<T, size_t> &mapping, std::vector<T2> &vec, const T &object) {
	if (mapping.find(object) == mapping.end()) {
		mapping.insert({object, vec.size()});
		vec.emplace_back();
	}
	return mapping[object];
}

template <typename T, size_t N>
void sortThis(std::array<T, N> &t) {
	std::sort(t.begin(), t.end());
}

template <typename T>
using myPair = std::array<T, 2>;

template <typename T>
using myTriple = std::array<T, 3>;

class myEdge {
	size_t data[2];

  public:
	myEdge(size_t a, size_t b) : data{a, b} {
		if (data[0] > data[1]) std::swap(data[0], data[1]);
	}
	size_t &operator[](int index) {
		if (index < 0 || index > 1)
			throw std::out_of_range("Index tried to access=" + index);
		return data[index];
	}
	const size_t &operator[](int index) const {
		if (index < 0 || index > 1)
			throw std::out_of_range("Index tried to access=" + index);
		return data[index];
	}
};

bool operator<(const myEdge &lhs, const myEdge &rhs) {
	if (lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]))
		return true;
	return false;
}

bool operator!=(const myEdge &lhs, const myEdge &rhs) {
	return lhs[0] != rhs[0] || lhs[1] != rhs[1];
}

std::ostream &operator<<(std::ostream &os, const myEdge &e) {
	return os << "(" << e[0] << " " << e[1] << ")";
}

class myFace {
	size_t data[3];

  public:
	myFace(size_t a, size_t b, size_t c) : data{a, b, c} {
		if (data[0] > data[1]) std::swap(data[0], data[1]);
		if (data[1] > data[2]) std::swap(data[1], data[2]);
		if (data[0] > data[1]) std::swap(data[0], data[1]);
	}
	size_t &operator[](int index) {
		if (index < 0 || index > 2)
			throw std::out_of_range("Index tried to access=" + index);
		return data[index];
	}
	const size_t &operator[](int index) const {
		if (index < 0 || index > 2)
			throw std::out_of_range("Index tried to access=" + index);
		return data[index];
	}
};

bool operator<(const myFace &lhs, const myFace &rhs) {
	if (lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]) || (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] < rhs[2]))
		return true;
	return false;
}

bool operator!=(const myFace &lhs, const myFace &rhs) {
	return lhs[0] != rhs[0] || lhs[1] != rhs[1] || lhs[2] != rhs[2];
}

std::ostream &operator<<(std::ostream &os, const myFace &f) {
	return os << "(" << f[0] << " " << f[1] << " " << f[2] << ")";
}

#endif // UTIL_H
