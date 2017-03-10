#ifndef UTIL_H
#define UTIL_H

#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <queue>

// Uncovering some classes from std
using std::vector;
using std::string;
using std::array;
using std::pair;
using std::cerr;
using std::cout;
using std::map;
using std::set;
using std::stack;
using std::make_tuple;
using std::make_pair;

#define val(x) #x "=", x

#define fprint(os, ...) ;
#define fprintln(os, ...) ;
// #define print(os, ...) ;
// #define println(os, ...) ;

// #define fprint(os, ...) debug_print(os, __VA_ARGS__)
// #define fprintln(os, ...)             \
// 	{                                 \
// 		debug_print(os, __VA_ARGS__); \
// 		os << '\n';                   \
// 	}

#define print(...) debug_print(__VA_ARGS__)
#define println(...)              \
	{                             \
		debug_print(__VA_ARGS__); \
		std::cout << std::endl;   \
	}

void debug_print() {
}
template <typename T, typename... Args>
void debug_print(T t, Args... args) {
	cout << t << " ";
	debug_print(args...);
}

template <class T>
int getIndex(const std::vector<T> &v, const T &elem) {
	auto ptr = lower_bound(v.begin(), v.end(), elem);
	if (ptr == v.end())
		throw std::runtime_error("World is coming to an end with elem");
	if (elem != *ptr)
		throw std::runtime_error("World is coming to an end with no matching element");
	return ptr - v.begin();
}

template <class T>
void deleteFromVector(std::vector<T> &v, const T &elem) {
	auto itr = std::find(v.begin(), v.end(), elem);
	if (itr != v.end())
		v.erase(itr);
}

template <class T>
void sortAndRemoveDuplicate(std::vector<T> &v) {
	std::sort(v.begin(), v.end());			// vector may have repeated elements like 1 1 2 2 3 3 3 4 4 5 5 6 7
	auto last = unique(v.begin(), v.end()); // vector now holds {1 2 3 4 5 6 7 x x x x x x}, where 'x' is indeterminate
	v.erase(last, v.end());
}

template <typename T, typename T2>
int addToMap(std::map<T, int> &mapping, std::vector<T2> &vec, const T &object) {
	if (mapping.find(object) == mapping.end()) {
		mapping.insert({object, vec.size()});
		vec.emplace_back();
	}
	return mapping[object];
}

template <typename T, int N>
void sortThis(std::array<T, N> &t) {
	std::sort(t.begin(), t.end());
}

template <typename T>
using myPair = std::array<T, 2>;

template <typename T>
using myTriple = std::array<T, 3>;

class ourEdge {
	int data[2];

  public:
	ourEdge(int a, int b) : data{a, b} {
		if (data[0] > data[1]) std::swap(data[0], data[1]);
	}
	int &operator[](int index) {
		//if (index < 0 || index > 1)
		//throw std::out_of_range("Index tried to access=" + index);
		return data[index];
	}
	const int &operator[](int index) const {
		//if (index < 0 || index > 1)
		//throw std::out_of_range("Index tried to access=" + index);
		return data[index];
	}
	int other(int u) {
		if (u == data[0]) return data[1];
		if (u == data[1]) return data[0];
		return -1;
	}
};

bool operator<(const ourEdge &lhs, const ourEdge &rhs) {
	if (lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]))
		return true;
	return false;
}

bool operator!=(const ourEdge &lhs, const ourEdge &rhs) {
	return lhs[0] != rhs[0] || lhs[1] != rhs[1];
}

bool operator==(const ourEdge &lhs, const ourEdge &rhs) {
	return lhs[0] == rhs[0] && lhs[1] == rhs[1];
}

std::ostream &operator<<(std::ostream &os, const ourEdge &e) {
	return os << "( " << e[0] << " " << e[1] << " )";
}

class ourFace {
	int data[3];

  public:
	ourFace(int a, int b, int c) : data{a, b, c} {
		if (data[0] > data[1]) std::swap(data[0], data[1]);
		if (data[1] > data[2]) std::swap(data[1], data[2]);
		if (data[0] > data[1]) std::swap(data[0], data[1]);
	}
	int &operator[](int index) {
		//if (index < 0 || index > 2)
		//throw std::out_of_range("Index tried to access=" + index);
		return data[index];
	}
	const int &operator[](int index) const {
		//if (index < 0 || index > 2)
		//throw std::out_of_range("Index tried to access=" + index);
		return data[index];
	}
};

bool operator<(const ourFace &lhs, const ourFace &rhs) {
	if (lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]) || (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] < rhs[2]))
		return true;
	return false;
}

bool operator!=(const ourFace &lhs, const ourFace &rhs) {
	return lhs[0] != rhs[0] || lhs[1] != rhs[1] || lhs[2] != rhs[2];
}

std::ostream &operator<<(std::ostream &os, const ourFace &f) {
	return os << "( " << f[0] << " " << f[1] << " " << f[2] << " )";
}

#endif // UTIL_H
