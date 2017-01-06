#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <array>
#include <vector>
#include <map>

#define val(x) cout << #x "=" << x << "\n"

#define print(...) debug_print(__VA_ARGS__)
#define println(...) {debug_print(__VA_ARGS__);debug_print('\n');}

void debug_print() {
}
template <typename T, typename... Args>
void debug_print(T t, Args... args) {
	std::cout << t << " ";
	debug_print(args...);
}

template <class T>
size_t getIndex(const std::vector<T> &v, const T &elem) {
	return lower_bound(v.begin(), v.end(), elem) - v.begin();
}

template <class T>
void sortAndRemoveDuplicate(std::vector<T> &v) {
	sort(v.begin(), v.end());				  // vector may have repeated elements like 1 1 2 2 3 3 3 4 4 5 5 6 7
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
	sort(t.begin(), t.end());
}

template <typename T>  
using myPair = std::array<T, 2>;

template <typename T>
using myTriple = std::array<T, 3>;

#endif // UTIL_H

