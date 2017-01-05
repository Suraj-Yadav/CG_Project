#ifndef UTIL_H
#define UTIL_H

#include <array>
#include <vector>

#define val(x) cout << #x "=" << x << "\n"

#define print(...) debug_print(__VA_ARGS__)
#define println(...) debug_print(__VA_ARGS__ ,"\n")

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

template <typename T, size_t N>
void sortThis(std::array<T, N> &t) {
	sort(t.begin(), t.end());
}

template <typename T>  
using myPair = std::array<T, 2>;

template <typename T>
using myTriple = std::array<T, 3>;

#endif // UTIL_H

