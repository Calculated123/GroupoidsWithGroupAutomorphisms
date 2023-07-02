#include "math_funcs.h"

int gcd(int a, int b){
	int result = ((a < b) ? a : b);
	while (result > 0) {
		if (a % result == 0 && b % result == 0) {
			break;
        }
		result--;
    }
    return result; // return gcd of a nd b
}

int lcm(int a, int b){
	return a * b / gcd(a, b);
}
