#include "checking.h"

bool is_latin_square(int *table, int n){
	int i, j, k, i2, j2, k2;
	for(k = 0, k2 = 0; k < n; k++, k2 += n)
		for(i = 0, i2 = 0; i < n; i++, i2 += n)
			for(j = i+1, j2 = n*(i+1); j < n; j++, j2 += n)
				if(table[k2 + i] == table[k2 + j]\
                || table[i2 + k] == table[j2 + k])
					return false;
	return true;
}

bool is_left_semiquasigroup(int *table, int n){
	for(int k = 0; k < n*n; k += n)
		for(int i = 0; i < n; i++)
			for(int j = i+1; j < n; j++)
				if(table[k + i] == table[k + j])
					return false;
	return true;
}

bool is_right_semiquasigroup(int *table, int n){
	for(int k = 0; k < n; k++)
		for(int i = 0; i < n*n; i += n)
			for(int j = i+n; j < n*n; j += n)
				if(table[i + k] == table[j + k])
					return false;
	return true;
}

bool is_associative(int *table, int n){
	int i, j, k, j2;
	for(i = 0; i < n*n; i += n)
		for(j = 0, j2 = 0; j < n; j++, j2 += n)
			for(k = 0; k < n; k++)
				if(table[table[i + j]*n + k]\
                   != table[i + table[j2 + k]])
					return false;
	return true;
}

bool is_identity(int *table, int n){
	bool is_id;
	int i, i2, j, j2;
	for(i = 0, i2 = 0; i < n; i++, i2 += n){
		is_id = true;
		for(j = 0, j2 = 0; j < n; j++, j2 += n){
			if(table[i2+j] != j || table[j2+i] != j){
				is_id = false;
				break;
			}
		}
		if(is_id) return true;
	}
	return false;
}

bool is_commutative(int *table, int n){
	int i, j, i2, j2;
	for(i = 0, i2 = 0; i < n; i++, i2 += n)
		for(j = 0, j2 = 0; j < n; j++, j2 += n)
			if(table[i2 + j] != table[j2 + i])
                return false;
	return true;
}

bool is_left_identity(int *table, int n){
	int count, i, j;
	for(i = 0; i < n*n; i += n){
		count = 0;
		for(j = 0; j < n; j++)
			if(table[i+j] == j) count++;
		if(count == n) return true;
	}
	return false;
}

bool is_right_identity(int *table, int n){
	int count, i, j, j2;
	for(i = 0; i < n; i++){
		count = 0;
		for(j = 0, j2 = 0; j < n; j++, j2 += n)
			if(table[j2 + i] == j) count++;
		if(count == n) return true;
	}
	return false;
}
