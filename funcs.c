#include "funcs.h"
#include "math_funcs.h"

int *mult(int *perm1, int *perm2, int n){
	int *res = calloc(n, sizeof(int));
	for(int i = 0; i < n; i++)
		res[i] = perm1[perm2[i]];
	return res;
}

bool is_identity_perm(int *perm, int n){
	for(int i = 0; i < n; i++)
		if(perm[i] != i) return false;
	return true;
}

int power_of_perm(int *perm, int n){
	int p = 1;
	int *temp = calloc(n, sizeof(int));
	for(int i = 0; i < n; temp[i] = perm[i], i++);
	while(!is_identity_perm(temp, n)){
		temp = mult(temp, perm, n);
		p++;
	}
	free(temp);
	return p;
}

int *get_cyclic_group(int *perm, int n, int size){
	int *group = calloc(n*size, sizeof(int));
	int *temp = calloc(n, sizeof(int));
	for(int i = 0; i < n; temp[i] = i, i++);
	for(int i = 0; i < n*size; i += n){
		for(int j = 0; j < n; j++)
			group[i+j] = temp[j];
		temp = mult(temp, perm, n);
	}
	free(temp);
	return group;
}

int number_of_cycles(int *perm, int n){
	int number = 0;
	int *temp = calloc(n, sizeof(int));
	for(int i = 0; i < n; temp[i] = perm[i], i++);
	for(int i = 0; i < n; i++){
		if(temp[i] == -1) continue;
		int k = i;
		while(temp[k] != -1){
			temp[k] = -1;
			k = perm[k];
		}
		number++;
	}
	free(temp);
	return number;
}

// Размер массива, содержащий длины циклов имеет длину: num_of_cycles (из другой функции)
int *length_of_cycles(int *perm, int n){
	int num_of_cycles = number_of_cycles(perm, n), cycle_count = 0;
	int *temp = calloc(n, sizeof(int));
	for(int i = 0; i < n; temp[i] = -1, i++);
	int *res = calloc(num_of_cycles, sizeof(int));
	for(int k = 0; k < n; k++){
		bool is_in = false;
		for(int i = 0; i < n; i++)
			if(temp[i] == k) is_in = true;
		if(is_in) continue;
		int count = 1;
		int fix = k;
		temp[fix] = perm[fix];
		fix = perm[fix];
		while(fix != k){
			temp[fix] = perm[fix];
			fix = perm[fix];
			count++;
		}
		res[cycle_count++] = count;
	}
	free(temp);
	return res;
}

int **get_cycle_representation(int *perm, int n){
	int count = 0, cycle_count = 0, num_of_cycles = number_of_cycles(perm, n);
	int *under_res = calloc(n, sizeof(int));
	int **res = calloc(num_of_cycles, sizeof(int*));
	for(int i = 0; i < n; under_res[i] = -1, i++);
	for(int k = 0; k < n; k++){
		bool is_in = false;
		for(int i = 0; i < count; i++)
			if(k == under_res[i]) is_in = true;
		if(is_in) continue;
		res[cycle_count++] = &under_res[count];
		under_res[count++] = k;
		int fix = perm[k];
		while(fix != k){
			under_res[count++] = fix;
			fix = perm[fix];
		}
	}
	return res;
}

int number_of_orbits(int *len_of_cycles, int num_of_cycles){
	int s = 0;
	for(int i = 0; i < num_of_cycles; i++)
		for(int j = 0; j < num_of_cycles; j++)
			s += gcd(len_of_cycles[i], len_of_cycles[j]);
	return s;
}

int *get_orbit_matrix(int *perm, int n){
    int *res = calloc(n*n, sizeof(int));
    int count = 1;
    for(int i = 0; i < n*n; res[i] = 0, i++);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(res[i*n+j] == 0){
                res[i*n + j] = count;
                int t_i = perm[i], t_j = perm[j];
                while(t_i != i || t_j != j){
                    res[t_i*n + t_j] = count;
                    t_i = perm[t_i];
                    t_j = perm[t_j];
                }
                count++;
            }
        }
    }
    return res;
}

void next_iteration(int *counter, int number, int *sizes, int step){
	int k = number-1;
	counter[k] += step;
	while(k > 0){
		if(counter[k] >= sizes[k]){
			counter[k-1] += counter[k] / sizes[k];
			counter[k] %= sizes[k];
			k--;
		}
		else break;
	}
}

void print_table(int *table, int n){
	for(int i = 0; i < n*n; i += n){
		for(int j = 0; j < n; j++){
			printf("%i, ", table[i+j]);
		}
		printf("\n");
	}
	printf("\n");
}

void tech_print(int *table, int n){
    printf("[\n");
    for(int i = 0; i < n; i++){
        printf("[");
        for(int j = 0; j < n-1; j++){
            printf("%i, ", table[i*n + j]);
        }
        printf("%i],\n", table[i*n+n-1]);
    }
    printf("],");
}

void next_perm(int *perm, int n){
    int temp, fix_j = -1;
    for(int j = n-2; j >= 0; j--){
        if(perm[j] < perm[j+1]){
            fix_j = j;
            break;
        }
    }
    if(fix_j == -1){
        for(int i = 0; i < n/2; i++){
            temp = perm[i];
            perm[i] = perm[n-i-1];
            perm[n-i-1] = temp;
        }
        return;
    }
    for(int l = n-1; l > fix_j; l--){
        if(perm[l] > perm[fix_j]){
            temp = perm[l];
            perm[l] = perm[fix_j];
            perm[fix_j] = temp;
            break;
        }
    }
    for(int i = 0; i < (n - fix_j + 1)/2 - 1; i++){
        temp = perm[fix_j+i+1];
        perm[fix_j+i+1] = perm[n-i-1];
        perm[n-i-1] = temp;
    }
}
