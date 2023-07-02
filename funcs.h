#include <malloc.h>
#include <stdbool.h>
#include "math_funcs.h"

int *mult(int *perm1, int *perm2, int n);
	
bool is_identity_perm(int *perm, int n);

int power_of_perm(int *perm, int n);

int *get_cyclic_group(int *perm, int n, int size);

int number_of_cycles(int *perm, int n);

int *length_of_cycles(int *perm, int n);

int **get_cycle_representation(int *perm, int n);

int number_of_orbits(int *len_of_cycles, int num_of_cycles);

void next_iteration(int *counter, int number, int *sizes, int step);

void print_table(int *table, int n);
