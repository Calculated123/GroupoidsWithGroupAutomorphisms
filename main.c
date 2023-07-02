#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>
#include "funcs.h"
#include "checking.h"

// ��� ���������
#include <windows.h>

int main(){
    // ��� ���������
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
	int n;
	printf("������� ������ ����������� w: ");
	scanf("%i", &n);
	int perm[n];
	printf("������� ����������� w (����� �� 1 �� %i):\n", n);
	for(int i = 0; i < n; i++){
		scanf("%i", &perm[i]);
		perm[i]--;
	}
	int num_of_cycles = number_of_cycles(perm, n);
	int *len_of_cycles = length_of_cycles(perm, n);
	int **cycle_repr = get_cycle_representation(perm, n);
	printf("���������� w � ������������ ������� ���������������� ������:\n");
	for(int i = 0; i < num_of_cycles; i++){
		printf("(");
		for(int j = 0; j < len_of_cycles[i]-1; j++){
			printf("%i ", cycle_repr[i][j]+1);
		}
		printf("%i)", cycle_repr[i][len_of_cycles[i]-1]+1);
	}
	printf("\n������ H, ���������� �������� ������������ w,\n�������� %i ���������.\n\n", power_of_perm(perm, n));

	int num_of_orbits = number_of_orbits(len_of_cycles, num_of_cycles);
	printf("��� �������� ���������� ������\n�� ��������� N^2_n, ��� ���������� �� %i �����.\n\n", num_of_orbits);
	int *orbit_matrix = calloc(n*n, sizeof(int));
    int count = 1;
    for(int i = 0; i < n*n; orbit_matrix[i] = 0, i++);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(orbit_matrix[i*n+j] == 0){
                orbit_matrix[i*n + j] = count;
                int t_i = perm[i], t_j = perm[j];
                while(t_i != i || t_j != j){
                    orbit_matrix[t_i*n + t_j] = count;
                    t_i = perm[t_i];
                    t_j = perm[t_j];
                }
                count++;
            }
        }
    }
    printf("������� ����� OH:\n");
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%2i ", orbit_matrix[i*n+j]);
        }
        printf("\n");
    }
	int count_for_orbits = 0;
	int *beg_of_orbits = calloc(2 * num_of_orbits, sizeof(int));
	int *len_of_orbits = calloc(num_of_orbits, sizeof(int));
	for(int i = 0; i < num_of_cycles; i++){
		for(int j = 0; j < num_of_cycles; j++){
			for(int k = 0; k < gcd(len_of_cycles[i], len_of_cycles[j]); k++){
				beg_of_orbits[2*count_for_orbits] = cycle_repr[i][0];
				beg_of_orbits[2*count_for_orbits + 1] = cycle_repr[j][k];
				len_of_orbits[count_for_orbits] = lcm(len_of_cycles[i], len_of_cycles[j]);
				count_for_orbits++;
			}
		}
	}
    printf("������� BH:\n");
    for(int i = 0; i < num_of_orbits; i++){
        printf("%2i) (%i, %i)\n", i+1, beg_of_orbits[2*i]+1, beg_of_orbits[2*i+1]+1, len_of_orbits[i]);
    }
	int number_of_allowable_symbols = 0;
	int *len_of_allowable_symbols = calloc(num_of_orbits, sizeof(int));
	for(int i = 0; i < num_of_orbits; i++){
		int t_n_o_a_f = 0;
		for(int j = 0; j < num_of_cycles; j++){
			if(len_of_orbits[i] % len_of_cycles[j] == 0)
                t_n_o_a_f += len_of_cycles[j];
		}
		len_of_allowable_symbols[i] = t_n_o_a_f;
		number_of_allowable_symbols += t_n_o_a_f;
	}
	int temp_count_for_pointers = 0;
	int *under_allowable_symbols = calloc(number_of_allowable_symbols, sizeof(int));
	int **allowable_symbols = calloc(num_of_orbits, sizeof(int*));
	for(int i = 0; i < num_of_orbits; i++){
		allowable_symbols[i] = &under_allowable_symbols[temp_count_for_pointers];
		temp_count_for_pointers += len_of_allowable_symbols[i];
	}
	for(int i = 0; i < num_of_orbits; i++){
		int temp_count = 0;
		for(int j = 0; j < num_of_cycles; j++){
			if(len_of_orbits[i] % len_of_cycles[j] == 0){
				for(int k = 0; k < len_of_cycles[j]; k++, temp_count++){
					allowable_symbols[i][temp_count] = cycle_repr[j][k];
				}
			}
		}
	}
    printf("\n��������� ���������� �������� LH\n��� Orb_1,...,Orb_%i:\n", num_of_orbits);
    for(int i = 0; i < num_of_orbits; i++){
        printf("LH(%i): ", i+1);
        for(int j = 0; j < len_of_allowable_symbols[i]; j++)
            printf("%i ", allowable_symbols[i][j]+1);
        printf("\n");
    }
	unsigned long long num_of_variants = 1;
	for(int i = 0; i < num_of_orbits; i++)
		num_of_variants *= len_of_allowable_symbols[i];
	// �������� ��������������
	omp_set_num_threads(1);
		// �������� ��������� ����� ����������
		// [0] - ������� ������
		// [1] - ������
		// [2] - ������������� ����
		// [3] - ����
		// [4] - ������������� �����������
		// [5] - �����������
		// [6] - ������������� ������
		// [7] - ������
		// [8] - ������������� ����������
		// [9] - ����������
	int number_of_structures_par[10][omp_get_max_threads()];
	for(int i = 0; i < 10; i++)
		for(int j = 0; j < omp_get_max_threads(); j++)
			number_of_structures_par[i][j] = 0;
		// �������� ������� ���������
		// �������� ������� ��������� ��� ������� ������
		// � ������ �������� ��� ������� �� �������!!!
	int num_of_threads = omp_get_max_threads();
	int counters[num_of_threads][num_of_orbits];
	for(int i = 0; i < num_of_threads; i++){
		for(int j = 0; j < num_of_orbits; j++){
			counters[i][j] = 0;
		}
		next_iteration(counters[i], num_of_orbits, len_of_allowable_symbols, i);
	}
	int tables[num_of_threads][n*n];
	double start, finish;
		// ������ �������� �������
		// ��� ��� ���������� ��������� �� ������ ������� �� ���������� �������, ������� ����� (�������������)
		// ��������� �������� ������ (������� �� ������� ���-�� ��������� �� ���-�� �����),
		// ����� ����� �������� ������������ ���� (� ������� ��� ���������� �������� ���������� ���������)
	start = omp_get_wtime();
	int excess_variants = num_of_variants % num_of_threads;
	for(int k = 0; k < excess_variants; k++){
		int thread_num = k;
		int *cur_table = tables[thread_num];
		int *cur_counters = counters[thread_num];
		for(int i = 0, i_2 = 0; i < num_of_orbits; i++, i_2 += 2){
			int fix_i = beg_of_orbits[i_2], fix_j = beg_of_orbits[i_2+1], fix_k = allowable_symbols[i][cur_counters[i]];
			cur_table[fix_i*n + fix_j] = fix_k;
			for(int j = 0; j < len_of_orbits[i]-1; j++){
				fix_i = perm[fix_i]; fix_j = perm[fix_j]; fix_k = perm[fix_k];
				cur_table[fix_i*n + fix_j] = fix_k;
			}
		}
		bool is_lc = is_latin_square(cur_table, n),is_assoc = is_associative(cur_table, n),\
			 is_id = is_identity(cur_table, n), is_comm = is_commutative(cur_table, n);
		if(is_lc && is_assoc && is_comm) number_of_structures_par[0][thread_num]++;
		if(is_lc && is_assoc) number_of_structures_par[1][thread_num]++;
		if(is_lc && is_id && is_comm) number_of_structures_par[2][thread_num]++;
		if(is_lc && is_id) number_of_structures_par[3][thread_num]++;
		if(is_lc && is_comm) number_of_structures_par[4][thread_num]++;
		if(is_lc) number_of_structures_par[5][thread_num]++;
        printf("!\n");
		if(is_assoc && is_id && is_comm) number_of_structures_par[6][thread_num]++;
		if(is_assoc && is_id) number_of_structures_par[7][thread_num]++;
		if(is_assoc && is_comm) number_of_structures_par[8][thread_num]++;
		if(is_assoc) number_of_structures_par[9][thread_num]++;
		next_iteration(cur_counters, num_of_orbits, len_of_allowable_symbols, num_of_threads);
	}
	int kolvo = (num_of_variants - excess_variants)/num_of_threads;
#pragma omp parallel for schedule(static, kolvo)
	for(unsigned long long k = 0; k < num_of_variants - excess_variants; k++){
		int thread_num = omp_get_thread_num();
		int *cur_table = tables[thread_num];
		int *cur_counters = counters[thread_num];
		int fix_i, fix_j, fix_k;
		for(int i = 0, i_2 = 0; i < num_of_orbits; i++, i_2 += 2){
			fix_i = beg_of_orbits[i_2], fix_j = beg_of_orbits[i_2+1], fix_k = allowable_symbols[i][cur_counters[i]];
			for(int j = 0; j < len_of_orbits[i]; j++, fix_i = perm[fix_i], fix_j = perm[fix_j], fix_k = perm[fix_k]){
				cur_table[fix_i*n+fix_j] = fix_k;
			}
		}
		bool is_lc = is_latin_square(cur_table, n),is_assoc = is_associative(cur_table, n),\
			 is_id = is_identity(cur_table, n), is_comm = is_commutative(cur_table, n);
		if(is_lc && is_assoc && is_comm) number_of_structures_par[0][thread_num]++;
		if(is_lc && is_assoc) number_of_structures_par[1][thread_num]++;
		if(is_lc && is_id && is_comm) number_of_structures_par[2][thread_num]++;
		if(is_lc && is_id) number_of_structures_par[3][thread_num]++;
		if(is_lc && is_comm) number_of_structures_par[4][thread_num]++;
		if(is_lc) number_of_structures_par[5][thread_num]++;
		if(is_assoc && is_id){
            tech_print(cur_table, n);
            printf("\n");
		}
		if(is_assoc && is_id && is_comm) number_of_structures_par[6][thread_num]++;
		if(is_assoc && is_id) number_of_structures_par[7][thread_num]++;
		if(is_assoc && is_comm) number_of_structures_par[8][thread_num]++;
		if(is_assoc) number_of_structures_par[9][thread_num]++;
		int temp_k = num_of_orbits-1;
		cur_counters[temp_k] += num_of_threads;
		while(temp_k > 0){
			if(cur_counters[temp_k] >= len_of_allowable_symbols[temp_k]){
				cur_counters[temp_k-1] += cur_counters[temp_k] / len_of_allowable_symbols[temp_k];
				cur_counters[temp_k] %= len_of_allowable_symbols[temp_k];
				temp_k--;
			}
			else break;
		}
	}
	finish = omp_get_wtime();
	printf("\n�������� ���������� ����������� w:\n");
	int *number_of_structures = calloc(10, sizeof(int));
	for(int i = 0; i < 10; i++){
		number_of_structures[i] = 0;
		for(int j = 0; j < num_of_threads; j++)
			number_of_structures[i] += number_of_structures_par[i][j];
	}
	printf("�������� �����: %i\n", number_of_structures[0]);
	printf("�����: %i\n", number_of_structures[1]);
	printf("������������� ���: %i\n", number_of_structures[2]);
	printf("���: %i\n", number_of_structures[3]);
	printf("������������� ����������: %i\n", number_of_structures[4]);
	printf("����������: %i\n", number_of_structures[5]);
	printf("������������� ��������: %i\n", number_of_structures[6]);
	printf("��������: %i\n", number_of_structures[7]);
	printf("������������� ���������: %i\n", number_of_structures[8]);
	printf("���������: %i\n", number_of_structures[9]);
	printf("����������: %lld\n", num_of_variants);
	printf("����� ����������: %g\n", finish - start);
	return 0;
}
