#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

typedef struct sol_node {
    int row;
    struct sol_node *prev;
} sol_node_t, *sol_t;

// Check if possition can be attacked (1) by currently positioned queens.
// The check is performed column-wise from right to left, checking the latest
// positioned queen, then the previous, etc.
static inline int check_attack(const int col, const int row, sol_t sol)
{
    int j;
    for (j = 0; j < col; j++) {
        const int tmp = abs(sol->row - row);
        if (tmp == 0 || tmp == j + 1)
            return 1;

        sol = sol->prev;
    }
    return  0;
}

int final_depth = 3;

static void solve(int n, const int col, sol_node_t& sol, int& result)
{
    if (col == n) {
        #pragma oss task reduction(+: result)
        result++;
    } else {
        for (int row = 0; row < n; row++) {
            if (!check_attack(col, row, &sol)) {
                sol_node_t new_sol;
                new_sol.prev = &sol;
                new_sol.row = row;
                #pragma oss task weakreduction(+: result) final(col >= final_depth)
                solve(n, col + 1, new_sol, result);
            }
        }
    }
}

int main(int argc, char **argv)
{
	int n = 14;

	if (argc > 3 || (argc >= 2 && strcmp(argv[1], "-h") == 0)) {
		fprintf(stderr, "Error: invalid number of arguments\n");
		fprintf(stderr, "Usage: nqueens [-h] n rows_in_parallel'\n");
		return 1;
	}

	if (argc >= 2)
		n = atoi(argv[1]);

	if (argc >= 3)
		final_depth = atoi(argv[2]);


	// Warmup iteration
	sol_node_t initial_node = {/* row */ -1, /* prev */ NULL};
	int count_main = 0;
	solve(n, 0, initial_node, count_main);
	#pragma oss taskwait

	initial_node = {/* row */ -1, /* prev */ NULL};
	count_main = 0;

	struct timeval tval_before, tval_after, tval_result;
	gettimeofday(&tval_before, NULL);

	solve(n, 0, initial_node, count_main);

	#pragma oss taskwait

	gettimeofday(&tval_after, NULL);

	timersub(&tval_after, &tval_before, &tval_result);

	double t = tval_result.tv_sec + 1e-6 * tval_result.tv_usec;

	printf("%14e %14d %14d %s\n", t, n, count_main, BENCH6_NAME);

	return 0;
}
