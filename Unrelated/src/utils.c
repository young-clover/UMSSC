#include "def.h"

void usage(const char* prog_name) {
    fprintf(stderr, "\nUsage: %s <input_file> <P> <Method> [mode]\n", prog_name);

    fprintf(stderr, "\nArguments:\n");
    fprintf(stderr, "  <input_file> : Path to the file containing instance names.\n");
    fprintf(stderr, "  <P>          : The P parameter (double value).\n");
    fprintf(stderr, "  <Method>     : 'SOCP' or 'PBD'.\n");
    fprintf(stderr, "  [mode]       : (Optional) 'single' or 'multi'. Defaults to 'single' if omitted.\n"); 

    fprintf(stderr, "\nExamples:\n");
    fprintf(stderr, "  %s list.txt 2.0 SOCP\n", prog_name);
    fprintf(stderr, "  %s list.txt 2.0 PBD          (Runs as single)\n", prog_name); 
    fprintf(stderr, "  %s list.txt 2.0 PBD multi    (Runs as multi)\n\n", prog_name);

    exit(1);
}

FILE* open_file(const char* name, const char* mode) {
    FILE* file = fopen(name, mode);
    if (!file) {
        fprintf(stderr, "Error: Failed to open file %s\n", name);
        exit(8);
    }
    return file;
}

void read_instance(const char* name, const char* type) {
    char path[300];
    sprintf(path, "%s%s", instances_path, name);

    FILE* inst = open_file(path, "r");

    fscanf(inst, "%d", &M);
    fscanf(inst, "%d", &N);

	// Allocate memory for instance data after reading M and N  
    initialize_memory();


    // fixed_charge
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
			fscanf(inst, "%lf", &fixed_charge[i][j]);//%lf for reading double values
    // weight
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            fscanf(inst, "%lf", &weight[i][j]);
    // capacity
    for (int i = 0; i < M; i++)
        fscanf(inst, "%lf", &capacity[i]);
    // lower_bound
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            fscanf(inst, "%lf", &lower_bound[i][j]);
    // upper_bound
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            fscanf(inst, "%lf", &upper_bound[i][j]);


    fclose(inst);
}

void initialize_memory(void) {
	// reset parameters
    LowerBound = -GRB_INFINITY;
	LP_UpperBound = GRB_INFINITY;
	UpperBound = GRB_INFINITY;
	num_of_unfixed = M * N;
	pre_fix_sum = 0.0;

	// Allocate memory for master problem data structures
    Stabilizing_Point = create_double_vector(M * N);
    ReducedCost = create_double_vector(M * N);
	EqualityPi = create_double_vector(N);
    

    num_of_unfixed = M * N;// initially all variables are unfixed
    unfixed_index = create_int_vector(M * N);
    fixed_point = create_double_vector(M * N);
    for (int i = 0; i < M * N; i++) {
        fixed_point[i] = -1.0;// -1 means unfixed        
        unfixed_index[i] = i;
    }


    benders_cut = (BENDERSCUT*)malloc(sizeof(BENDERSCUT));
    if (pbd_mode == MODE_SINGLE) {
		// single-cut
        Y = create_double_vector(M * N + 1);
        Separating_Point = create_double_vector(M * N + 1);

        theta_try_solution = create_double_vector(1);

        best_try_solution = create_double_vector(M * N + 1);
        rounding_solution = create_double_vector(M * N + 1);
        greedy_solution = create_double_vector(M * N + 1);

        benders_cut->beg = create_int_vector(1);
        benders_cut->ind = create_int_vector(M * N + 1);
        benders_cut->val = create_double_vector(M * N + 1);
        benders_cut->sense = create_char_vector(1);
        benders_cut->rhs = create_double_vector(1);

        // Initialize Benders cut structure
        benders_cut->sense[0] = GRB_LESS_EQUAL;
        benders_cut->beg[0] = 0;
        for (int i = 0; i < M * N; i++) {
            benders_cut->ind[i] = i;
        }
        benders_cut->ind[M * N] = M * N;
        benders_cut->val[M * N] = -1;
    }

    if (pbd_mode == MODE_MULTI) {
        // multi-cut
        Y = create_double_vector(M * N + M);
        Separating_Point = create_double_vector(M * N + M);

        theta_try_solution = create_double_vector(M);

        best_try_solution = create_double_vector(M * N + M);
        rounding_solution = create_double_vector(M * N + M);
        greedy_solution = create_double_vector(M * N + M);

        benders_cut->beg = create_int_vector(M);
        benders_cut->ind = create_int_vector(M * N + M);
        benders_cut->val = create_double_vector(M * N + M);
        benders_cut->sense = create_char_vector(M);
        benders_cut->rhs = create_double_vector(M);

        // Initialize Benders cut structure
        for (int i = 0; i < M; i++) {        
            benders_cut->beg[i] = i * (N + 1);
            for (int j = 0; j < N; j++) {
                benders_cut->ind[j + i * (N + 1)] = i * N + j;
            }
            benders_cut->ind[N + i * (N + 1)] = M * N + i;
            benders_cut->val[N + i * (N + 1)] = -1;
            benders_cut->sense[i] = GRB_LESS_EQUAL;
		}
        num_of_coeff_per_cut = create_int_vector(M);
        for (int i = 0; i < M; i++) {
            num_of_coeff_per_cut[i] = N + 1;
        }
    }
    
	// Allocate memory for instance data
    fixed_charge = create_double_matrix(M, N);
    weight = create_double_matrix(M, N);
    capacity = create_double_vector(M);
    lower_bound = create_double_matrix(M, N);
    upper_bound = create_double_matrix(M, N);

	// Allocate memory for subproblem data structures
    S = create_int_vector(N);
	L = create_double_vector(N);
	U = create_double_vector(N);
	W = create_double_vector(N);
	alpha = create_double_vector(N);
	beta = create_double_vector(N);
    x = create_double_vector(M * N); // solution vector for continuous variables
    subgradient = create_double_vector(M * N);

	// Allocate memory for greedy heuristic
    remainingCap = (double*)malloc(M * sizeof(double));
    jobs = (JobRcDiff*)malloc(N * sizeof(JobRcDiff));
    
}

void free_memory(void) {

    free_double_matrix(fixed_charge, M);
    free_double_matrix(weight, M);
    free(capacity);
    free_double_matrix(lower_bound, M);
    free_double_matrix(upper_bound, M);


    free(Stabilizing_Point);
    free(ReducedCost);
	free(EqualityPi);

    free(theta_try_solution);

    free(Separating_Point);
    free(Y);

	free(best_try_solution);
	free(rounding_solution);
	free(greedy_solution);

	free(fixed_point);
	free(unfixed_index);

    free(benders_cut->rhs);
    free(benders_cut->sense);
    free(benders_cut->beg);
    free(benders_cut->ind);
    free(benders_cut->val);
    free(benders_cut);     

    if (pbd_mode == MODE_MULTI)
        free(num_of_coeff_per_cut);

    free(S);
    free(x);
    free(L);
    free(U);
    free(W);
    free(alpha);
    free(beta);
	free(subgradient);

    free(remainingCap);
    free(jobs);
}

int varindex(int m, int n) {
    return m * N + n;
}

int* create_int_vector(int n) {
    int* ptr = (int*)calloc(n, sizeof(int));
    if (!ptr) { printf("Memory error\n"); exit(8); }
    return ptr;
}

double* create_double_vector(int n) {
    double* ptr = (double*)calloc(n, sizeof(double));
    if (!ptr) { printf("Memory error\n"); exit(8); }
    return ptr;
}

double** create_double_matrix(int m, int n) {
    double** ptr = (double**)calloc(m, sizeof(double*));
    if (!ptr) { printf("Memory error\n"); exit(8); }
    for (int i = 0; i < m; i++)
        ptr[i] = create_double_vector(n);
    return ptr;
}

void free_double_matrix(double** ptr, int m) {
    for (int i = 0; i < m; i++) free(ptr[i]);
    free(ptr);
}

char* create_char_vector(int n) {
    char* ptr = (char*)calloc(n, sizeof(char));
    if (!ptr) { printf("Memory error\n"); exit(8); }
    return ptr;
}

char** create_stringarray(int m, int n) {
    char** ptr = (char**)calloc(m, sizeof(char*));
    if (!ptr) { printf("Memory error\n"); exit(8); }
    for (int i = 0; i < m; i++) ptr[i] = create_char_vector(n);
    return ptr;
}

void free_stringarray(char** ptr, int m) {
    for (int i = 0; i < m; i++) free(ptr[i]);
    free(ptr);
}
