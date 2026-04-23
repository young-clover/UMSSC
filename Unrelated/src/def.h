#ifndef DEF_H
#define DEF_H

#define char_array_size 30

#define MODE_SINGLE 0 // single cut mode
#define MODE_MULTI  1 // multi cut mode

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h> // Sleep function
#include "gurobi_c.h"
#include <time.h>
#include <math.h> // pow function
#include <stdbool.h> // bool type

//Benders cut data structure
typedef struct {
    int* beg;
    int* ind;
    double* val;
    char* sense;
    double* rhs;
} BENDERSCUT;

// structure for greedy heuristic
typedef struct {
    int job;
    double diff; //difference in reduced cost
} JobRcDiff;


// Global variables
extern char input_text[256];
extern char output_text[256];
extern char instances_path[256];
extern char instance_name[256];
extern FILE* input;
extern FILE* output;

extern double Gap;

extern int pbd_mode; // global variable: 0=SINGLE, 1=MULTI

extern int CONTINUE; // whether to continue to the next module

extern double LowerBound;
extern double LP_UpperBound;
extern double UpperBound;

extern double RootRuntime;
extern double HeurRuntime;
extern double FixingRuntime;


extern double P;
extern double P_array[4];
extern int M;
extern int N;

extern double* Y; // Solution vector from the master problem
extern double* Stabilizing_Point;
extern double* Separating_Point;
extern BENDERSCUT* benders_cut;// Benders cut data structure

extern double* ReducedCost;
extern double* EqualityPi;// dual variables for equality constraints in master

extern double* best_try_solution;
extern double* theta_try_solution;

extern double* rounding_solution;
extern double* greedy_solution;


extern int num_of_unfixed;
extern double pre_fix_sum;
extern double* fixed_point;
extern int* unfixed_index;
extern int* num_of_coeff_per_cut;// only for multi-cut


extern int* S;
extern double* x;
extern double* L;
extern double* U;
extern double* W; // temp_array for calculations
extern double* alpha;
extern double* beta;
extern double* subgradient;

extern double* remainingCap;
extern JobRcDiff* jobs;


// Data for instances
extern double** fixed_charge;
extern double** weight;
extern double* capacity;
extern double** lower_bound;
extern double** upper_bound;

// Function prototypes
void usage(const char* progname);
FILE* open_file(const char* name, const char* mode);
void read_instance(const char* name, const char* type);
void initialize_memory(void);
void free_memory(void);
int varindex(int m, int n);

int* create_int_vector(int n);
double* create_double_vector(int n);
char* create_char_vector(int n);
char** create_stringarray(int m, int n);
void free_stringarray(char** ptr, int m);
double** create_double_matrix(int m, int n);
void free_double_matrix(double** ptr, int m);

void PBD(void);
void MasterProblem(GRBenv* env, GRBmodel* model);
void CutSeparator(void);
int cmpJobRcDescendingDiff(const void* a, const void* b);
void RunHeuristics(void);
void checkSolutionFeasibility(void);

void RootNodeSolving(GRBenv* env, GRBmodel* model);
void LocalSearch(void);
void ReducedCostFixing(GRBenv* env, GRBmodel* model);
void MIPSolving(GRBenv* env, GRBmodel* model);
int __stdcall SingleBendersCutCallback(GRBmodel* model, void* cbdata, int where, void* usrdata);
int __stdcall MultiBendersCutCallback(GRBmodel* model, void* cbdata, int where, void* usrdata);

double MKP(GRBenv* env, int m1, int m2, int* lns_sub_tasks, int sub_task_count);

void SOCP(void);
#endif
