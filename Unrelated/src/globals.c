#include "def.h"

char input_text[256];
char output_text[256];
char instances_path[256];
char instance_name[256];
FILE* input;
FILE* output;

double Gap = 1.0;

int pbd_mode = -1; // global variable: 0=SINGLE, 1=MULTI

int CONTINUE = 0; // whether to continue to the next module

double LowerBound;
double LP_UpperBound;
double UpperBound;

double RootRuntime;
double HeurRuntime;
double FixingRuntime;

double P;
double P_array[4] = { 1.5, 2.0, 2.5, 3.0 };
int M;
int N;

double* Y;
double* ReducedCost;
double* EqualityPi;// dual variables for equality constraints in master

double* Stabilizing_Point;
double* Separating_Point;
BENDERSCUT* benders_cut;

double* best_try_solution;
double* theta_try_solution;

double* rounding_solution;
double* greedy_solution;

int num_of_unfixed;
double pre_fix_sum;
double* fixed_point;
int* unfixed_index;
int* num_of_coeff_per_cut;

int* S;
double* x;
double* L;
double* U;
double* W;//temp_array for calculations
double* alpha;
double* beta;
double* subgradient;

// Data for instances
double** fixed_charge;
double** weight;
double* capacity;
double** lower_bound;
double** upper_bound;

// Variables for greedy heuristic
double* remainingCap;
JobRcDiff* jobs;




