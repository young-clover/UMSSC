#include "def.h"

void SOCP(void) {
    GRBenv* env = NULL;
    GRBmodel* model = NULL;
    int       error = 0;

    // Create environment
    error = GRBemptyenv(&env);
    if (error) goto QUIT;
    error = GRBsetintparam(env, "OutputFlag", 0);
    if (error) goto QUIT;   
    error = GRBstartenv(env);
    if (error) goto QUIT;

    //Create an empty model and  Add variables
    int numvars = N * M + 2 * N * M;
    double* obj = create_double_vector(numvars);
    double* lb = create_double_vector(numvars);
    double* ub = create_double_vector(numvars);
    char* vtype = create_char_vector(numvars);
    char** varnames = create_stringarray(numvars, char_array_size);
    if (obj == NULL || lb == NULL || ub == NULL || vtype == NULL || varnames == NULL) {
        printf("Memory allocation failed for variables.\n");
        exit(1);
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            obj[varindex(i, j)] = fixed_charge[i][j];
            lb[varindex(i, j)] = 0;
            ub[varindex(i, j)] = 1;
            vtype[varindex(i, j)] = GRB_BINARY;// GRB_BINARY;
            sprintf(varnames[varindex(i, j)], "y_(%d,%d)", i, j);

            obj[varindex(i, j) + M * N] = 0.0;
            lb[varindex(i, j) + M * N] = 0;
            ub[varindex(i, j) + M * N] = GRB_INFINITY;
            vtype[varindex(i, j) + M * N] = GRB_CONTINUOUS;
            sprintf(varnames[varindex(i, j) + M * N], "x_(%d,%d)", i, j);

            obj[varindex(i, j) + N * M * 2] = pow(weight[i][j], P);   // weight^P: all machines use same P
            if (P == 0)// hibrid case with different P values
                obj[varindex(i, j) + N * M * 2] = pow(weight[i][j], P_array[i % 4]); // different P for each machine
            lb[varindex(i, j) + N * M * 2] = 0;
            ub[varindex(i, j) + N * M * 2] = GRB_INFINITY;
            vtype[varindex(i, j) + N * M * 2] = GRB_CONTINUOUS;
            sprintf(varnames[varindex(i, j) + N * M * 2], "z_(%d,%d)", i, j);
        }
    }

    // new model 
    error = GRBnewmodel(env, &model, "misocp", numvars, obj, lb, ub, vtype, varnames);
    if (error) goto QUIT;

    // Free allocated memory for variables
    free(obj);
    free(lb);
    free(ub);
    free(vtype);
    free_stringarray(varnames, char_array_size);

    /* Change sense to minimization */
    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);
    if (error) goto QUIT;

    // Add linear constraints to the model
    int numconstrs = 2 * N * M + M + N;
    int numnz = 6 * N * M;
    int* cbeg = create_int_vector(numconstrs);
    int* cind = create_int_vector(numnz);
    double* cval = create_double_vector(numnz);
    char* sense = create_char_vector(numconstrs);
    double* rhs = create_double_vector(numconstrs);
    char** constrnames = create_stringarray(numconstrs, char_array_size);
    if (cbeg == NULL || cind == NULL || cval == NULL || sense == NULL || rhs == NULL || constrnames == NULL) {
        printf("Memory allocation failed for constraints.\n");
        exit(1);
    }

    // --- Upper and lower bounds constraints ---
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            int idx = varindex(i, j);

            // x_ij - upper_bound * y_ij <= 0
            cbeg[idx] = 2 * idx;
            cind[2 * idx] = M * N + idx;    cval[2 * idx] = 1.0;
            cind[2 * idx + 1] = idx;        cval[2 * idx + 1] = -upper_bound[i][j];
            sense[idx] = 'L';
            rhs[idx] = 0;
            sprintf(constrnames[idx], "upper_bound(%d,%d)", i, j);

            // lower_bound * y_ij - x_ij <= 0
            int idx2 = M * N + idx;
            cbeg[idx2] = 2 * M * N + 2 * idx;
            cind[2 * M * N + 2 * idx] = idx; cval[2 * M * N + 2 * idx] = lower_bound[i][j];
            cind[2 * M * N + 2 * idx + 1] = M * N + idx; cval[2 * M * N + 2 * idx + 1] = -1.0;
            sense[idx2] = 'L';
            rhs[idx2] = 0;
            sprintf(constrnames[idx2], "lower_bound(%d,%d)", i, j);
        }
    }

    // --- Assignment constraints sum_i y_ij = 1 ---
    for (int j = 0; j < N; j++) {
        int row = 2 * N * M + j;
        cbeg[row] = 4 * N * M + j * M;
        for (int i = 0; i < M; i++) {
            cind[4 * N * M + j * M + i] = varindex(i, j);
            cval[4 * N * M + j * M + i] = 1.0;
        }
        rhs[row] = 1;
        sense[row] = 'E';
        sprintf(constrnames[row], "assign(%d)", j);
    }

    // --- Capacity constraints sum_j x_ij <= capacity_i ---
    for (int i = 0; i < M; i++) {
        int row = 2 * N * M + N + i;
        cbeg[row] = 5 * N * M + N * i;
        for (int j = 0; j < N; j++) {
            cind[5 * N * M + N * i + j] = M * N + varindex(i, j);
            cval[5 * N * M + N * i + j] = 1.0;
        }
        rhs[row] = capacity[i];
        sense[row] = 'L';
        sprintf(constrnames[row], "capacity(%d)", i);
    }

    //Add new linear constraints to a model.
    error = GRBaddconstrs(model, numconstrs, numnz, cbeg, cind, cval, sense, rhs, constrnames);
    if (error) goto QUIT;

    //Free allocated memory for constraints
    free(rhs);
    free(sense);
    free(cbeg);
    free(cind);
    free(cval);
    free_stringarray(constrnames, char_array_size);

    if (P == 1.5) {
        // ==========================
        // add w variables
        // ==========================
        int num_new_cols = M * N;
        double* new_obj = create_double_vector(num_new_cols);
        double* new_lb = create_double_vector(num_new_cols);
        double* new_ub = create_double_vector(num_new_cols);
        char* new_xctype = create_char_vector(num_new_cols);
        char** new_colname = create_stringarray(num_new_cols, char_array_size);

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {

                int idx_w = varindex(i, j);

                new_obj[idx_w] = 0.0;
                new_lb[idx_w] = 0.0;
                new_ub[idx_w] = GRB_INFINITY;
                new_xctype[idx_w] = GRB_CONTINUOUS;
                sprintf(new_colname[idx_w], "w_(%d,%d)", i, j);
            }
        }
        error = GRBaddvars(model, num_new_cols, 0, NULL, NULL, NULL, new_obj, new_lb, new_ub, new_xctype, new_colname);
        if (error) goto QUIT;

        free(new_obj); free(new_lb); free(new_ub); free(new_xctype); free_stringarray(new_colname, char_array_size);

        // ==========================
        // SOC constraints
        // w^2 <= x*y, y^2 <= w*z
        // ==========================
        char qsense = GRB_LESS_EQUAL;
        double qrhs = 0.0;
        int qnumnz = 2;
        int qindrow[2], qindcol[2];
        double qval[2] = { 1.0, -1.0 };
        char qconstrname[char_array_size];

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {

                int idx_y = varindex(i, j);
                int idx_x = M * N + idx_y;
                int idx_z = 2 * M * N + idx_y;
                int idx_w = 3 * M * N + idx_y;

                // w^2 - x*y <= 0
                qindrow[0] = idx_w; qindcol[0] = idx_w;   // w^2
                qindrow[1] = idx_x; qindcol[1] = idx_y;   // -x*y
                sprintf(qconstrname, "socp1(%d,%d)", i, j);
                //Add a new quadratic constraint to a model.
                error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                if (error) goto QUIT;


                // y^2 - w*z <= 0
                qindrow[0] = idx_y; qindcol[0] = idx_y;   // y^2
                qindrow[1] = idx_w; qindcol[1] = idx_z;   // -w*z
                sprintf(qconstrname, "socp2(%d,%d)", i, j);
                //Add a new quadratic constraint to a model.
                error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                if (error) goto QUIT;
            }
        }
    }


    if (P == 2.0) {
        // ==========================
        // Quadratic constraints y^2 <= x * z
        // ==========================
        char qsense = GRB_LESS_EQUAL;
        double qrhs = 0.0;
        int qnumnz = 2;
        int qindrow[2], qindcol[2];
        double qval[2] = { 1.0, -1.0 };
        char qconstrname[char_array_size];

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                int idx_y = varindex(i, j);
                int idx_x = M * N + idx_y;
                int idx_z = 2 * M * N + idx_y;
                qindrow[0] = idx_y; qindcol[0] = idx_y; // y^2
                qindrow[1] = idx_x; qindcol[1] = idx_z; // -xz
                sprintf(qconstrname, "socp(%d,%d)", i, j);
                //Add a new quadratic constraint to a model.
                error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                if (error) goto QUIT;
            }
        }
    }

    if (P == 2.5) {
        // ==========================
        // add w, u, v variables
        // ==========================
        int num_new_cols = 3 * M * N;  // w, u, v
        double* new_obj = create_double_vector(num_new_cols);
        double* new_lb = create_double_vector(num_new_cols);
        double* new_ub = create_double_vector(num_new_cols);
        char* new_xctype = create_char_vector(num_new_cols);
        char** new_colname = create_stringarray(num_new_cols, char_array_size);

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {

                int idx_w = varindex(i, j);
                int idx_u = M * N + varindex(i, j);
                int idx_v = 2 * M * N + varindex(i, j);

                // w_ij
                new_obj[idx_w] = 0.0;
                new_lb[idx_w] = 0.0;
                new_ub[idx_w] = GRB_INFINITY;
                new_xctype[idx_w] = GRB_CONTINUOUS;
                sprintf(new_colname[idx_w], "w_(%d,%d)", i, j);

                // u_ij
                new_obj[idx_u] = 0.0;
                new_lb[idx_u] = 0.0;
                new_ub[idx_u] = GRB_INFINITY;
                new_xctype[idx_u] = GRB_CONTINUOUS;
                sprintf(new_colname[idx_u], "u_(%d,%d)", i, j);

                // v_ij
                new_obj[idx_v] = 0.0;
                new_lb[idx_v] = 0.0;
                new_ub[idx_v] = GRB_INFINITY;
                new_xctype[idx_v] = GRB_CONTINUOUS;
                sprintf(new_colname[idx_v], "v_(%d,%d)", i, j);
            }
        }

        error = GRBaddvars(model, num_new_cols, 0, NULL, NULL, NULL, new_obj, new_lb, new_ub, new_xctype, new_colname);
        if (error) goto QUIT;

        free(new_obj); free(new_lb); free(new_ub); free(new_xctype); free_stringarray(new_colname, char_array_size);

        // ==========================
        // Add SOC constraints
        // w^2 ≤ x*y, u^2 ≤ x*w, v^2 ≤ y*z, y^2 ≤ u*v
        // ==========================
        char qsense = GRB_LESS_EQUAL;
        double qrhs = 0.0;
        int qnumnz = 2;
        int qindrow[2], qindcol[2];
        double qval[2] = { 1.0, -1.0 };
        char qconstrname[char_array_size];

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {

                int idx_y = varindex(i, j);
                int idx_x = M * N + idx_y;
                int idx_z = 2 * M * N + idx_y;

                int idx_w = 3 * M * N + idx_y;
                int idx_u = 4 * M * N + idx_y;
                int idx_v = 5 * M * N + idx_y;

                // w^2 ≤ x*y
                qindrow[0] = idx_w; qindcol[0] = idx_w;
                qindrow[1] = idx_x; qindcol[1] = idx_y;
                sprintf(qconstrname, "socp1(%d,%d)", i, j);
                //Add a new quadratic constraint to a model.
                error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                if (error) goto QUIT;

                // u^2 ≤ x*w
                qindrow[0] = idx_u; qindcol[0] = idx_u;
                qindrow[1] = idx_x; qindcol[1] = idx_w;
                sprintf(qconstrname, "socp2(%d,%d)", i, j);
                //Add a new quadratic constraint to a model.
                error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                if (error) goto QUIT;

                // v^2 ≤ y*z
                qindrow[0] = idx_v; qindcol[0] = idx_v;
                qindrow[1] = idx_y; qindcol[1] = idx_z;
                sprintf(qconstrname, "socp3(%d,%d)", i, j);
                //Add a new quadratic constraint to a model.
                error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                if (error) goto QUIT;

                // y^2 ≤ u*v
                qindrow[0] = idx_y; qindcol[0] = idx_y;
                qindrow[1] = idx_u; qindcol[1] = idx_v;
                sprintf(qconstrname, "socp4(%d,%d)", i, j);
                //Add a new quadratic constraint to a model.
                error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                if (error) goto QUIT;
            }
        }
    }


    if (P == 3.0) {
        // ==========================
        // add w variables
        // ==========================
        int num_new_cols = M * N;
        double* new_obj = create_double_vector(num_new_cols);
        double* new_lb = create_double_vector(num_new_cols);
        double* new_ub = create_double_vector(num_new_cols);
        char* new_xctype = create_char_vector(num_new_cols);
        char** new_colname = create_stringarray(num_new_cols, char_array_size);

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {

                int idx_w = varindex(i, j);

                new_obj[idx_w] = 0.0;
                new_lb[idx_w] = 0.0;
                new_ub[idx_w] = GRB_INFINITY;
                new_xctype[idx_w] = GRB_CONTINUOUS;
                sprintf(new_colname[idx_w], "w_(%d,%d)", i, j);
            }
        }

        error = GRBaddvars(model, num_new_cols, 0, NULL, NULL, NULL, new_obj, new_lb, new_ub, new_xctype, new_colname);
        if (error) goto QUIT;

        free(new_obj); free(new_lb); free(new_ub); free(new_xctype); free_stringarray(new_colname, char_array_size);

        // ==========================
        // SOC constraints
        // w^2 <= y*z, y^2 <= w*x
        // ==========================
        char qsense = GRB_LESS_EQUAL;
        double qrhs = 0.0;
        int qnumnz = 2;
        int qindrow[2], qindcol[2];
        double qval[2] = { 1.0, -1.0 };
        char qconstrname[char_array_size];

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {

                int idx_y = varindex(i, j);
                int idx_x = M * N + idx_y;
                int idx_z = 2 * M * N + idx_y;
                int idx_w = 3 * M * N + idx_y;

                // w^2 - y*z <= 0
                qindrow[0] = idx_w; qindcol[0] = idx_w;   // w^2
                qindrow[1] = idx_y; qindcol[1] = idx_z;   // -y*z
                sprintf(qconstrname, "socp1(%d,%d)", i, j);
                //Add a new quadratic constraint to a model.
                error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                if (error) goto QUIT;

                // y^2 - w*x <= 0
                qindrow[0] = idx_y; qindcol[0] = idx_y;   // y^2
                qindrow[1] = idx_w; qindcol[1] = idx_x;   // -w*x
                sprintf(qconstrname, "socp2(%d,%d)", i, j);
                //Add a new quadratic constraint to a model.
                error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                if (error) goto QUIT;
            }
        }
    }

    //hibrid case with different P values
    if (P == 0) {
        // ==========================
        // add w, u, v variables
        // ==========================
        int num_new_cols = 3 * M * N;  // w, u, v
        double* new_obj = create_double_vector(num_new_cols);
        double* new_lb = create_double_vector(num_new_cols);
        double* new_ub = create_double_vector(num_new_cols);
        char* new_xctype = create_char_vector(num_new_cols);
        char** new_colname = create_stringarray(num_new_cols, char_array_size);

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {

                int idx_w = varindex(i, j);
                int idx_u = M * N + varindex(i, j);
                int idx_v = 2 * M * N + varindex(i, j);

                // w_ij
                new_obj[idx_w] = 0.0;
                new_lb[idx_w] = 0.0;
                new_ub[idx_w] = GRB_INFINITY;
                new_xctype[idx_w] = GRB_CONTINUOUS;
                sprintf(new_colname[idx_w], "w_(%d,%d)", i, j);

                // u_ij
                new_obj[idx_u] = 0.0;
                new_lb[idx_u] = 0.0;
                new_ub[idx_u] = GRB_INFINITY;
                new_xctype[idx_u] = GRB_CONTINUOUS;
                sprintf(new_colname[idx_u], "u_(%d,%d)", i, j);

                // v_ij
                new_obj[idx_v] = 0.0;
                new_lb[idx_v] = 0.0;
                new_ub[idx_v] = GRB_INFINITY;
                new_xctype[idx_v] = GRB_CONTINUOUS;
                sprintf(new_colname[idx_v], "v_(%d,%d)", i, j);
            }
        }

        error = GRBaddvars(model, num_new_cols, 0, NULL, NULL, NULL, new_obj, new_lb, new_ub, new_xctype, new_colname);
        if (error) goto QUIT;

        free(new_obj); free(new_lb); free(new_ub); free(new_xctype); free_stringarray(new_colname, char_array_size);


        // ==========================
        // SOC constraints according to P_current
        // ==========================

        char qsense = GRB_LESS_EQUAL;
        double qrhs = 0.0;
        int qnumnz = 2;
        int qindrow[2], qindcol[2];
        double qval[2] = { 1.0, -1.0 };
        char qconstrname[char_array_size];

        double P_current;
        for (int i = 0; i < M; i++) {
            P_current = P_array[i % 4]; // different P for each machine
            for (int j = 0; j < N; j++) {
                int idx_y = varindex(i, j);
                int idx_x = M * N + idx_y;
                int idx_z = 2 * M * N + idx_y;
                int idx_w = 3 * M * N + idx_y;
                int idx_u = 3 * M * N + M * N + idx_y;
                int idx_v = 3 * M * N + 2 * M * N + idx_y;

                if (P_current == 1.5) {
                    // w^2 - x*y <= 0
                    qindrow[0] = idx_w; qindcol[0] = idx_w;   // w^2
                    qindrow[1] = idx_x; qindcol[1] = idx_y;   // -x*y
                    sprintf(qconstrname, "socp1(%d,%d)", i, j);
                    //Add a new quadratic constraint to a model.
                    error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                    if (error) goto QUIT;

                    // y^2 - w*z <= 0
                    qindrow[0] = idx_y; qindcol[0] = idx_y;   // y^2
                    qindrow[1] = idx_w; qindcol[1] = idx_z;   // -w*z
                    sprintf(qconstrname, "socp2(%d,%d)", i, j);
                    //Add a new quadratic constraint to a model.
                    error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                    if (error) goto QUIT;
                }
                else if (P_current == 2.0) {
                    // y^2 - x * z <= 0
                    qindrow[0] = idx_y; qindcol[0] = idx_y; // y^2
                    qindrow[1] = idx_x; qindcol[1] = idx_z; // -xz
                    sprintf(qconstrname, "socp(%d,%d)", i, j);
                    //Add a new quadratic constraint to a model.
                    error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                    if (error) goto QUIT;
                }
                else if (P_current == 2.5) {
                    // w^2 ≤ x*y
                    qindrow[0] = idx_w; qindcol[0] = idx_w;
                    qindrow[1] = idx_x; qindcol[1] = idx_y;
                    sprintf(qconstrname, "socp1(%d,%d)", i, j);
                    //Add a new quadratic constraint to a model.
                    error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                    if (error) goto QUIT;

                    // u^2 ≤ x*w
                    qindrow[0] = idx_u; qindcol[0] = idx_u;
                    qindrow[1] = idx_x; qindcol[1] = idx_w;
                    sprintf(qconstrname, "socp2(%d,%d)", i, j);
                    //Add a new quadratic constraint to a model.
                    error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                    if (error) goto QUIT;

                    // v^2 ≤ y*z
                    qindrow[0] = idx_v; qindcol[0] = idx_v;
                    qindrow[1] = idx_y; qindcol[1] = idx_z;
                    sprintf(qconstrname, "socp3(%d,%d)", i, j);
                    //Add a new quadratic constraint to a model.
                    error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                    if (error) goto QUIT;

                    // y^2 ≤ u*v
                    qindrow[0] = idx_y; qindcol[0] = idx_y;
                    qindrow[1] = idx_u; qindcol[1] = idx_v;
                    sprintf(qconstrname, "socp4(%d,%d)", i, j);
                    //Add a new quadratic constraint to a model.
                    error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                    if (error) goto QUIT;
                }
                else if (P_current == 3.0) {
                    // w^2 - y*z <= 0
                    qindrow[0] = idx_w; qindcol[0] = idx_w;   // w^2
                    qindrow[1] = idx_y; qindcol[1] = idx_z;   // -y*z
                    sprintf(qconstrname, "socp1(%d,%d)", i, j);
                    //Add a new quadratic constraint to a model.
                    error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                    if (error) goto QUIT;

                    // y^2 - w*x <= 0
                    qindrow[0] = idx_y; qindcol[0] = idx_y;   // y^2
                    qindrow[1] = idx_w; qindcol[1] = idx_x;   // -w*x
                    sprintf(qconstrname, "socp2(%d,%d)", i, j);
                    //Add a new quadratic constraint to a model.
                    error = GRBaddqconstr(model, 0, NULL, NULL, qnumnz, qindrow, qindcol, qval, qsense, qrhs, qconstrname);
                    if (error) goto QUIT;
                }
            }
        }
    }

    //// Write model to 'miscop.lp'
    //error = GRBwrite(model, "miscop.lp");
    //if (error) goto QUIT;

    // ==========================
    // Set solver parameters
    // ==========================
    
    clock_t start, end;
    double NodesRuntime, gap, nodecount;
    

    // ==========================
    // Solve root node relaxation
    // ==========================
    
	// Limit to 1 node
    error = GRBsetdblparam(GRBgetenv(model), "NodeLimit", 1);
	if (error) goto QUIT;
    // Set number of threads to 1
    error = GRBsetintparam(GRBgetenv(model), "Threads", 1);
    if (error) goto QUIT;
    // Set feasibility tolerance
    error = GRBsetdblparam(GRBgetenv(model), "FeasibilityTol", 1e-9);
    if (error) goto QUIT;
    // Set optimality tolerance
    error = GRBsetdblparam(GRBgetenv(model), "OptimalityTol", 1e-9);
    if (error) goto QUIT;
    // Set integer feasibility tolerance
    error = GRBsetdblparam(GRBgetenv(model), "IntFeasTol", 1e-9);
    if (error) goto QUIT;
	// Set time limit for root node relaxation
    error = GRBsetdblparam(GRBgetenv(model), "TimeLimit", 1000.0);
    if (error) goto QUIT;

    start = clock();
    /* Optimize model */
    error = GRBoptimize(model);
    if (error) goto QUIT;
    end = clock();
    RootRuntime = ((double)(end - start)) / CLOCKS_PER_SEC;

    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &UpperBound);
    if (error) goto QUIT;
    error = GRBgetdblattr(model, "ObjBoundC", &LowerBound);   // 当前最好下界
    if (error) goto QUIT;
    error = GRBgetdblattr(model, GRB_DBL_ATTR_MIPGAP, &gap);
	if (error) goto QUIT;
	printf("\nRoot Node Relaxation complete!\n");
	printf("RootLowerBound:%.6f;RootUpperBound:%.6f;RootGap:%.6f;RootRuntime:%.4f;\n",
        LowerBound, UpperBound, gap, RootRuntime);
    fprintf(output, "RootLowerBound:%.6f;RootUpperBound:%.6f;RootGap:%.6f;RootRuntime:%.4f;\n",
        LowerBound, UpperBound, gap, RootRuntime);

    if (RootRuntime < 1000.0) {

        // ==========================
        // Solve full MIP
        // ==========================    
        // Remove node limit
        error = GRBsetdblparam(GRBgetenv(model), "NodeLimit", GRB_INFINITY);
        if (error) goto QUIT;
        // Set time limit for full MIP solve
        error = GRBsetdblparam(GRBgetenv(model), "TimeLimit", 1000.0 - RootRuntime);
        if (error) goto QUIT;


        start = clock();
        /* Optimize model */
        error = GRBoptimize(model);
        if (error) goto QUIT;
        end = clock();
        NodesRuntime = ((double)(end - start)) / CLOCKS_PER_SEC;

        /* Capture solution information */
        int       optimstatus; 
        error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
        if (error) goto QUIT;

        printf("\nOptimization complete!\n");
        if (optimstatus == GRB_OPTIMAL) {
            // Get the objective value       
            error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &UpperBound);
            if (error) goto QUIT;
            GRBgetdblattr(model, "ObjBoundC", &LowerBound);   // 当前最好下界
            GRBgetdblattr(model, GRB_DBL_ATTR_MIPGAP, &gap);
            GRBgetdblattr(model, GRB_DBL_ATTR_NODECOUNT, &nodecount);


            printf("MIP_ObjVal: %f\n", UpperBound);
            printf("MIP_Gap: %f\n", gap);
            printf("MIP_Nodes: %f\n", nodecount);
            printf("MIP_Runtime: %f\n", NodesRuntime);

            fprintf(output, "LowerBound:%.6f;UpperBound:%.6f;Gap:%.6f;BranchTime:%.4f;Nodes:%f;Status:%d;\n",
                LowerBound, UpperBound, gap, NodesRuntime, nodecount, optimstatus);
        }
        else if (optimstatus == GRB_INF_OR_UNBD) {
            printf("Model is infeasible or unbounded\n");
        }
        else {
            printf("Optimization was stopped early\n");
            // Get the objective value       
            error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &UpperBound);
            if (error) goto QUIT;
            GRBgetdblattr(model, "ObjBoundC", &LowerBound);   // 当前最好下界
            GRBgetdblattr(model, GRB_DBL_ATTR_MIPGAP, &gap);
            GRBgetdblattr(model, GRB_DBL_ATTR_NODECOUNT, &nodecount);


            printf("MIP_ObjVal: %f\n", UpperBound);
            printf("MIP_Gap: %f\n", gap);
            printf("MIP_Nodes: %f\n", nodecount);
            printf("MIP_Runtime: %f\n", NodesRuntime);

            fprintf(output, "LowerBound:%.6f;UpperBound:%.6f;Gap:%.6f;BranchTime:%.4f;Nodes:%f;Status:%d;\n",
                LowerBound, UpperBound, gap, NodesRuntime, nodecount, optimstatus);
        }

        // Check if a feasible solution was found
        int solcount;
        error = GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &solcount);
        if (error) goto QUIT;

        if (solcount > 0) {
            printf("Found feasible solution.\n");
            // check the solution
            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, M * N, Separating_Point);
            if (error) goto QUIT;
            // Calculate objective value from fixed charges   
            double obj_val = 0.0;
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    obj_val += fixed_charge[i][j] * Separating_Point[i * N + j];
                }
            }
            // Solve the subproblem
            CutSeparator();
            obj_val = obj_val + theta_try_solution[0];
            printf("\nChecking Final Objective Value: %f\n", obj_val);
        }
        else {
            printf("No feasible solution found.\n");
        }
    }

QUIT:
    /* Error reporting */
    if (error) {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }

    /* Free model */
    GRBfreemodel(model);

    /* Free environment */
    GRBfreeenv(env);
}
