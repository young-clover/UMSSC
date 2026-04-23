#include "def.h"

typedef struct {
    int m1;
    int m2;
    int* lns_sub_tasks;
    int sub_task_count;
} CallbackData;

void CutSeparatorLNS(int m1, int m2, const int* lns_sub_tasks, int k) {
    double voptsum = 0, all_subgradient_sum = 0;

    // 清零当前子问题的割系数，并绑定正确的变量索引，注意留出空间给最后一个 theta
    for (int step = 0; step < k; step++) {
        benders_cut->val[step] = 0.0;
        benders_cut->ind[step] = step; 
    }
    benders_cut->ind[k] = k;
    benders_cut->val[k] = -1.0;
    // 只遍历给定的两台机器: m1 和 m2
    int machines[2] = { m1, m2 };

    for (int idx = 0; idx < 2; idx++) {
        int i = machines[idx];

        // Preprocessing 
        if (P == 0)
            P = P_array[i % 4]; // different P for each machine    

        double vopt_i = 0, sum_x = 0;

        // compute temp variables only for the sub-problem variables in [k]
        for (int step = 0; step < k; step++) {
            int j = lns_sub_tasks[step]; // 获取真实的任务索引 

            // 取决于是机器 m1 还是 m2 来确定 Separation Point
            double sp = Separating_Point[step];
            if (i == m2) {
                sp = 1.0 - sp; // m2 的参数取 1 - sp
            }

            L[step] = lower_bound[i][j] * sp;
            U[step] = upper_bound[i][j] * sp;
            W[step] = weight[i][j] * sp;
            sum_x += U[step];
        }

        // ad-hoc solver
        double lambda = 0;
        int active_count = 0;

        // if total sum greater than capacity, use ad-hoc solver
        if (sum_x > capacity[i]) {
            double tau = 0;
            for (int step = 0; step < k; step++) {
                int j = lns_sub_tasks[step];
                double sp = (i == m2) ? (1.0 - Separating_Point[step]) : Separating_Point[step];

                if (sp == 0) {
                    x[i * N + j] = 0;
                }
                else {
                    S[active_count++] = step; // 记录局部索引 step
                    tau += W[step];
                }
            }
            lambda = capacity[i] / tau;

            while (active_count > 0) {
                double delta1 = 0, delta2 = 0;
                // compute x_j and deltas
                for (int m = 0; m < active_count; m++) {
                    int step = S[m];
                    int j = lns_sub_tasks[step];
                    x[i * N + j] = lambda * W[step];

                    if (x[i * N + j] <= L[step])
                        delta1 += (L[step] - x[i * N + j]);
                    else if (x[i * N + j] >= U[step])
                        delta2 += (x[i * N + j] - U[step]);
                }

                // if deltas are equal, we have found the solution
                int new_count = 0;
                if (delta1 == delta2) {
                    for (int m = 0; m < active_count; m++) {
                        int step = S[m];
                        int j = lns_sub_tasks[step];
                        if (x[i * N + j] <= L[step]) {
                            x[i * N + j] = L[step];
                        }
                        else if (x[i * N + j] >= U[step]) {
                            x[i * N + j] = U[step];
                        }
                        else {
                            S[new_count++] = step; // retain active elements                           
                        }
                    }
                    break;
                }
                // adjust lambda and update active set
                else if (delta1 > delta2) {
                    for (int m = 0; m < active_count; m++) {
                        int step = S[m];
                        int j = lns_sub_tasks[step];
                        if (x[i * N + j] <= L[step]) {
                            x[i * N + j] = L[step]; // fixed lower bound                            
                            tau -= W[step];
                        }
                        else {
                            S[new_count++] = step; // retain active elements                      
                        }
                    }
                    lambda -= delta1 / tau; // update lambda					
                }
                else {
                    for (int m = 0; m < active_count; m++) {
                        int step = S[m];
                        int j = lns_sub_tasks[step];
                        if (x[i * N + j] >= U[step]) {
                            x[i * N + j] = U[step]; // fixed upper bound                            
                            tau -= W[step];
                        }
                        else {
                            S[new_count++] = step; // retain active elements
                        }
                    }
                    lambda += delta2 / tau; // update lambda
                }
                active_count = new_count;
            }
        }
        else {
            for (int step = 0; step < k; step++) {
                int j = lns_sub_tasks[step];
                x[i * N + j] = U[step];
            }
        }

        // recover lambda
        if (active_count > 0) {
            int step = S[0];
            int j = lns_sub_tasks[step];
            lambda = (P - 1) * pow(W[step] / x[i * N + j], P);
        }
        else {
            lambda = 0;
        }

        // compute vopt_i and subgradient       
        for (int step = 0; step < k; step++) {
            int j = lns_sub_tasks[step];
            double sp = (i == m2) ? (1.0 - Separating_Point[step]) : Separating_Point[step];
            int sign_sp = (i == m2) ? -1 : 1; // 对 m2 时求关于原 Separating_Point 的导数需乘 -1

            alpha[j] = 0;
            beta[j] = 0;

            if (sp == 0) {
                if (lambda >= 0)
                    alpha[j] = lambda;
                else
                    beta[j] = -lambda;

                subgradient[step] = lower_bound[i][j] * alpha[j] - upper_bound[i][j] * beta[j];
            }
            else {
                double tempValue = pow(W[step], P) * pow(x[i * N + j], 1 - P);
                vopt_i += tempValue;

                if (x[i * N + j] == L[step])
                    alpha[j] = (1.0 - P) * pow(weight[i][j] / lower_bound[i][j], P) + lambda;
                if (x[i * N + j] == U[step])
                    beta[j] = -lambda - (1.0 - P) * pow(weight[i][j] / upper_bound[i][j], P);

                subgradient[step] = P * tempValue / sp + lower_bound[i][j] * alpha[j] - upper_bound[i][j] * beta[j];
            }
            
			benders_cut->val[step] += sign_sp * subgradient[step];// 注意 m2 的 subgradient 需要乘 -1
            all_subgradient_sum += sign_sp * subgradient[step] * Separating_Point[step];

        }
        voptsum += vopt_i;
    }

    benders_cut->rhs[0] = all_subgradient_sum - voptsum;
    // 显式指定符号为 GRB_LESS_EQUAL (也就是 <=)
    benders_cut->sense[0] = GRB_LESS_EQUAL; 
}

/* Define my callback function for reduced model*/
int __stdcall subBendersCutCallback(GRBmodel* model, void* cbdata, int where, void* usrdata) {
    int error = 0;
    // 取出传入的参数
    CallbackData* d = (CallbackData*)usrdata;
    int m1 = d->m1;
    int m2 = d->m2;
    int* tasks = d->lns_sub_tasks;
    int count = d->sub_task_count;

    if (where == GRB_CB_MIPSOL) {
        /* MIP solution callback */
        error = GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, Separating_Point);
        for (int i = 0; i < count; i++) {
            if (Separating_Point[i] < 0) Separating_Point[i] = 0;
            if (Separating_Point[i] > 1) Separating_Point[i] = 1;
        }
        // Solve the subproblem
        CutSeparatorLNS(m1, m2, tasks, count);

        error = GRBcblazy(cbdata, count + 1, benders_cut->ind, benders_cut->val, benders_cut->sense[0], benders_cut->rhs[0]);

    }
    else if (where == GRB_CB_MIPNODE) {
        int status;
        /* MIP node callback */
        error = GRBcbget(cbdata, where, GRB_CB_MIPNODE_STATUS, &status);
        if (status == GRB_OPTIMAL) {
            error = GRBcbget(cbdata, where, GRB_CB_MIPNODE_REL, Separating_Point);
            for (int i = 0; i < count; i++) {
                if (Separating_Point[i] < 0) Separating_Point[i] = 0;
                if (Separating_Point[i] > 1) Separating_Point[i] = 1;
            }
            // Solve the subproblem
            CutSeparatorLNS(m1, m2, tasks, count);

            error = GRBcbcut(cbdata, count + 1, benders_cut->ind, benders_cut->val, benders_cut->sense[0], benders_cut->rhs[0]);

        }
    }
    return error;
}

double MKP(GRBenv* env, int m1, int m2, int* lns_sub_tasks, int sub_task_count) {
    GRBmodel* model = NULL;
    int error = 0;

    // 总变量数: y_j 的数量
    int num_vars = sub_task_count;

    // 用于创建变量的数组
    double* obj = (double*)malloc(num_vars * sizeof(double));
    char* vtype = (char*)malloc(num_vars * sizeof(char));
    double* lb = (double*)malloc(num_vars * sizeof(double));
    double* ub = (double*)malloc(num_vars * sizeof(double));

    // 用于创建约束的数组
    int* ind = (int*)malloc(sub_task_count * sizeof(int));
    double* val = (double*)malloc(sub_task_count * sizeof(double));

    if (!obj || !vtype || !lb || !ub || !ind || !val) {
        printf("Memory allocation failed!\n");
        goto QUIT;
    }

    double const_obj = 0.0;             // 目标函数的常数项 \sum f_{2j}
    double sum_x2 = 0.0;                // 用于计算约束2的右侧常量

    // 1. 初始化每个 y_j 变量
    for (int k = 0; k < sub_task_count; ++k) {
        int j = lns_sub_tasks[k];

        // 目标函数系数组项: (f_1j - f_2j)
        obj[k] = fixed_charge[m1][j] - fixed_charge[m2][j];
        vtype[k] = GRB_CONTINUOUS; // 设为连续变量
        lb[k] = 0.0;            // 下界为 0
        ub[k] = 1.0;            // 上界为 1

        const_obj += fixed_charge[m2][j];
        sum_x2 += lower_bound[m2][j];
    }
	
    // 3. 创建模型
    error = GRBnewmodel(env, &model, NULL, num_vars, obj, lb, ub, vtype, NULL);
    if (error) goto QUIT;

    // 4. 将常数项叠加进目标函数
    error = GRBsetdblattr(model, GRB_DBL_ATTR_OBJCON, const_obj);
    if (error) goto QUIT;

    // 设定目标是最小化
    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);
    if (error) goto QUIT;

    // 5. 添加约束 1: \sum_{j=1}^n \underline{x}_{1j} y_j <= c_1
    for (int k = 0; k < sub_task_count; ++k) {
        int j = lns_sub_tasks[k];
        ind[k] = k;                      // 模型中的 y_j 变量索引
        val[k] = lower_bound[m1][j];     // \underline{x}_{1j}
    }
    // 注意：theta不参与该约束，所以约束变量个数传入 sub_task_count 即可
    error = GRBaddconstr(model, sub_task_count, ind, val, GRB_LESS_EQUAL, capacity[m1], "Capacity_m1");
    if (error) goto QUIT;

    // 6. 添加约束 2: \sum_{j=1}^n \underline{x}_{2j} y_j >= \sum_{j=1}^n \underline{x}_{2j} - c_2
    for (int k = 0; k < sub_task_count; ++k) {
        int j = lns_sub_tasks[k];
        val[k] = lower_bound[m2][j];     // \underline{x}_{2j}
    }
    double rhs_c2 = sum_x2 - capacity[m2];

    error = GRBaddconstr(model, sub_task_count, ind, val, GRB_GREATER_EQUAL, rhs_c2, "Capacity_m2");
    if (error) goto QUIT;

    // ==========================
    // parameter settings
    // ==========================
    // Turn off Gurobi output
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    if (error) goto QUIT;
    // Set feasibility tolerance
    error = GRBsetdblparam(GRBgetenv(model), "FeasibilityTol", 1e-9);
    if (error) goto QUIT;
    // Set optimality tolerance
    error = GRBsetdblparam(GRBgetenv(model), "OptimalityTol", 1e-9);
    if (error) goto QUIT;
    // Set number of threads to 1
    error = GRBsetintparam(GRBgetenv(model), "Threads", 1);
    if (error) goto QUIT;

    // ----- 利用 best_try_solution 对 PBD 进行 Warm Start -----
    double* start_vals = (double*)malloc(sub_task_count * sizeof(double));
    for (int k = 0; k < sub_task_count; ++k) {
        int j = lns_sub_tasks[k];
        // 取解在 m1 上的分配结果作为 y_k 的初值 
        start_vals[k] = best_try_solution[m1 * N + j];
		Stabilizing_Point[k] = start_vals[k]; // 同时初始化稳定化点
    }
    error = GRBsetdblattrarray(model, GRB_DBL_ATTR_START, 0, sub_task_count, start_vals);
    if (error) goto QUIT;
    free(start_vals);    

    // Optimize the model
    error = GRBoptimize(model);
    if (error) goto QUIT;

    // Get the solution
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, sub_task_count, Y);
    if (error) goto QUIT;

    // Get the objective value
	double subLowerBound;
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &subLowerBound);
    if (error) goto QUIT;

    // Add single theta variable to master problem
    error = GRBaddvar(model, 0, NULL, NULL, 1.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS, "theta");
    if (error) goto QUIT;

    //配置 theta 的系数   
    benders_cut->val[sub_task_count] = -1.0;           // theta的系数为-1 

    //Stabilization parameters
    double sigma = 0.6, gamma = 0.9, terminate = 0, last_lowerbound;
    do {
        last_lowerbound = subLowerBound;

        for (int i = 0; i < sub_task_count; i++) {
            // update separating point
            Separating_Point[i] = (1 - sigma) * Y[i] + sigma * Stabilizing_Point[i];
            // update stabilizing point
            Stabilizing_Point[i] = (1 - gamma) * Y[i] + gamma * Stabilizing_Point[i];
            if (Separating_Point[i] < 0)
                Separating_Point[i] = 0.0;
            if (Separating_Point[i] > 1)
                Separating_Point[i] = 1.0;
        }

        // Solve the subproblem
        CutSeparatorLNS(m1, m2, lns_sub_tasks, sub_task_count);
        
        // Add Benders cut to the master problem 
        error = GRBaddconstr(model, sub_task_count + 1, benders_cut->ind, benders_cut->val, benders_cut->sense[0], benders_cut->rhs[0], NULL);
        if (error) goto QUIT;        

        // Optimize the model
        error = GRBoptimize(model);
        if (error) goto QUIT;
        
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, sub_task_count, Y);
        if (error) goto QUIT;

        // Get the objective value
        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &subLowerBound);
        if (error) goto QUIT;

        if (subLowerBound - last_lowerbound > 1e-3) {
            terminate = 0;
            last_lowerbound = subLowerBound;
        }
        else {
            terminate += 1;
        }

    } while (terminate <= 3);
    // end of stabilization

    // Set the variables to binary type    
    char* vtype_array = create_char_vector(sub_task_count);
    for (int i = 0; i < sub_task_count; i++) {
        vtype_array[i] = GRB_BINARY;
    }
    error = GRBsetcharattrarray(model, GRB_CHAR_ATTR_VTYPE, 0, sub_task_count, vtype_array);
    if (error) goto QUIT;
    free(vtype_array);
    
    //// write model to disk for checking
    //error = GRBwrite(model, "submaster.lp");
    //if (error) goto QUIT;

	// Set Time limit
	error = GRBsetdblparam(GRBgetenv(model), "TimeLimit", 3.0);
	if (error) goto QUIT; 
    // Set integer feasibility tolerance
    error = GRBsetdblparam(GRBgetenv(model), "IntFeasTol", 1e-9);
    if (error) goto QUIT;
    // Enable lazy constraints
    error = GRBsetintparam(GRBgetenv(model), "LazyConstraints", 1);
    if (error) goto QUIT;

    CallbackData data;
    data.m1 = m1;
    data.m2 = m2;
    data.lns_sub_tasks = lns_sub_tasks;  // 指针直接传
    data.sub_task_count = sub_task_count;

    // callback function to add Benders cuts    
    error = GRBsetcallbackfunc(model, subBendersCutCallback, (void*)&data);
    if (error) goto QUIT;

    // Optimize the model
    error = GRBoptimize(model);
    if (error) goto QUIT;

    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, sub_task_count, Y);
    if (error) goto QUIT;

    for (int i = 0; i < sub_task_count; i++) {
        // update separating point
        Separating_Point[i] = Y[i] ;
        // update stabilizing point
        Stabilizing_Point[i] = Y[i] ;
        if (Separating_Point[i] < 0)
            Separating_Point[i] = 0.0;
        if (Separating_Point[i] > 1)
            Separating_Point[i] = 1.0;
    }

	// Solve the subproblem update x and subgradient for the final solution
    CutSeparatorLNS(m1, m2, lns_sub_tasks, sub_task_count);

	double ObjVal;
    // Get the objective value
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &ObjVal);
    if (error) goto QUIT;

QUIT:
    if (obj) free(obj);
    if (vtype) free(vtype);
    if (lb) free(lb);
    if (ub) free(ub);
    if (ind) free(ind);
    if (val) free(val);

    /* Error reporting */
    if (error) {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }
    /* Free model */
    GRBfreemodel(model);

    return ObjVal;
}