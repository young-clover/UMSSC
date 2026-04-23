#include "def.h"

void PBD(void) {
    GRBenv* env = NULL;
    GRBmodel* model = NULL;
    int       error = 0;
	CONTINUE = 0; // reset global variable in case of multiple runs

	// Create environment
    error = GRBemptyenv(&env);
	if (error) goto QUIT;
    error = GRBsetintparam(env, "OutputFlag", 0);
	if (error) goto QUIT;
    error = GRBstartenv(env);
	if (error) goto QUIT;

    // Create new model 
    error = GRBnewmodel(env, &model, "MasterProblem", 0, NULL, NULL, NULL, NULL, NULL);
    if (error) goto QUIT;

    // ==========================
    // Build the master problem
    // ==========================
    MasterProblem(env, model);

	// ==========================
	// Solve the root node
	// ==========================
    RootNodeSolving(env, model);

	// =======================================================
	// Apply LocalSearch to get a high-quality feasible solution
	// =======================================================
    if (CONTINUE >= 1)
        LocalSearch();

    // ==========================
    // Apply reduced cost fixing
    // ==========================
    if (CONTINUE >= 1)
        ReducedCostFixing(env, model);

	// ==========================
	// Solve the MIP
	// ==========================
    if (CONTINUE >= 0)
        MIPSolving(env, model);

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

void MasterProblem(GRBenv* env, GRBmodel* model) {
    int       error = 0;
    int numvars = M * N;
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
            vtype[varindex(i, j)] = GRB_CONTINUOUS;
            sprintf(varnames[varindex(i, j)], "y_(%d,%d)", i, j);
        }
    }

    // Add variables to the model
    error = GRBaddvars(model, numvars, 0, NULL, NULL, NULL, obj, lb, ub, vtype, varnames);
    if (error) goto QUIT;

    // Free allocated memory for variables
    free(obj);
    free(lb);
    free(ub);
    free(vtype);
    free_stringarray(varnames, numvars);

    // Add linear constraints to the model
    int numconstrs = N + M;
    int numnz = 2 * N * M;
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
	// Assignment constraints: each task is assigned to exactly one machine
    for (int j = 0; j < N; j++) {
        cbeg[j] = j * M;
        for (int i = 0; i < M; i++) {
            cind[j * M + i] = varindex(i, j);
            cval[j * M + i] = 1.0;
        }
        rhs[j] = 1;
        sense[j] = GRB_EQUAL;
        sprintf(constrnames[j], "assign(%d)", j);
    }
	// Capacity constraints: total load on each machine does not exceed its capacity
    for (int i = 0; i < M; i++) {
        cbeg[N + i] = M * N + N * i;
        for (int j = 0; j < N; j++) {
            cind[M * N + N * i + j] = varindex(i, j);
            cval[M * N + N * i + j] = lower_bound[i][j];
        }
        rhs[N + i] = capacity[i];
        sense[N + i] = GRB_LESS_EQUAL;
        sprintf(constrnames[N + i], "cap(%d)", i);
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
    free_stringarray(constrnames, numconstrs);

QUIT:
    /* Error reporting */
    if (error) {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }
}



void RootNodeSolving(GRBenv* env, GRBmodel* model) {
    clock_t start, end;
    int error = 0;   

	// ==========================
	// parameter settings
	// ==========================   
	// Set feasibility tolerance
    error = GRBsetdblparam(GRBgetenv(model), "FeasibilityTol", 1e-9);
	if (error) goto QUIT;
	// Set optimality tolerance
    error = GRBsetdblparam(GRBgetenv(model), "OptimalityTol", 1e-9);
    if (error) goto QUIT;		
	// Set number of threads to 1
    error = GRBsetintparam(GRBgetenv(model), "Threads", 1);
	if (error) goto QUIT;


    start = clock();
    /* Change sense to maximization to construct the stabilizing point */
    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
    if (error) goto QUIT;

    // Optimize the model
    error = GRBoptimize(model);
    if (error) goto QUIT;

	// initialize the stabilizing point
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, M * N, Stabilizing_Point);
    if (error) goto QUIT;

    /* Change sense to minimization */
    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);
    if (error) goto QUIT;

    // Optimize the model
    error = GRBoptimize(model);
    if (error) goto QUIT;

    // Get the solution
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, M * N, Y);
    if (error) goto QUIT;

    // Get the objective value
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &LowerBound);
    if (error) goto QUIT;

    if (pbd_mode == MODE_SINGLE)
        printf("\n>>> Start Single-Cut Iteration <<<\n");
    else if (pbd_mode == MODE_MULTI)
        printf("\n>>> Start Multi-Cut Iteration <<<\n");
    // printf("iter(%d)-LP_Lowerbound: %.6f \n", 0, LowerBound);

	// Add theta variable(s) to the master problem
    if (pbd_mode == MODE_SINGLE){
        // Add single theta variable to master problem
        error = GRBaddvar(model, 0, NULL, NULL, 1.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS, "theta");
        if (error) goto QUIT;
	}	
    else if (pbd_mode == MODE_MULTI) {
        char new_colname[char_array_size];
        // add multiple theta variables to master problem objective: minimize sum(theta_i)       
        for (int i = 0; i < M; i++) {
            sprintf(new_colname, "theta_(%d)", i);
            error = GRBaddvar(model, 0, NULL, NULL, 1.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS, new_colname);
            if (error) goto QUIT;
        }
    }
	
	//Stabilization parameters
    double sigma = 0.6, gamma = 0.9, LpGap = 1.0;
      
    int iteration_count = 0;

    char cutname[char_array_size];
    do {
        for (int i = 0; i < M * N; i++){
		    // update separating point
            Separating_Point[i] = (1 - sigma) * Y[i] + sigma * Stabilizing_Point[i];
            // update stabilizing point
            Stabilizing_Point[i] = (1 - gamma) * Y[i] + gamma * Stabilizing_Point[i];  
            if (Separating_Point[i] < 0)
                Separating_Point[i] = 0.0;
            if (Separating_Point[i] > 1)
				Separating_Point[i] = 1.0;
        }

        // Get the new objective value to update lp upper bound
        LP_UpperBound = 0.0;
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                LP_UpperBound += fixed_charge[i][j] * Separating_Point[varindex(i, j)];

        // Solve the subproblem
        CutSeparator();
        //printf("iter(%d)-LP_UpperBound: %.6f; \n", iteration_count, LP_UpperBound);

        // print the lp gap
        LpGap = fabs(LP_UpperBound - LowerBound) / fabs(LP_UpperBound);
        //printf("iter(%d)-LP_Gap: %.4f;\n", iteration_count, LpGap);
        
        iteration_count++; 
        int total_nz = M * (N + 1);
		// Add Benders cut to the master problem
        if (pbd_mode == MODE_SINGLE) {
            sprintf(cutname, "cut(%d)", iteration_count);
            error = GRBaddconstr(model, M * N + 1, benders_cut->ind, benders_cut->val, benders_cut->sense[0], benders_cut->rhs[0], cutname);
            if (error) goto QUIT;
        }
        else if (pbd_mode == MODE_MULTI) {
            //for (int i = 0; i < M; i++) {
            //    sprintf(cutname, "cut_(%d)_(%d)", iteration_count, i);
            //    error = GRBaddconstr(model, N + 1, &(benders_cut->ind[i * (N + 1)]), &(benders_cut->val[i * (N + 1)]), benders_cut->sense[i], benders_cut->rhs[i], cutname);
            //    if (error) goto QUIT;
            //}
            error = GRBaddconstrs(
                model,
                M,               // 添加的约束数量
                total_nz,        // 非零元素总数
                benders_cut->beg,// 每条约束的起始位置索引数组
                benders_cut->ind,// 所有约束的列索引合集
                benders_cut->val,// 所有约束的系数合集
                benders_cut->sense, // 符号数组 (GRB_LESS_EQUAL)
                benders_cut->rhs,   // 右端项数组
                NULL            // 约束名称
            );
		}
        
        // Optimize the model
        error = GRBoptimize(model);
        if (error) goto QUIT;

        // Get the new solution
        if (pbd_mode == MODE_SINGLE) {
            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, M * N + 1, Y);
            if (error) goto QUIT;
        }
        else if (pbd_mode == MODE_MULTI) {
            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, M * N + M, Y);
            if (error) goto QUIT;
        }	

        // Get the Reduced Cost
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_RC, 0, M * N, ReducedCost);
        if (error) goto QUIT;
		// Get the dual variables for assignment constraints
        error = GRBgetdblattrarray(model, "Pi", 0, N, EqualityPi);
        if (error) goto QUIT;

        if (LpGap < 0.1) {
            // Apply Greedy Heuristic to update upper bound
			RunHeuristics();
        }
        
		if (M * N > 1000 && iteration_count % 5 == 0) {
            // Delete non-binding constraints
            int start_index, numconstrs;
            start_index = N + M;

            // 获取约束数量
            error = GRBgetintattr(model, "NumConstrs", &numconstrs);
            if (error) goto QUIT;

            double* pi = create_double_vector(numconstrs);
            int* del_indices = create_int_vector(numconstrs);
            if (pi == NULL || del_indices == NULL) {
    	        printf("Memory allocation failed for dual variables.\n");
    	        exit(1);
            }

            // 获取所有约束对应的对偶变量值 Pi
            error = GRBgetdblattrarray(model, "Pi", 0, numconstrs, pi);
            if (error) goto QUIT;

            int del_count = 0;
            for (int i = start_index; i < numconstrs; i++) {
    	        if (pi[i] == 0) {
    		        del_indices[del_count++] = i;
    	        }
            }

            if (del_count > 0) {
    	        error = GRBdelconstrs(model, del_count, del_indices);
    	        if (error) goto QUIT;
    	        // printf("\nRemove %d cuts!\n", del_count);
            }
            // Free allocated memory
            free(pi);
            free(del_indices);
		}       
        
		// Get the new objective value to update lower bound
        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &LowerBound);
        if (error) goto QUIT;
        // printf("iter(%d)-LP_Lowerbound: %.6f \n", iteration_count, LowerBound);  	

        // Get the new objective value to update lp upper bound
        LP_UpperBound = 0.0;
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                LP_UpperBound += fixed_charge[i][j] * Y[varindex(i, j)];

    } while (LpGap >= 1e-5 && iteration_count <= 1000);
	// end of stabilization

    end = clock();
    
    //// write model to disk for checking
    //error = GRBwrite(model, "master.lp");
    //if (error) goto QUIT;

    // record root node results
    RootRuntime = ((double)end - (double)start) / CLOCKS_PER_SEC;
    printf("RootRuntime: %f\n", RootRuntime);

    int numCuts;
    // 获取约束数量
    error = GRBgetintattr(model, "NumConstrs", &numCuts);
    if (error) goto QUIT;
    numCuts = numCuts - N - M; // 减去原始约束数量，得到添加的 Benders cut 数量

	Gap = fabs(UpperBound - LowerBound) / fabs(UpperBound);
	printf("RootUpperBound: %.6f; RootLowerBound: %.6f; RootGap: %.6f;\n", UpperBound, LowerBound, Gap);
    fprintf(output, "RootLowerBound:%.6f;RootUpperBound:%.6f;RootGap:%.6f;RootRuntime:%.4f;NumCuts:%d;NumIters:%d;\n", LowerBound, UpperBound, Gap, RootRuntime, numCuts, iteration_count);    

    if (Gap < 0.0001)
		CONTINUE = -1; // 如果根节点 gap 已经非常小了，就不继续后续的 MIP 求解了
QUIT:
    /* Error reporting */
    if (error) {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }

}


// Local Search 辅助结构体与排序函数
typedef struct {
    int id;
    double cost;
} MachineInfo;

int cmp_machine_info(const void* a, const void* b) {
    double diff = ((MachineInfo*)a)->cost - ((MachineInfo*)b)->cost;
    if (diff < 0) return -1;
    if (diff > 0) return 1;
    return 0;
}

void LocalSearch() {
    // 全局启发式计时起点
    clock_t start, end;
    start = clock();

    // Local Search 总预算：最多 500s，且不超过剩余总时限
    double total_remaining_time = 1000.0 - RootRuntime;
    const double TIME_LIMIT_SECONDS = (total_remaining_time < 500.0) ? total_remaining_time : 500.0;

    // 缓存当前最优连续变量解x
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            Separating_Point[i * N + j] = best_try_solution[i * N + j];
            if (Separating_Point[i * N + j] < 0)
                Separating_Point[i * N + j] = 0.0;
            if (Separating_Point[i * N + j] > 1)
                Separating_Point[i * N + j] = 1.0;
        }
    }
    CutSeparator();

    // ----------------------------------------------
    // Local Search: Systematic 2-Machine Improvement
    // ----------------------------------------------
    printf("\n>>> Start Local Search <<<\n");

    // ==== 创建全局静默环境 ====
    GRBenv* ls_env2 = NULL;
    GRBemptyenv(&ls_env2);
    GRBsetintparam(ls_env2, "OutputFlag", 0);
    // 设置1个线程以避免多核频繁上下文切换
    GRBsetintparam(ls_env2, "Threads", 1);
    GRBstartenv(ls_env2);

    int** banned_pair2 = (int**)malloc(M * sizeof(int*));
    MachineInfo* m_info2 = (MachineInfo*)malloc(M * sizeof(MachineInfo));

    if (banned_pair2 == NULL || m_info2 == NULL) {
        printf("Memory allocation failed in LocalSearch.\n");
        if (ls_env2) GRBfreeenv(ls_env2);
        if (banned_pair2) free(banned_pair2);
        if (m_info2) free(m_info2);
        exit(1);
    }

    // 计算当前最优解中每台机器确切的成本并初始化禁忌表
    for (int i = 0; i < M; i++) {
        banned_pair2[i] = (int*)calloc(M, sizeof(int));
        if (banned_pair2[i] == NULL) {
            printf("Memory allocation failed for banned_pair2 row.\n");
            for (int k = 0; k < i; k++) free(banned_pair2[k]);
            free(banned_pair2);
            free(m_info2);
            if (ls_env2) GRBfreeenv(ls_env2);
            exit(1);
        }

        m_info2[i].id = i;
        m_info2[i].cost = 0.0;

        double P_val = (P == 0) ? P_array[i % 4] : P;
        for (int j = 0; j < N; j++) {
            if (best_try_solution[i * N + j] > 0.5) {
                m_info2[i].cost += fixed_charge[i][j] + pow(weight[i][j], P_val) * pow(x[i * N + j], 1 - P_val);
            }
        }
    }

    int iter_count2 = 0;
    int max_iters2 = 1000 * M;
    int* lns_sub_tasks2 = (int*)malloc(N * sizeof(int));
    if (lns_sub_tasks2 == NULL) {
        printf("Memory allocation failed for lns_sub_tasks2.\n");
        for (int i = 0; i < M; i++) free(banned_pair2[i]);
        free(banned_pair2);
        free(m_info2);
        if (ls_env2) GRBfreeenv(ls_env2);
        exit(1);
    }

    // 循环外部仅排序一次，发生更新时才重新排序
    qsort(m_info2, M, sizeof(MachineInfo), cmp_machine_info);

    // 定义外部游标保存查找状态，避免重复O(M^2)扫描
    int current_len = M - 1;
    int current_left = 0;

    while (iter_count2++ < max_iters2) {
        double elapsed_seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
        if (elapsed_seconds >= TIME_LIMIT_SECONDS) break;

        int m1 = -1, m2 = -1;
        int found_pair = 0;

        // 使用游标继续寻找差距最大的未禁忌机器对
        while (current_len > 0 && !found_pair) {
            while (current_left < M - current_len && !found_pair) {
                int id1 = m_info2[current_left].id;
                int id2 = m_info2[current_left + current_len].id;
                if (!banned_pair2[id1][id2]) {
                    m1 = id1;
                    m2 = id2;
                    found_pair = 1;
                }
                if (!found_pair) {
                    current_left++; // 步进游标，向后扫描
                }
            }
            if (!found_pair) {
                current_len--;   // 降低组合跨度
                current_left = 0; // 重置起始点
            }
        }

        // 无更多可用机器对跳出
        if (!found_pair) break;

        // 提取分配给这两台机器上的任务
        int sub_task_count2 = 0;
        for (int j = 0; j < N; j++) {
            if (best_try_solution[m1 * N + j] > 0.5 || best_try_solution[m2 * N + j] > 0.5) {
                lns_sub_tasks2[sub_task_count2++] = j;
            }
        }

        if (sub_task_count2 <= 1) {
            current_left++; // 步进，忽略该对            
            continue;
        }

        // 获取更新前的老机器整体成本
        double old_cost_m1 = 0.0, old_cost_m2 = 0.0;
        int idx_m1 = -1, idx_m2 = -1;
        for (int k = 0; k < M; k++) {
            if (m_info2[k].id == m1) { old_cost_m1 = m_info2[k].cost; idx_m1 = k; }
            if (m_info2[k].id == m2) { old_cost_m2 = m_info2[k].cost; idx_m2 = k; }
        }
        if (idx_m1 < 0 || idx_m2 < 0) {
            current_left++;
            continue;
        }

        double old_obj_pair_val = old_cost_m1 + old_cost_m2;
        double new_obj_pair = MKP(ls_env2, m1, m2, lns_sub_tasks2, sub_task_count2);
        // 规则1：该机器对一旦被评估到，立即加入禁忌（无论后续是否改进）
        banned_pair2[m1][m2] = 1;
        banned_pair2[m2][m1] = 1;

        // 判断是否有提升
        if (new_obj_pair < old_obj_pair_val - 0.1) {
            // 根据 Y 提取结果，更新当前最优解二分类分配
            for (int k = 0; k < sub_task_count2; k++) {
                int j = lns_sub_tasks2[k];
                // Y 为 MKP 调用的最优二元分配结果
                if (Y[k] > 0.5) { // 重新分配给 m1
                    best_try_solution[m1 * N + j] = 1.0;
                    best_try_solution[m2 * N + j] = 0.0;
                }
                else { // 重新分配给 m2
                    best_try_solution[m1 * N + j] = 0.0;
                    best_try_solution[m2 * N + j] = 1.0;
                }
            }

            // 直接根据全局变量best_try_solution和x显式计算最新成本
            double new_cost_m1 = 0.0, new_cost_m2 = 0.0;
            double P_val1 = (P == 0) ? P_array[m1 % 4] : P;
            double P_val2 = (P == 0) ? P_array[m2 % 4] : P;
            for (int j = 0; j < N; j++) {
                if (best_try_solution[m1 * N + j] > 0.5) {
                    new_cost_m1 += fixed_charge[m1][j] + pow(weight[m1][j], P_val1) * pow(x[m1 * N + j], 1 - P_val1);
                }
                if (best_try_solution[m2 * N + j] > 0.5) {
                    new_cost_m2 += fixed_charge[m2][j] + pow(weight[m2][j], P_val2) * pow(x[m2 * N + j], 1 - P_val2);
                }
            }

            m_info2[idx_m1].cost = new_cost_m1;
            m_info2[idx_m2].cost = new_cost_m2;

            UpperBound = UpperBound - old_obj_pair_val + new_obj_pair;
            //printf("LocalSearch iter(%d): Improved pair (%d, %d) with cost %.6f -> %.6f; New UpperBound: %.6f;\n", iter_count2, m1, m2, old_obj_pair_val, new_obj_pair, UpperBound);

            // 规则2：若改进，仅解禁 m1/m2 与其它机器的组合；当前对(m1,m2)保持禁忌
            for (int i = 0; i < M; i++) {
                if (i != m1 && i != m2) {
                    banned_pair2[m1][i] = 0; banned_pair2[i][m1] = 0;
                    banned_pair2[m2][i] = 0; banned_pair2[i][m2] = 0;
                }
            }

            // 发生了确切提升：重新排序并重置游标
            qsort(m_info2, M, sizeof(MachineInfo), cmp_machine_info);
            current_len = M - 1;
            current_left = 0;
        }
        else {
            // 已在本轮开头禁忌该对，这里仅步进游标
            current_left++;
        }
    }

    // 释放分配的内存
    for (int i = 0; i < M; i++)
        free(banned_pair2[i]);
    free(banned_pair2);
    free(m_info2);
    free(lns_sub_tasks2);

    // ==== 结束后再释放环境 ====
    if (ls_env2) GRBfreeenv(ls_env2);

    end = clock();
    HeurRuntime = (double)(end - start) / CLOCKS_PER_SEC;
    printf("HeurRuntime: %.6f\n", HeurRuntime);

    // 获取其成本
    double obj_val = 0.0;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            Separating_Point[i * N + j] = best_try_solution[i * N + j]; // pass data for separating cut
            obj_val += fixed_charge[i][j] * best_try_solution[i * N + j];
            if (Separating_Point[i * N + j] < 0)
                Separating_Point[i * N + j] = 0.0;
            if (Separating_Point[i * N + j] > 1)
                Separating_Point[i * N + j] = 1.0;
        }
    }

    // Solve the subproblem
    CutSeparator();

    // update best solution after local search
    if (pbd_mode == MODE_SINGLE) {
        best_try_solution[M * N] = theta_try_solution[0];
        // final objective value
        obj_val += theta_try_solution[0];
    }
    else if (pbd_mode == MODE_MULTI) {
        for (int i = 0; i < M; i++) {
            best_try_solution[M * N + i] = theta_try_solution[i];
            // final objective value
            obj_val += theta_try_solution[i];
        }
    }

    UpperBound = obj_val;
    Gap = fabs(UpperBound - LowerBound) / fabs(UpperBound);
    printf("RootUpperBound: %.6f; RootLowerBound: %.6f; HeurGap: %.6f;\n", UpperBound, LowerBound, Gap);
    printf("Checking Final Objective Value: %f\n", obj_val);

    // 最终写文件：启发式时间
    fprintf(output, "HeurUpperBound:%.6f;HeurGap:%.6f;HeurRuntime:%.4f;\n", UpperBound, Gap, HeurRuntime);

    if (Gap < 0.0001)
        CONTINUE = -1; // 如果根节点 gap 已经非常小了，就不继续后续的 MIP 求解了

    if ((LowerBound - UpperBound) > 1e-6) {
        CONTINUE = -1;
        printf("\nSomething Wrong !\n");
        fprintf(output, "Something Wrong !\n");
    }
}


void ReducedCostFixing(GRBenv* env, GRBmodel* model) {
    clock_t start, end;
    int error = 0;
	//防御编程 处理数值误差导致的上下界过于接近的情况
    if ( (LowerBound - UpperBound) > 0 && (LowerBound - UpperBound) < 1e-6)
		UpperBound = LowerBound;
    printf("\n>>> Start Variable Fixing <<<\n");
	start = clock();

    // load start solution
    if (pbd_mode == MODE_SINGLE) {
        error = GRBsetdblattrarray(model, "Start", 0, M * N + 1, best_try_solution);
        if (error) goto QUIT;
    }
    else if (pbd_mode == MODE_MULTI) {
        error = GRBsetdblattrarray(model, "Start", 0, M * N + M, best_try_solution);
        if (error) goto QUIT;
    }
    
    

    //reduced_cost_fixing
    int fixed_sum = 0;
    for (int i = 0; i < M * N; i++) {
        if (ReducedCost[i] > UpperBound - LowerBound) {
            error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, i, 0);  // 设置上界值为0
            if (error) goto QUIT;
            fixed_sum += 1;
            fixed_point[i] = 0;//0代表固定为0
        }
        if (ReducedCost[i] < LowerBound - UpperBound) {
            error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, i, 1);  // 设置下界值为1
            if (error) goto QUIT;
            fixed_sum += 1;
            fixed_point[i] = 1;//1代表固定为1
        }
    }
    int num_need_more_try = 0;
    for (int i = 0; i < M * N; i++) {
        if (fixed_point[i] == -1 && (Y[i] < 0.1 || Y[i] > 0.9)) {
            num_need_more_try++;
        }
    }
    // try to fix more variables based on rlexation solution
    if (num_need_more_try < 5000) {
        int optimstatus;
        double lowerboud_try;
        for (int i = 0; i < M * N; i++) {
            if (fixed_point[i] == -1) {//从没有固定的里面下手
                if (Y[i] < 0.1) {
                    error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, i, 1);  // 尝试固定为1
                    if (error) goto QUIT;

                    error = GRBoptimize(model);
                    if (error) goto QUIT;
                    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
                    if (error) goto QUIT;
                    if (optimstatus == GRB_OPTIMAL) {
                        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &lowerboud_try);
                        if (error) goto QUIT;
                    }
                    else if (optimstatus == GRB_INFEASIBLE) {
                        lowerboud_try = UpperBound + 1;
                    }

                    if (lowerboud_try > UpperBound) {
                        error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, i, 0);  // 尝试成功上下界均固定为0
                        if (error) goto QUIT;
                        error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, i, 0);  // 尝试成功上下界均固定为0
                        if (error) goto QUIT;

                        fixed_point[i] = 0;//0代表固定为0
                        fixed_sum += 1;
                    }
                    else {
                        error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, i, 0);  //尝试失败 下界还原为0
                        if (error) goto QUIT;
                    }
                }
                if (Y[i] > 0.9) {
                    error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, i, 0);//尝试固定为0
                    if (error) goto QUIT;

                    error = GRBoptimize(model);
                    if (error) goto QUIT;
                    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
                    if (error) goto QUIT;
                    if (optimstatus == GRB_OPTIMAL) {
                        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &lowerboud_try);
                        if (error) goto QUIT;
                    }
                    else if (optimstatus == GRB_INFEASIBLE) {
                        lowerboud_try = UpperBound + 1;
                    }

                    if (lowerboud_try > UpperBound) {
                        error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, i, 1);  // 尝试成功上下界均固定为1
                        if (error) goto QUIT;
                        error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, i, 1);  // 尝试成功上下界均固定为1
                        if (error) goto QUIT;

                        fixed_point[i] = 1;//1代表固定为1
                        fixed_sum += 1;
                    }
                    else {
                        error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, i, 1);  //尝试失败 上界还原为1
                        if (error) goto QUIT;
                    }

                }
            }
        }
    }    
    printf("Fixing %d varaibles!\n", fixed_sum);

    // remove fixed variables from the model	
    if (fixed_sum > 0) {
        num_of_unfixed = num_of_unfixed - fixed_sum;
        if (num_of_unfixed == 0)
            printf("Eliminate all varaibles！\n");
        else {
            int* delvarind = create_int_vector(fixed_sum);
            if (delvarind == NULL) {
                printf("Memory allocation failed for delvarind.\n");
                exit(1);
            }


            for (int k1 = 0, k2 = 0, i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    if (fixed_point[i * N + j] >= 0) {
                        pre_fix_sum += fixed_charge[i][j] * fixed_point[i * N + j];
                        delvarind[k1] = i * N + j; // 记录需要删除的变量索引
                        k1++;
                    }
                    else {
                        unfixed_index[k2] = i * N + j;// 记录未固定变量的索引
                        k2++;
                    }
                }
            }

            //change the rhs of the constraints
            int numconstrs;
            error = GRBgetintattr(model, "NumConstrs", &numconstrs);
            if (error) goto QUIT;
            double coef, rhs;
            for (int i = 0; i < numconstrs; i++) {
                error = GRBgetdblattrelement(model, GRB_DBL_ATTR_RHS, i, &rhs);
                if (error) goto QUIT;
                for (int j = 0; j < M * N; j++) {
                    if (fixed_point[j] == 1.0) {
                        error = GRBgetcoeff(model, i, j, &coef);
                        if (error) goto QUIT;
                        rhs -= coef;
                    }
                }
                error = GRBsetdblattrelement(model, GRB_DBL_ATTR_RHS, i, rhs);
                if (error) goto QUIT;
            }

            //remove fixed variables from the model			
            error = GRBdelvars(model, fixed_sum, delvarind);
            if (error) goto QUIT;
            printf("Remove %d varaibles!\n", fixed_sum);
        }
    } 

    // update model
    error = GRBupdatemodel(model);
    if (error) goto QUIT;

    // adjust Benders cut coefficients for theta variables
    if (num_of_unfixed < M * N) {
        //printf("Current pre_fix_sum: %.6f\n", pre_fix_sum);
        //固定分离点中的变量
        for (int i = 0; i < M * N; i++) {
            if (fixed_point[i] >= 0) {
                Separating_Point[i] = fixed_point[i];
            }
        }

        // reset coefficients
        if (pbd_mode == MODE_SINGLE) {
            benders_cut->val[num_of_unfixed] = -1; // theta的系数为-1 
        }
        else if (pbd_mode == MODE_MULTI) {
            int unfixed_idx = 0;
            for (int i = 0; i < M; i++) { 
                int num_unfixed_in_row = 0;
                for (int j = 0; j < N; j++) {
                    if (fixed_point[i * N + j] == -1) {
                        benders_cut->ind[unfixed_idx] = unfixed_idx - i;
                        unfixed_idx++;
                        num_unfixed_in_row++;
                    }
                }
                benders_cut->ind[unfixed_idx] = num_of_unfixed + i;
                benders_cut->val[unfixed_idx] = -1; // coefficient for theta_i
                unfixed_idx++;
                num_unfixed_in_row++;
                // record number of coefficients per cut
                num_of_coeff_per_cut[i] = num_unfixed_in_row;
            }
        }

    }
	end = clock();
	FixingRuntime = (double)(end - start) / CLOCKS_PER_SEC;
    
	printf("FixingRuntime: %.6f\n", FixingRuntime);    
    fprintf(output, "FixingRuntime:%.4f;FixedNum:%d\n", FixingRuntime, fixed_sum);

QUIT:
    /* Error reporting */
    if (error) {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }
}

// Comparator function for sorting jobs by reduced cost difference in descending order
int cmpJobRcDescendingDiff(const void* a, const void* b) {// diff 大的先
    double da = ((JobRcDiff*)a)->diff;
    double db = ((JobRcDiff*)b)->diff;
    return (da < db) - (da > db); // 降序
}


void UpdateGlobalBestSolution(int* current_sol, const char* label) {
    double obj_val = 0.0;
    // 1. 设置 Separating_Point
    for (int i = 0; i < M * N; i++) {
        Separating_Point[i] = (double)current_sol[i];
        obj_val += fixed_charge[i / N][i % N] * current_sol[i];
    }

    // 2. 调用 Benders 子问题计算 theta
    CutSeparator();

    // 3. 累加 theta
    if (pbd_mode == MODE_SINGLE) obj_val += theta_try_solution[0];
    else for (int i = 0; i < M; i++) obj_val += theta_try_solution[i];

    // 4. 更新全局最优
    if (obj_val < UpperBound) {
        UpperBound = obj_val;
        CONTINUE = 1;
        for (int i = 0; i < M * N; i++) best_try_solution[i] = current_sol[i];

        int offset = (pbd_mode == MODE_SINGLE) ? 1 : M;
        for (int i = 0; i < offset; i++)
            best_try_solution[M * N + i] = theta_try_solution[i];

        // printf("[%s Heuristic] New UpperBound: %f\n", label, UpperBound);
    }
}

void Heuristic_Rounding(void) {
    double* remainingCap = (double*)malloc(M * sizeof(double));
    int* temp_sol = (int*)calloc(M * N, sizeof(int));
    int feasible = 1;

    for (int i = 0; i < M; i++) remainingCap[i] = capacity[i];

    for (int j = 0; j < N; j++) {
        int best_m = -1;
        double max_y = -1.0;
        // 找 LP 值最大的机器
        for (int i = 0; i < M; i++) {
            if (Y[i * N + j] > max_y) {
                max_y = Y[i * N + j];
                best_m = i;
            }
        }
        // 检查物理容量
        if (best_m != -1 && remainingCap[best_m] >= lower_bound[best_m][j]) {
            temp_sol[best_m * N + j] = 1;
            remainingCap[best_m] -= lower_bound[best_m][j];
        }
        else {
            feasible = 0; break;
        }
    }

    if (feasible) {
        UpdateGlobalBestSolution(temp_sol, "Rounding");
    }

    free(remainingCap); free(temp_sol);
}

void Heuristic_Regret(void) {
    double* remainingCap = (double*)malloc(M * sizeof(double));
    int* temp_sol = (int*)calloc(M * N, sizeof(int));
    int numSuccess = 0;

    for (int i = 0; i < M; i++) remainingCap[i] = capacity[i];

    // Step 1: 计算 Regret (min2 - min1)
    for (int j = 0; j < N; j++) {
        double m1 = GRB_INFINITY, m2 = GRB_INFINITY;
        for (int i = 0; i < M; i++) {
            double rc = ReducedCost[i * N + j];
            if (rc < m1) { m2 = m1; m1 = rc; }
            else if (rc < m2) { m2 = rc; }
        }
        jobs[j].job = j;
        jobs[j].diff = m2 - m1;
    }

    // Step 2: 降序排序 (后悔值大者优先)
    qsort(jobs, N, sizeof(JobRcDiff), cmpJobRcDescendingDiff);

    // Step 3: 分配
    for (int k = 0; k < N; k++) {
        int j = jobs[k].job;
        int best_m = -1; double best_rc = GRB_INFINITY;
        for (int i = 0; i < M; i++) {
            if (remainingCap[i] >= lower_bound[i][j] && ReducedCost[i * N + j] < best_rc) {
                best_rc = ReducedCost[i * N + j];
                best_m = i;
            }
        }
        if (best_m != -1) {
            temp_sol[best_m * N + j] = 1;
            remainingCap[best_m] -= lower_bound[best_m][j];
            numSuccess++;
        }
    }

    if (numSuccess == N) UpdateGlobalBestSolution(temp_sol, "Regret");

    free(remainingCap); free(temp_sol);
}

void RunHeuristics() {
    Heuristic_Rounding();
    Heuristic_Regret();
}

void CutSeparator(void) {
    double voptsum = 0, all_subgradient_sum = 0;
    int unfixed_idx = 0;
    for (int i = 0; i < M; i++) {
        // Preprocessing 
        if (P == 0)
            P = P_array[i % 4]; // different P for each machine    
        double vopt_i = 0, sum_x = 0;
        // compute temp variables
        for (int j = 0; j < N; j++) {
            double sp = Separating_Point[i * N + j];
            L[j] = lower_bound[i][j] * sp;
            U[j] = upper_bound[i][j] * sp;
            W[j] = weight[i][j] * sp;
            sum_x += U[j];
        }

        // ad-hoc solver
        double lambda = 0;
        int active_count = 0;
        // if total sum greater than capacity, use ad-hoc solver
        if (sum_x > capacity[i]) {
            double tau = 0;
            for (int j = 0; j < N; j++) {
                if (Separating_Point[i * N + j] == 0) {
                    x[i * N + j] = 0;
                }
                else {
                    S[active_count++] = j; // active set increment
                    tau += W[j];
                }
            }
            lambda = capacity[i] / tau;

            while (active_count > 0) {
                double delta1 = 0, delta2 = 0;
                // conpute x_j and deltas
                for (int k = 0; k < active_count; k++) {
                    int j = S[k];
                    x[i * N + j] = lambda * W[j];

                    if (x[i * N + j] <= L[j])
                        delta1 += (L[j] - x[i * N + j]);
                    else if (x[i * N + j] >= U[j])
                        delta2 += (x[i * N + j] - U[j]);
                }
                // if deltas are equal, we have found the solution
                int new_count = 0;
                if (delta1 == delta2) {
                    for (int k = 0; k < active_count; k++) {
                        int j = S[k];
                        if (x[i * N + j] <= L[j]) {
                            x[i * N + j] = L[j];
                        }
                        else if (x[i * N + j] >= U[j]) {
                            x[i * N + j] = U[j];
                        }
                        else {
                            S[new_count++] = j; // retain active elements                           
                        }
                    }
                    break;
                }
                // adjust lambda and update active set
                else if (delta1 > delta2) {
                    for (int k = 0; k < active_count; k++) {
                        int j = S[k];
                        if (x[i * N + j] <= L[j]) {
                            x[i * N + j] = L[j]; // fixed lower bound                            
                            tau -= W[j];
                        }
                        else {
                            S[new_count++] = j; // retain active elements                      
                        }
                    }
                    lambda -= delta1 / tau; // update lambda					
                }
                else {
                    for (int k = 0; k < active_count; k++) {
                        int j = S[k];
                        if (x[i * N + j] >= U[j]) {
                            x[i * N + j] = U[j]; // fixed upper bound                            
                            tau -= W[j];
                        }
                        else {
                            S[new_count++] = j; // retain active elements
                        }
                    }
                    lambda += delta2 / tau; // update lambda
                }
                active_count = new_count;
            }
            //printf("Ad-hoc solver used for subproblem %d\n", i);
        }
        else {
            for (int j = 0; j < N; j++) {
                x[i * N + j] = U[j];
            }
            //printf("Direct computation used for subproblem %d\n", i);
        }

        // recover lambda
        if (active_count > 0) {
            int j = S[0];
            lambda = (P - 1) * pow(W[j] / x[i * N + j], P);
        }
        else
            lambda = 0;

        // compute vopt_i and subgradient
		double i_subgradient_sum = 0;
        for (int j = 0; j < N; j++) {
            alpha[j] = 0;
            beta[j] = 0;
            if (Separating_Point[i * N + j] == 0) {
                if (lambda >= 0)
                    alpha[j] = lambda;
                else
                    beta[j] = -lambda;
                // 缓存子梯度
                subgradient[i * N + j] = lower_bound[i][j] * alpha[j] - upper_bound[i][j] * beta[j];
                
            }
            else {
                double tempValue = pow(W[j], P) * pow(x[i * N + j], 1 - P);				
                vopt_i += tempValue;

                if (x[i * N + j] == L[j])
                    alpha[j] = (1.0 - P) * pow(weight[i][j] / lower_bound[i][j], P) + lambda;
                if (x[i * N + j] == U[j])
                    beta[j] = -lambda - (1.0 - P) * pow(weight[i][j] / upper_bound[i][j], P);
                // 缓存子梯度
                subgradient[i * N + j] = P * tempValue / Separating_Point[i * N + j]  + lower_bound[i][j] * alpha[j] - upper_bound[i][j] * beta[j];         
            }
            
            // cut structure
            if (pbd_mode == MODE_SINGLE) {
                if (num_of_unfixed == M * N) {
                    benders_cut->val[i * N + j] = subgradient[i * N + j];
                    all_subgradient_sum += subgradient[i * N + j] * Separating_Point[i * N + j];
                }
            }
            else if (pbd_mode == MODE_MULTI) {
                if (num_of_unfixed == M * N) {
                    benders_cut->val[j + i * (N + 1)] = subgradient[i * N + j];
                    i_subgradient_sum += subgradient[i * N + j] * Separating_Point[i * N + j];
                }
                else {// adjust for reduced model
                    if (fixed_point[i * N + j] == -1) {
                        i_subgradient_sum += subgradient[i * N + j] * Separating_Point[i * N + j];
                        benders_cut->val[unfixed_idx] = subgradient[i * N + j];
                        unfixed_idx++;
                    }
                }
                theta_try_solution[i] = vopt_i;
                benders_cut->rhs[i] = i_subgradient_sum - vopt_i;
            }            
        }
        unfixed_idx++; // for theta_i coeff only for MODE_MULTI reduced model        
        
        voptsum += vopt_i;
        //printf("Subproblem %d: vopt_i=%.6f\n", i, vopt_i);        
    }
   
    // Get the new objective value to update lp upper bound
    LP_UpperBound += voptsum;

    // cut structure
    if (pbd_mode == MODE_SINGLE) {
        if (num_of_unfixed < M * N) {// adjust for reduced model
            for (int i = 0; i < num_of_unfixed; i++) {
                benders_cut->val[i] = subgradient[unfixed_index[i]];
                all_subgradient_sum += benders_cut->val[i] * Y[i];
            }
        }
        theta_try_solution[0] = voptsum;
        benders_cut->rhs[0] = all_subgradient_sum - voptsum;
    }
}

void MIPSolving(GRBenv* env, GRBmodel* model) {

    printf("\n>>> Start Branching <<<\n");
    int error = 0;

    // Set the variables to binary type
    if (num_of_unfixed == M * N) {
        char* vtype_array = create_char_vector(M * N);
        for (int i = 0; i < M * N; i++) {
            vtype_array[i] = GRB_BINARY;
        }
        error = GRBsetcharattrarray(model, GRB_CHAR_ATTR_VTYPE, 0, M * N, vtype_array);
        if (error) goto QUIT;
        free(vtype_array);
    }
    else {
        char* vtype_array_partial = create_char_vector(num_of_unfixed);
        for (int i = 0; i < num_of_unfixed; i++) {
            vtype_array_partial[i] = GRB_BINARY;
        }
        error = GRBsetcharattrarray(model, GRB_CHAR_ATTR_VTYPE, 0, num_of_unfixed, vtype_array_partial);
        if (error) goto QUIT;
        free(vtype_array_partial);
    }

    //// Turn off Gurobi output
    //error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 1);
    //if (error) goto QUIT;

    // Set integer feasibility tolerance
    error = GRBsetdblparam(GRBgetenv(model), "IntFeasTol", 1e-9);
    if (error) goto QUIT;
    // Enable lazy constraints
    error = GRBsetintparam(GRBgetenv(model), "LazyConstraints", 1);
    if (error) goto QUIT; 
    // set time limit for MIP solving
	double time_limit = 1000 - RootRuntime - HeurRuntime - FixingRuntime;
    if (time_limit > 0.1) {
        error = GRBsetdblparam(GRBgetenv(model), "TimeLimit", time_limit);
        if (error) goto QUIT;
    }else
        goto QUIT;
    

    // callback function to add Benders cuts
    if (pbd_mode == MODE_SINGLE) {
        error = GRBsetcallbackfunc(model, SingleBendersCutCallback, NULL);
        if (error) goto QUIT;
    }
    else if (pbd_mode == MODE_MULTI) {
        error = GRBsetcallbackfunc(model, MultiBendersCutCallback, NULL);
        if (error) goto QUIT;
    }

    printf("Optimizing the MIP model with %d unfixed variables...\n", num_of_unfixed);
    //printf("Current pre_fix_sum: %.6f\n", pre_fix_sum);
    // 设置模型的 ObjConstant 属性
    error = GRBsetdblattr(model, GRB_DBL_ATTR_OBJCON, pre_fix_sum);
	if (error) goto QUIT;

    clock_t start, end;
    double NodesRuntime, nodecount;

    // Optimize the model
    start = clock();
    error = GRBoptimize(model);
    if (error) goto QUIT;
    end = clock();


	NodesRuntime = (double)(end - start) / CLOCKS_PER_SEC;
    // Get the optimization status
    int optimstatus;
    error = GRBgetintattr(model, "Status", &optimstatus);
    if (error) goto QUIT;

    printf("\nOptimization complete\n");
    //check the optimization status
    if (optimstatus == GRB_OPTIMAL) {
        // Get the objective value       
        GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &UpperBound);
        GRBgetdblattr(model, "ObjBoundC", &LowerBound);   // 当前最好下界
        GRBgetdblattr(model, GRB_DBL_ATTR_MIPGAP, &Gap);
        GRBgetdblattr(model, GRB_DBL_ATTR_NODECOUNT, &nodecount);


        printf("MIP_ObjVal: %f\n", UpperBound);
        printf("MIP_Gap: %f\n", Gap);
        printf("MIP_Nodes: %f\n", nodecount);
        printf("MIP_Runtime: %f\n", NodesRuntime);

        fprintf(output, "LowerBound:%.6f;UpperBound:%.6f;Gap:%.6f;BranchTime:%.4f;Nodes:%f;Status:%d;\n",
            LowerBound, UpperBound, Gap, NodesRuntime, nodecount, optimstatus);
    }
    else if (optimstatus == GRB_INF_OR_UNBD) {
        printf("Model is infeasible or unbounded.\n");
    }
    else if (optimstatus == GRB_INFEASIBLE) {
        printf("Model is infeasible.\n");
    }
    else if (optimstatus == GRB_UNBOUNDED) {
        printf("Model is unbounded.\n");
    }
    else if (optimstatus == GRB_TIME_LIMIT) {
        printf("Time limit reached.\n");

        GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &UpperBound); 
        GRBgetdblattr(model, "ObjBoundC", &LowerBound);   // 当前最好下界
        GRBgetdblattr(model, GRB_DBL_ATTR_MIPGAP, &Gap);
        GRBgetdblattr(model, GRB_DBL_ATTR_NODECOUNT, &nodecount);

        printf("MIP_ObjVal: %f\n", UpperBound);
        printf("MIP_Gap: %f\n", Gap);
        printf("MIP_Nodes: %f\n", nodecount);
        printf("MIP_Runtime: %f\n", NodesRuntime);

        fprintf(output, "LowerBound:%.6f;UpperBound:%.6f;Gap:%.6f;BranchTime:%.4f;Nodes:%f;Status:%d;\n",
            LowerBound, UpperBound, Gap, NodesRuntime, nodecount, optimstatus);
    }
    else {
        printf("Optimization ended with status %d\n", optimstatus);
    }

    //check the solution
    if (pbd_mode == MODE_SINGLE) {
        if (num_of_unfixed == M * N) {
            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, M * N + 1, Separating_Point);
            if (error) goto QUIT;
        }
        else {
            double* Xsol = create_double_vector(num_of_unfixed + 1);
            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, num_of_unfixed + 1, Xsol);
            if (error) goto QUIT;
            for (int i = 0; i < num_of_unfixed; i++) {
                Separating_Point[unfixed_index[i]] = Xsol[i];
            }
            Separating_Point[M * N] = Xsol[num_of_unfixed];
            free(Xsol);
        }
    }
    else if (pbd_mode == MODE_MULTI) {
        if (num_of_unfixed == M * N) {
            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, M * N + M, Separating_Point);
            if (error) goto QUIT;
        }
        else {
            double* Xsol = create_double_vector(num_of_unfixed + M);
            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, num_of_unfixed + M, Xsol);
            if (error) goto QUIT;
            for (int i = 0; i < num_of_unfixed; i++) {
                Separating_Point[unfixed_index[i]] = Xsol[i];
            }
            for (int i = 0; i < M; i++)
                Separating_Point[M * N + i] = Xsol[num_of_unfixed + i];
            free(Xsol);
        }
    }

    // 获取其成本    
    double obj_val = 0.0;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
			if (Separating_Point[i * N + j] < 0) Separating_Point[i * N + j] = 0;
			if (Separating_Point[i * N + j] > 1) Separating_Point[i * N + j] = 1;
            obj_val += fixed_charge[i][j] * Separating_Point[i * N + j];
            best_try_solution[i * N + j] = Separating_Point[i * N + j];
        }
    }
    // Solve the subproblem
    CutSeparator();
    if (pbd_mode == MODE_SINGLE) {
        obj_val += theta_try_solution[0];
    }
    else if (pbd_mode == MODE_MULTI) {
        for (int i = 0; i < M; i++)
            obj_val += theta_try_solution[i];
    }
    printf("\nChecking Final Objective Value: %f\n", obj_val);

    // check feasibility
    checkSolutionFeasibility();

QUIT:
    /* Error reporting */
    if (error) {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }
}


/* Define my callback function for reduced model*/
int __stdcall SingleBendersCutCallback(GRBmodel* model, void* cbdata, int where, void* usrdata) {
    int error = 0;
    if (where == GRB_CB_MIPSOL) {
        /* MIP solution callback */
        error = GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, Y);
        for (int i = 0; i < num_of_unfixed; i++) {
            if (Y[i] < 0) Y[i] = 0;
            if (Y[i] > 1) Y[i] = 1;
            Separating_Point[unfixed_index[i]] = Y[i];//map back to the original variable index
        }
        // Solve the subproblem
        CutSeparator();
        // assess whether the cut is violated
        if (Y[num_of_unfixed] + 1e-5 < theta_try_solution[0]) {
        error = GRBcblazy(cbdata, num_of_unfixed + 1, benders_cut->ind, benders_cut->val, benders_cut->sense[0], benders_cut->rhs[0]);
        }
    }
    else if (where == GRB_CB_MIPNODE) {
        int status;
        /* MIP node callback */
        error = GRBcbget(cbdata, where, GRB_CB_MIPNODE_STATUS, &status);
        if (status == GRB_OPTIMAL) {
            error = GRBcbget(cbdata, where, GRB_CB_MIPNODE_REL, Y);
            for (int i = 0; i < num_of_unfixed; i++) {
                if (Y[i] < 0) Y[i] = 0;
                if (Y[i] > 1) Y[i] = 1;
                Separating_Point[unfixed_index[i]] = Y[i];//map back to the original variable index
            }
            // Solve the subproblem
            CutSeparator();
            // assess whether the cut is violated
            if (Y[num_of_unfixed] + 1e-5 < theta_try_solution[0]) {
            error = GRBcbcut(cbdata, num_of_unfixed + 1, benders_cut->ind, benders_cut->val, benders_cut->sense[0], benders_cut->rhs[0]);
            }
        }
    }
    return error;
}

int __stdcall MultiBendersCutCallback(GRBmodel* model, void* cbdata, int where, void* usrdata) {
    int error = 0;
    if (where == GRB_CB_MIPSOL) {
        /* MIP solution callback */
        error = GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, Y);
        for (int i = 0; i < num_of_unfixed; i++) {
            if (Y[i] < 0) Y[i] = 0;
            if (Y[i] > 1) Y[i] = 1;
            Separating_Point[unfixed_index[i]] = Y[i];//map back to the original variable index
        }
        // Solve the subproblem
        CutSeparator();
        int idx = 0;
        for (int i = 0; i < M; i++) {
            // assess whether the cut is violated
            if (Y[num_of_unfixed + i] + 1e-5 < theta_try_solution[i]) {
                error = GRBcblazy(cbdata, num_of_coeff_per_cut[i], &(benders_cut->ind[idx]), &(benders_cut->val[idx]), benders_cut->sense[i], benders_cut->rhs[i]);
            }
            idx += num_of_coeff_per_cut[i];
        }
    }
    else if (where == GRB_CB_MIPNODE) {
        int status;
        /* MIP node callback */
        error = GRBcbget(cbdata, where, GRB_CB_MIPNODE_STATUS, &status);
        if (status == GRB_OPTIMAL) {
            error = GRBcbget(cbdata, where, GRB_CB_MIPNODE_REL, Y);
            for (int i = 0; i < num_of_unfixed; i++) {
                if (Y[i] < 0) Y[i] = 0;
                if (Y[i] > 1) Y[i] = 1;
                Separating_Point[unfixed_index[i]] = Y[i];//map back to the original variable index
            }
            // Solve the subproblem
            CutSeparator();
            int idx = 0;
            for (int i = 0; i < M; i++) {
                // assess whether the cut is violated
                if (Y[num_of_unfixed + i] + 1e-5 < theta_try_solution[i]) {
                error = GRBcbcut(cbdata, num_of_coeff_per_cut[i], &(benders_cut->ind[idx]), &(benders_cut->val[idx]), benders_cut->sense[i], benders_cut->rhs[i]);
                }
                idx += num_of_coeff_per_cut[i];
            }
        }
    }
    return error;
}

void checkSolutionFeasibility(void) {    
    // 1. 检查每个任务是否分配给且仅分配给一台机器
    for (int j = 0; j < N; j++) {
        int assigned = 0;
        for (int i = 0; i < M; i++) {
            if (best_try_solution[i * N + j] >= 1 - 1e-5 && best_try_solution[i * N + j] <= 1 + 1e-5)
                assigned++;
        }
        if (assigned != 1) {
            printf("Task %d assigned to %d machines (should be 1)\n", j, assigned);
        }
    }

    // 2. 检查每台机器容量
    for (int i = 0; i < M; i++) {
        double used = 0.0;
        for (int j = 0; j < N; j++) {
            if (best_try_solution[i * N + j] >= 1 - 1e-5 && best_try_solution[i * N + j] <= 1 + 1e-5)
                used += lower_bound[i][j];
        }
        if (used > capacity[i] + 1e-5) {
            printf("Machine %d overloaded: used %.2f, capacity %.2f\n", i, used, capacity[i]);
        }
    }
    // 3. 检查每台机器的 x 值是否在上下界内
    for (int i = 0; i < M; i++) {
        double used = 0.0;
        for (int j = 0; j < N; j++) {
            if (best_try_solution[i * N + j] == 1 && x[i * N + j] < lower_bound[i][j])
                printf("Task %d on Machine %d exceeds lower bound: x %.2f, lb %.2f\n", j, i, x[i * N + j], lower_bound[i][j]);
            if (x[i * N + j] > upper_bound[i][j])
                printf("Task %d on Machine %d exceeds upper bound: x %.2f, ub %.2f\n", j, i, x[i * N + j], upper_bound[i][j]);
            used += x[i * N + j];
        }
        if (used > capacity[i] + 1e-5) {
            printf("Machine %d overloaded: used %.2f, capacity %.2f\n", i, used, capacity[i]);
        }
    }
}