#include "def.h"

int main(int argc, char* argv[]) {
    // argv[0]=exe, argv[1]=file, argv[2]=P, argv[3]=method, [argv[4]=mode]

    if (argc < 4) {
        usage(argv[0]);
        return 1;
    }

    // 定义后缀数组
    char mode_suffix[20] = "";

    // 参数个数判断          
    if (argc < 5) {
        // Case 1: 没有提供 mode，默认为以PBD single 分配内存, 也可用于 SOCP check solution        
        pbd_mode = MODE_SINGLE;
        if (strcmp(argv[3], "PBD") == 0)
            strcpy(mode_suffix, "_single");
    }
    else {
        // Case 2: 用户提供了参数，解析它
        if (strcmp(argv[4], "single") == 0) {
            pbd_mode = MODE_SINGLE;
            strcpy(mode_suffix, "_single");
        }
        else if (strcmp(argv[4], "multi") == 0) {
            pbd_mode = MODE_MULTI;
            strcpy(mode_suffix, "_multi");
        }
        else {
            fprintf(stderr, "Error: Invalid mode '%s' for PBD. Use 'single' or 'multi'.\n", argv[4]);
        }
    }

    /* --------------------------------------------------------
       Build paths
       -------------------------------------------------------- */
    snprintf(input_text, sizeof(input_text), "data/unrelated/%s", argv[1]);
    snprintf(instances_path, sizeof(instances_path), "data/unrelated/");

    char base[100];
    snprintf(base, sizeof(base), "%.*s",
        (int)(strrchr(argv[1], '.') ? strrchr(argv[1], '.') - argv[1] : strlen(argv[1])),
        argv[1]);
    strcat(instances_path, base);
    strcat(instances_path, "/");

    /* --------------------------------------------------------
       Construct output file name
       -------------------------------------------------------- */
       // 格式: base_P_Method_Suffix.txt
    snprintf(output_text, sizeof(output_text),
        "results/unrelated/%s_%s_%s%s.txt",
        base,       // 实例名
        argv[2],    // P (放在中间)
        argv[3],    // Method (放在后面)
        mode_suffix // 后缀
    );

    /* --------------------------------------------------------
       Print info
       -------------------------------------------------------- */
    printf("Input: %s\nOutput: %s\n", input_text, output_text);

    P = atof(argv[2]); // 读取 P
    printf("Func_Arg (P): %s\n", (P == 0 ? "Hybrid" : argv[2]));

    printf("Method: %s\n", argv[3]); // 读取 Method
    if (strcmp(argv[3], "PBD") == 0) {
        printf("PBD Mode: %s\n", (pbd_mode == MODE_SINGLE ? "Single" : "Multi"));
    }

    /* --------------------------------------------------------
       Read number of instances
       -------------------------------------------------------- */
    input = open_file(input_text, "r");
    int num_inst;
    if (fscanf(input, "%d", &num_inst) != 1) {
        fprintf(stderr, "Error reading number of instances.\n");
        fclose(input);
        exit(1);
    }

    // Create output file
    output = open_file(output_text, "w+");
    fclose(output);

    /* --------------------------------------------------------
       Solve each instance
       -------------------------------------------------------- */
    //for (int i = 0; i < 1; i++) {
    for (int i = 0; i < num_inst; i++) {
        if (fscanf(input, "%49s", instance_name) != 1) {
            break;
        }

        printf("\n++++++++++++++++++++ Solving %s ++++++++++++++++++++\n", instance_name);
        read_instance(instance_name, "unrelated");

        output = open_file(output_text, "a+");
        fprintf(output, "Instance: %s;\n", instance_name);

        if (strcmp(argv[3], "PBD") == 0) {
            PBD();
        }
        else if (strcmp(argv[3], "SOCP") == 0) {
            SOCP();
        }

        fclose(output);
        printf("Finished solving %s\n", instance_name);
        free_memory();
        Sleep(5000);
    }

    fclose(input);
    return 0;
}