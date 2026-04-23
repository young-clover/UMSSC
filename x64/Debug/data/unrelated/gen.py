import random
import os
import glob

def generate_robust_instance(input_file, output_file):
    """
    读取原始算例文件，生成整数区间参数，输出新的算例文件
    
    参数:
        input_file: 原始算例文件名
        output_file: 输出文件名
    """
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # 移除空行
    lines = [line.rstrip('\n') for line in lines if line.strip()]
    
    # 读取第一行：m n
    m, n = map(int, lines[0].split())
    
    # 初始化数据结构
    costs = []
    resources = []
    
    # 当前读取行索引
    idx = 1
    
    # 读取cost矩阵
    for i in range(m):
        # 读取第i个agent的cost
        row = list(map(float, lines[idx].split()))
        while len(row) < n:
            idx += 1
            row.extend(map(float, lines[idx].split()))
        costs.append(row)
        idx += 1
    
    # 读取resource矩阵
    for i in range(m):
        # 读取第i个agent的resource consumption
        row = list(map(float, lines[idx].split()))
        while len(row) < n:
            idx += 1
            row.extend(map(float, lines[idx].split()))
        resources.append(row)
        idx += 1
    
    # 读取capacity    
    row = list(map(float, lines[idx].split()))
    while len(row) < m:
        idx += 1
        row.extend(map(float, lines[idx].split()))
    capacities = row
    
    
    # 生成整数区间参数
    lower_bounds = []
    upper_bounds = []
    
    for i in range(m):
        lower_row = []
        upper_row = []
        for j in range(n):
            w_ij = resources[i][j]
            # 生成下界: U[0.5, 0.8] * w_ij，然后四舍五入取整
            lower = round(random.uniform(0.5, 0.8) * w_ij)
            # 确保下界至少为1
            lower = max(1, lower)
            
            # 生成上界: U[1.2, 1.5] * w_ij，然后四舍五入取整
            upper = round(random.uniform(1.2, 1.5) * w_ij)
            # 确保上界至少比下界大1
            upper = max(lower + 1, upper)
            
            lower_row.append(int(lower))
            upper_row.append(int(upper))
        lower_bounds.append(lower_row)
        upper_bounds.append(upper_row)
    
    # 确保输出文件夹存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 写入新的算例文件，没有空行
    with open(output_file, 'w') as f:
        # 写入m n
        f.write(f"{m} {n}\n")
        
        # 写入cost矩阵（每个agent一行）
        for i in range(m):
            for j in range(n):
                f.write(f"{int(costs[i][j])} ")
            f.write("\n")
        
        # 写入原始resource矩阵（每个agent一行）
        for i in range(m):
            for j in range(n):
                f.write(f"{int(resources[i][j])} ")
            f.write("\n")
        
        # 写入capacity（每个agent一行）
        for i in range(m):
            f.write(f"{int(capacities[i])}\n")
        
        # 写入下界矩阵（追加到文件尾部，格式与resource矩阵相同）
        for i in range(m):
            for j in range(n):
                f.write(f"{lower_bounds[i][j]} ")
            f.write("\n")
        
        # 写入上界矩阵（追加到文件尾部，格式与resource矩阵相同）
        for i in range(m):
            for j in range(n):
                f.write(f"{upper_bounds[i][j]} ")
            f.write("\n")
    
    return output_file

def process_folder(folder_name):
    """
    处理指定文件夹中的所有算例文件
    
    参数:
        folder_name: 原始文件夹名（如'A', 'B'等）
    """
    # 创建输出文件夹名（添加下划线后缀）
    output_folder = f"{folder_name}_"
    
    # 检查原始文件夹是否存在
    if not os.path.exists(folder_name):
        print(f"警告: 文件夹 '{folder_name}' 不存在，跳过处理")
        return 0
    
    # 创建输出文件夹
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"创建输出文件夹: {output_folder}")
    
    # 获取文件夹中所有文件（假设都是算例文件）
    # 可以根据需要修改文件扩展名过滤
    input_files = glob.glob(os.path.join(folder_name, "*.txt"))
    
    if not input_files:
        # 如果没有.txt文件，尝试获取所有文件
        input_files = [f for f in glob.glob(os.path.join(folder_name, "*")) 
                      if os.path.isfile(f) and not f.endswith('.py')]
    
    if not input_files:
        print(f"警告: 文件夹 '{folder_name}' 中没有找到算例文件")
        return 0
    
    print(f"处理文件夹 '{folder_name}'，找到 {len(input_files)} 个文件")
    
    processed_count = 0
    for input_file in input_files:
        try:
            # 生成输出文件名
            file_name = os.path.basename(input_file)
            # 获取文件名（不含扩展名）并添加_bounds后缀
            base_name = os.path.splitext(file_name)[0]
            output_file = os.path.join(output_folder, f"{base_name}_bounds.txt")
            
            # 处理文件
            generate_robust_instance(input_file, output_file)
            processed_count += 1
            print(f"  已处理: {file_name} -> {os.path.basename(output_file)}")
        except Exception as e:
            print(f"  处理文件 '{input_file}' 时出错: {e}")
    
    print(f"文件夹 '{folder_name}' 处理完成: {processed_count}/{len(input_files)} 个文件\n")
    return processed_count

def batch_process_all_folders():
    """
    批量处理所有文件夹
    """
    print("=== 批量处理算例文件 ===")
    print("说明: 将处理当前目录下A, B, C, D, E文件夹中的所有算例文件")
    print("      生成的文件将保存在对应的A_, B_, C_, D_, E_文件夹中")
    print("=" * 50)
    
    # 设置随机种子以保证可重复性
    random.seed(42)
    
    # 要处理的文件夹列表
    folders = ['A', 'B', 'C', 'D', 'E']
    
    total_processed = 0
    total_files = 0
    
    for folder in folders:
        # 统计文件夹中的文件数
        if os.path.exists(folder):
            files_in_folder = len([f for f in os.listdir(folder) 
                                  if os.path.isfile(os.path.join(folder, f))])
            total_files += files_in_folder
        
        # 处理文件夹
        processed = process_folder(folder)
        total_processed += processed
    
    print("=" * 50)
    print(f"总共处理了 {total_processed} 个文件")
    
    # 显示生成的文件夹结构
    print("\n生成的文件夹结构:")
    for folder in folders:
        output_folder = f"{folder}_"
        if os.path.exists(output_folder):
            file_count = len([f for f in os.listdir(output_folder) 
                            if os.path.isfile(os.path.join(output_folder, f))])
            print(f"  {output_folder}/: {file_count} 个文件")
        else:
            print(f"  {output_folder}/: (文件夹不存在)")



# 主程序
if __name__ == "__main__":
    print("=== 算例文件批量处理工具 ===")
    print("功能: 批量处理A, B, C, D, E文件夹中的算例文件")
    print("      生成带区间参数的算例文件到A_, B_, C_, D_, E_文件夹")
    print("-" * 50)
    
    # 检查当前目录下是否有A, B, C, D, E文件夹
    existing_folders = [f for f in ['A', 'B', 'C', 'D', 'E'] if os.path.exists(f)]
    
    if not existing_folders:        
        print("请确保当前目录下有A, B, C, D, E文件夹，然后重新运行程序")
        exit(0)
    else:
        print(f"检测到以下文件夹: {', '.join(existing_folders)}")
        
    
    # 开始批量处理
    print("\n开始批量处理...")
    batch_process_all_folders()    