# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 15:25:09 2025

@author: Young
"""

import os

# 设置要统计的文件夹路径
folder_path = './E_/'  # 替换为你要操作的文件夹路径
output_file = 'E_.txt'  # 输出文件名

# 获取文件夹中的所有文件（不包含子目录）
file_list = [f for f in os.listdir(folder_path)
             if os.path.isfile(os.path.join(folder_path, f))]

# 按字典序排序
file_list.sort()

# 写入结果到 txt 文件
with open(output_file, 'w') as f:
    f.write(f"{len(file_list)}\n")
    for filename in file_list:
        f.write(f"{filename}\n")

print(f"Done. {len(file_list)} files written to {output_file}.")
