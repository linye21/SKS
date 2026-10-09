import subprocess, threading
import sys, shutil
import re, os, time
import numpy as np
import pyDOE2


def generate_lhs_samples(n_reactions, n_samples=50):
    lhs = pyDOE2.lhs(
        n_reactions,
        samples=n_samples,
        criterion='maximin',
        random_state=543
    )
    return lhs * 1.5


N_REACTIONS =94
print(f"Total reactions to modify: {N_REACTIONS}")
LHS_MATRIX = generate_lhs_samples(N_REACTIONS)  
print(f"LHS samples generated: {LHS_MATRIX.shape}")

def replace_value(line, sample_values, reaction_idx):
    pattern = re.compile(
        r'(\[)\s*'                  # 组1: 左方括号
        r'([+-]?\d+\.\d+?)'         # 组2: 能垒值（将被替换）
        r'(\s+)'                    # 组3: 空格
        r'([+-]?\d+\.\d+?)'         # 组4: 基准值（保留）
        r'(\s*\])'                  # 组5: 右方括号
    )
    match = pattern.search(line)
    if not match:
        return line
    barrier = sample_values[reaction_idx]
    re_value=float(match.group(2))
    on_value=float(match.group(4))

    if on_value>0:
        re_value=barrier + on_value
    else:
        re_value=barrier + float(0.0)

    new_line=(
        match.group(1) + 
        str(re_value) +          
        match.group(3) + 
        match.group(4) + 
        match.group(5)
    )
    new_line=line[:match.start()] + str(new_line) + line[match.end():]
    return new_line

def line_get_replace(file1_lines, sample_values):
    line_alls = []
    reaction_idx = 0  # 当前反应序号
    pattern = re.compile(
        r'(\[)\s*'                  # 组1: 左方括号
        r'([+-]?\d+\.\d+?)'         # 组2: 能垒值（将被替换）
        r'(\s+)'                    # 组3: 空格
        r'([+-]?\d+\.\d+?)'         # 组4: 基准值（保留）
        r'(\s*\])'                  # 组5: 右方括号
    )
    for line in file1_lines:
        match = pattern.search(line)
        if match:
            modified_line = replace_value(line, sample_values, reaction_idx)
            line_alls.append(modified_line)
            reaction_idx += 1
        else:
            continue
    return line_alls


def main(out_file_name, sample_id):
    with open('INCAR1.m', 'r', encoding="UTF-8") as f:
        file1_lines = f.readlines()
    match_pattern1 = "% examples :"
    match_pattern2 = "% #1i for site i, #11 for site 1"
    new_file1_lines = []
    new_equation_lines = []
    new_file1_after=[]
    matched_first_line = False
    matched_second_line = True
    # file2_content = line_get_replace('INCAR1.m')
    # print('file2_content:',file2_content)
    for line in file1_lines:
        if match_pattern1 in line:
            matched_first_line = True
            new_file1_lines.append(line+"\n")
            # new_file1_lines.extend(file2_content)
        if match_pattern2 in line:
            matched_second_line = False
            # matched_second_line = False  # 重置第二个标志
        if not matched_second_line:
            new_file1_after.append(line)

        elif matched_first_line and matched_second_line:
            #new_file1_lines.append("\n"+file2_content+"\n")
            new_equation_lines.append(line)
        elif not matched_first_line:
            new_file1_lines.append(line)

    sample_values = LHS_MATRIX[sample_id]
    modified_lines = line_get_replace(new_equation_lines, sample_values)
    new_file1_lines.extend(modified_lines)
    new_file1_lines.extend(new_file1_after)
    with open(out_file_name, 'w', encoding="UTF-8") as f:
        f.writelines(new_file1_lines)
    shutil.copy(out_file_name, 'INCAR.m')

# ---------- 主流程 ----------
if __name__ == "__main__":
    for sample_id in range(50):
        out_file = f'INCAR_{sample_id+1}.m'
        main(out_file, sample_id)
