import os, re
import numpy as np
from scipy.stats import pearsonr
import subprocess, threading
import sys, shutil
import re, os, time
import numpy as np
import pyDOE2

# 全局变量定义
drc_all = {}
path_all = {}
drc_number = {}

def submit_job():   ## asynchronous
    os.system('qsub kinetic.script')

def job_qsub(i):
    job_thread = threading.Thread(target=submit_job)
    job_thread.start()
    print("wait")
    start1=True
    while start1:
        time.sleep(5)
        command = "qstat -a | grep curve_CO | awk \'{print $10}\'"
        result = os.popen(command).read()
        print(result.strip())
        if not result:
            print("no")
            continue
        else:
            if "R" in result.strip():
                start1=False
            else:
                continue
            while True:
                command = "qstat -a | grep curve_CO | awk \'{print $10}\'"
                result = os.popen(command).read()
                if not result:
                    print("no")
                    break
                if "R" in result.strip():
                    print(result.strip())
                    time.sleep(60)
                    continue
                else:
                    break

    commands = [
            'id=$(qstat -a | grep curve_CO | awk \'{print $1}\')',
            'qdel $id'
    ]
    bash_script = '\n'.join(commands)
    subprocess.run(['bash', '-c', bash_script])    # 使用Linux命令，删除不行的文件
#    last_line=get_last_line("print-out")
#    print(last_line)
#    if "Run log in the file" not in last_line:
#        shutil.copy("print-out", 'last_line')
#        commands = [
#            'cd "result_1"',
#            'log=$(ls -1 | grep -E \'^log[0-9]+$\' | sort -V | tail -n 1)',
#            'echo "$log"',
#            'log_number=$(echo "$log" | awk -F"log|[^0-9]+" \'{{print $2}}\')',
#            'echo "$log_number"',
#            'rm "$log" "ORG$log_number.mat"'
#        ]
#        bash_script = '\n'.join(commands)
#        subprocess.run(['bash', '-c', bash_script])
    return i

def process_reaction_line(line):
    """处理包含化学反应的文本行，移除反应式中的空格"""
    pattern = re.compile(r'^(\s*.*?[^ ])(\s+\[.*)$') 
    match = pattern.match(line)
    
    if not match:
        print('process_reaction_line not match ',line)
        return line
    
    reaction = match.group(1).replace(' ', '')
    return reaction + match.group(2)

def parse_marked_line(line):
    line = line.strip()
    #pattern = r"^(#*)\s*(.*?)\s*\[\s*([-\d\.]+)\s+([-\d\.]+)\s*\]"
    pattern = r"^(#*)\s*(.*?)\s*\[\s*([^\s\]]+)\s+([^\s\]]+)\s*\]"

    match = re.match(pattern, line)
    if not match:
        print('parse_marked_line not match',line)
        return None
    hashes = match.group(1)
    reactants = match.group(2).strip()
    Ea = float(eval(match.group(3)))
    deltaE = float(eval(match.group(4)))
    mark_count = len(hashes)
    return {'reaction':reactants, 'Ea': Ea, 'deltaE': deltaE, 'mark_count': mark_count}

def replace_value(line):
    pattern = re.compile(
        r'(\[)\s*'                 # 组1: 左方括号
        r'([+-]?\d+\.\d+?)'        # 组2: 第一个数值
        r'(\s+)'                   # 组3: 数值间空格
        r'([+-]?\d+\.\d+?)'        # 组4: 第二个数值
        r'(\s*\])'                 # 组5: 右方括号
    )
    match = pattern.search(line)
    if not match:
        print('replace_value',line)
        return None
    re_value = float(match.group(2))
    return re_value

def get_added_E(filename, lin):
    add_E = None
    try:
        with open(filename, "r", encoding="UTF-8") as file:
            lines = file.readlines()
            for line in lines:
                line = process_reaction_line(line)
                if lin in line:
                    add_E = replace_value(line)
    except FileNotFoundError:
        return None
    return add_E

def drc_all_add(equation, log_n, drc, added_E, number):
    if equation not in drc_all:
        drc_all[equation] = {}
        drc_all[equation]['log_n'] = []
        drc_all[equation]['drc'] = []
        drc_all[equation]['added_E'] = []
        drc_all[equation]['number'] = 0
    drc_all[equation]['log_n'].append(log_n)
    drc_all[equation]['drc'].append(drc)
    drc_all[equation]['added_E'].append(added_E)
    drc_all[equation]['number'] = number

def get_tof_infulence(filename_log, n):
    replace_path_tem = []
    with open(filename_log, "r", encoding="UTF-8") as file:
        lines = file.readlines()
        data = []
        start_marker = "Base on reactants/products :"
        i = 0
        end_marker = "Simplification :"
        reading = False
        for line in lines:
            if start_marker in line:
                i += 1
                if i < 2:
                    continue
                reading = True
                data = []
            elif reading:
                data.append(line)
                if end_marker in line:
                    break
        if data: # 增加判断防止空数据报错
            data.pop(0)
            data.pop(-1)

        for line in data:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    value = abs(float(parts[-2])) # DRC值取绝对值通常更有意义
                    
                    reaction = parts[-1]
                    if reaction not in drc_number:
                        drc_number[reaction] = 0
                    if abs(value) > 0.0001:
                        drc_number[reaction] += 1
                    
                    # 获取 added_E
                    added_e_val = get_added_E(f'../../INCAR_{n}.m', reaction)
                    drc_all_add(reaction, filename_log, value, added_e_val, drc_number[reaction])
                    replace_path_tem.append((value, reaction))
                except ValueError:
                    continue
        replace_path_tem.sort(reverse=True, key=lambda x: abs(x[0]))
    return replace_path_tem


def natural_sort_key(filename):
    match = re.search(r'log(\d+)$', filename)
    if match:
        return int(match.group(1))
    return 0  

def process_pure_log_files(directory):
    log_pattern = re.compile(r'^log\d+$')
    all_files = os.listdir(directory)
    valid_files = [f for f in all_files if log_pattern.match(f)]
    sorted_files = sorted(valid_files, key=natural_sort_key)
    
    print(f"found {len(sorted_files)} log files.")
    print("-" * 30)
    for idx, filename in enumerate(sorted_files, 1):
        filepath = os.path.join(directory, filename)
        print('filepath',filepath)
        try:
            replace_paths = get_tof_infulence(filepath, natural_sort_key(filename))
            path_all[filepath] = replace_paths
        except Exception as e:
            print(f"Error in {filename}: {str(e)}")
            print("-" * 30)
    return path_all

def judge_k():
    cor_list = []
    print("\nStarting Statistical Analysis...")
    
    for eq_name, values in drc_all.items():
        drc_list = values['drc']
        added_E_list = values['added_E']
        
        valid_indices = [i for i, x in enumerate(added_E_list) if x is not None and drc_list[i] is not None]
        clean_drc = [drc_list[i] for i in valid_indices]
        clean_add_E = [added_E_list[i] for i in valid_indices]

        # 1. 计算基础统计量 (Mean & Std)
        if len(clean_drc) > 0:
            drc_mean = np.mean(clean_drc)
            drc_std = np.std(clean_drc)
            score = abs(drc_mean) + drc_std 
        else:
            drc_mean = 0
            drc_std = 0
            score = 0

        # 2. 检查数据有效性以计算相关性
        if len(clean_add_E) < 2 or len(clean_drc) < 2:
            continue
            
        add_E_std = np.std(clean_add_E)
        drc_data_std = np.std(clean_drc)
        
        if add_E_std == 0 or drc_data_std == 0:
            corr, slope, p_value = 0, 0, 1
        else:
            try:
                corr, p_value = pearsonr(clean_add_E, clean_drc)
                slope, intercept = np.polyfit(clean_add_E, clean_drc, deg=1)
            except Exception as e:
                print(f"Calc Error in {eq_name}: {e}")
                corr, slope, p_value = 0, 0, 1

        cor_list.append((eq_name, abs(corr), slope, p_value, values['number'], drc_mean, drc_std, score))
        
    return cor_list

def barrier_candidate(corr_list, can_threshold=0.5, out_file='INCAR.m'):
    """
    """
    num_tomask = int(can_threshold * len(corr_list))
    # 获取需要保留精确能垒的目标反应集合
    target_equations = {item[0].replace(' ', '') for item in corr_list[:num_tomask]}

    try:
        with open('INCAR1.m', 'r', encoding='UTF-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"[!] 找不到文件 INCAR1.m，跳过生成 {out_file}")
        return

    match_pattern1 = "% examples :"
    match_pattern2 = "% #1i for site i,"
    new_file1_lines = []
    new_equation = []
    new_file1_after = []
    matched_first = True
    matched_second = True
    
    for line in lines:
        if match_pattern1 in line:
            matched_first = False
            new_file1_lines.append(line + '\n')
        if match_pattern2 in line:
            matched_second = False
        if not matched_second:
            new_file1_after.append(line)
        if matched_first:
            new_file1_lines.append(line)
        if not matched_first and matched_second:
            new_equation.append(line)

    lines = new_equation
    updated_lines = []

    for line in lines:
        line_ = parse_marked_line(line)
        if line_:
            line_clean = (line_['reaction']).replace(' ', '')
            deltaE = line_['deltaE']

            if line_clean in target_equations:
                Ea = line_['Ea']
            else:
                Ea = float(max(deltaE, 0.0) + 0.75 )
                
            new_line = f"{line_['reaction']:<40} [ {Ea:<10.4f}  {deltaE:<10.4f} ] \n"
            updated_lines.append(new_line)

    new_file1_lines.extend(updated_lines)
    new_file1_lines.append('\n'*2)
    new_file1_lines.extend(new_file1_after)
    
    with open(out_file, 'w', encoding='UTF-8') as f:
        for line in new_file1_lines:
            f.writelines(line)


if __name__ == "__main__":
    #cwd = os.getcwd()
    cwd= '../../result_1'
    path_all = process_pure_log_files(cwd)
    
    corr_list = judge_k()
    print('corr_list',corr_list)
    
    corr_list.sort(key=lambda x: x[4], reverse=True) 
    
    print('\nTop 5 Critical Steps (Sorted by Count):')
    for item in corr_list[:5]:
        print(f"Rxn: {item[0]:<25} | Count: {item[4]:<5} | Score: {item[7]:.4f}")

    with open('corr_list.txt', 'w', encoding='UTF-8') as f:
        f.write(f"{'Reaction Name':<40} {'Count':<10} {'Score':<10} {'Mean_DRC':<10} {'Std_Dev':<10} {'Corr':<10} {'Slope':<10}\n")
        f.write("-" * 110 + "\n")
        for eq_name, corr, slope, p_value, number, mean, std, score in corr_list:
            line = f"{eq_name:<40} {number:<10} {score:<10.4f} {mean:<10.4f} {std:<10.4f} {corr:<10.4f} {slope:<10.4f}\n"
            f.write(line)
            
    drc_number = dict(sorted(drc_number.items(), key=lambda x: x[1], reverse=True))
    with open('drc_number.txt', 'w', encoding='UTF-8') as f:
        for eq_name, number in drc_number.items():
            f.write(f"{eq_name}: {number}\n")

    threshold_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
    
    print(f"\n开始根据比例 {threshold_list} 批量生成 INCAR 文件...")
    for ratio in threshold_list:
        output_filename = f'INCAR_{ratio}.m'
        barrier_candidate(corr_list, can_threshold=ratio, out_file=output_filename)
        shutil.copy(output_filename, 'INCAR.m')
        i=job_qsub(ratio)
        print(f"已生成: {output_filename} (保留了 Count 排名前 {int(ratio*100)}% 反应的精确能垒)")

    print("\nDone! Results saved.")
