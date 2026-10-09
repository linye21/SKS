import subprocess, threading
import sys, shutil
import re, os, time
import numpy as np
import pyDOE2


def submit_job():   ## asynchronous
    os.system('qsub kinetic.script')

def get_last_line(file_path):           # get the penultimate nonempty line
    with open(file_path, 'r') as file:
        lines = file.readlines()
        non_empty_lines_count = 0
        for line in reversed(lines):
            stripped_line = line.strip()
            if stripped_line:
                non_empty_lines_count += 1
                if non_empty_lines_count == 1:
                    return stripped_line
            else:
                continue

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
    subprocess.run(['bash', '-c', bash_script])
    last_line=get_last_line("print-out")
    print(last_line)
    if "Run log in the file" not in last_line:
        shutil.copy("print-out", 'last_line')
        commands = [
            'cd "result_1"',
            'log=$(ls -1 | grep -E \'^log[0-9]+$\' | sort -V | tail -n 1)',
            'echo "$log"',
            'log_number=$(echo "$log" | awk -F"log|[^0-9]+" \'{{print $2}}\')',
            'echo "$log_number"',
            'rm "$log" "ORG$log_number.mat"'
        ]
        bash_script = '\n'.join(commands)
        subprocess.run(['bash', '-c', bash_script])
    return last_line,i

def generate_lhs_samples(n_reactions, n_samples=100):
    lhs = pyDOE2.lhs(
        n_reactions,
        samples=n_samples,
        criterion='maximin',
        random_state=43
    )
    return lhs * 1.5


def replace_value(line, sample_values, reaction_idx):
    prefix, re_value, on_value = process_line(line)

    barrier = sample_values[reaction_idx]
    if prefix is None:
        return line
    if on_value>0:
        re_value=barrier + on_value
    else:
        re_value=barrier + float(0.0)

    new_line=(
        prefix + '[' +
        str(re_value) +          
        '   ' + str(on_value) +
        ']'
    )
    print('new_line:', new_line)
    return new_line

def process_line(line):
    if '[' not in line or ']' not in line:
        return None
    prefix,rest= line.split('[', 1)
    # print(prefix)
    inside,_= rest.split(']', 1)
    parts=re.findall(r'[\d\.\+\-\*/\(\)]+', inside)
    if len(parts) != 2:
        return None
    re_value= float(eval(parts[0]))
    on_value= float(eval(parts[1]))
    return prefix, re_value, on_value

def line_get_replace(file1_lines, sample_values):
    line_alls = []
    reaction_idx = 0 

    for line in file1_lines:
        match=True
        if '[' not in line or ']' not in line:
            match = None
        if match:
            print('line:', line)
            modified_line = replace_value(line, sample_values, reaction_idx)
            line_alls.append(modified_line+ "\n")
            reaction_idx += 1
        else:
            continue
    return line_alls


def main(out_file_name, sample_id):
    with open('INCAR1.m', 'r', encoding="UTF-8") as f:
        file1_lines = f.readlines()
    match_pattern1 = "% examples :"
    match_pattern2 = "Q0 = [1 1 1];"
    new_file1_lines = []
    new_equation_lines = []
    new_file1_after=[]
    matched_first_line = True
    matched_second_line = True
    # file2_content = line_get_replace('INCAR1.m')
    # print('file2_content:',file2_content)
    for line in file1_lines:
        if match_pattern1 in line:
            matched_first_line = False
            new_file1_lines.append(line+"\n")
            # new_file1_lines.extend(file2_content)
        if match_pattern2 in line:
            matched_second_line = False
        if not matched_second_line:
            new_file1_after.append(line)

        elif not matched_first_line and matched_second_line:
            #new_file1_lines.append("\n"+file2_content+"\n")
            new_equation_lines.append(line)
        elif not matched_first_line:
            new_file1_lines.append(line)

    sample_values = LHS_MATRIX[sample_id]
    modified_lines = line_get_replace(new_equation_lines, sample_values)
    modified_lines.append('\n'*2)
    new_file1_lines.extend(modified_lines)
    new_file1_lines.extend(new_file1_after)
    with open(out_file_name, 'w', encoding="UTF-8") as f:
        f.writelines(new_file1_lines)
    shutil.copy(out_file_name, 'INCAR.m')
    
N_REACTIONS =81
print(f"Total reactions to modify: {N_REACTIONS}")
LHS_MATRIX = generate_lhs_samples(N_REACTIONS)  
print(f"LHS samples generated: {LHS_MATRIX.shape}")

if __name__ == "__main__":
    for sample_id in range(100):
        out_file = f'INCAR_{sample_id+1}.m'
        main(out_file, sample_id)
        last_line,i=job_qsub(sample_id+1)
        print('last_line:',last_line,'log_n',i)

