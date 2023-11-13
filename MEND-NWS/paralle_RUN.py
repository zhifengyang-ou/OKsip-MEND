import os
import shutil
import sys
from itertools import combinations
import multiprocessing
import re
import time



def worker_function(lock, process_id):
    filename = 'MEND_namelist.nml'
    folder = ''
    folder_par='userio\\out\\c'+str(process_id)
    filename_par=os.path.join(folder_par, 'OKW_SIM_obs.out')
    with open(filename_par, 'r') as file:
        previous_line=""
        for line in file:
            if 'LF0' in previous_line:
                par_line = line
            previous_line=line
    with lock:
        print(f"Process {process_id} is modifying the file.")
        time.sleep(10)
        # Define the replacement value
        import re
        shutil.copyfile('MEND_namelist_base2.nml', 'MEND_namelist.nml')
        replacement = '1'  # Replace the second '0' with 'X'

        # Create a list to store the modified lines
        modified_lines = []
        with open(filename, 'r') as file:
            for line in file:
                if 'Pinitial =' in line:
                    modified_line = 'Pinitial ='+par_line
                    modified_lines.append(modified_line)
                        
                else:
                    modified_lines.append(line)
                # Modify the lines as needed
            modified_lines[28] = "    Dir_Output  = 'userio/out/cRUN"+str(process_id)+"'\n"
        # Write the modified lines back to the file
        with open(filename, 'w') as file:
            file.writelines(modified_lines)  
            print(f"Process {process_id} finished modifying the file.")
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if 'Dir_Output  = ' in line:
                folder = line.split('=')[1].strip().replace("'", '').replace('/', '\\')

    if folder:
        os.makedirs(folder, exist_ok=True)

    shutil.copy(filename, os.path.join(folder, 'MEND_namelist.nml'))
    os.system(".\\dist\\Debug\\Cygwin-Windows\\mendokw.exe")
    

if __name__ == "__main__":
    
    ## from sensitivity analysis
    #original_list = ['Yg',    'Vg'   , 'KD', 'alpha',   'pEP'   ,'rE','fpEM' ,'WPA2D', 'kYg', 'gamma','Q10','beta','fD','gD']
    original_list = [20,    17   , 19, 18,   13   ,12,14 ,25, 21, 23,22,24,15,16]
    ## from uncertainty analysis
    #original_list=["fD", "alpha", "WPA2D", "KD" ,   "Yg"  ,  "Vg" ,   "Q10" ,  "rE" ,   "pEP" ,  "gD",    "Ygsl" , "fpEM" , "beta",  "gamma"]


    n = 14  # Number of top elements to select
    combinations_list=[]
    for m in range(14):
        m = m+1  # Number of combinations to generate

        # Generate combinations of the top n values
        combinations_list_m = list(combinations(original_list[:n], m))[:3]
        combinations_list.extend(combinations_list_m)
    for j in range(4):
        num_processes = 10

        lock = multiprocessing.Lock()


        processes = []
        for i in range(num_processes):
            process = multiprocessing.Process(target=worker_function, args=(lock, i+j*num_processes))
            processes.append(process)
            process.start()

        for process in processes:
            process.join()


        print("Finished 10 processes..")

