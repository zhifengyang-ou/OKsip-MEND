import os
import shutil
import sys
from datetime import datetime

now_string = str(datetime.now())
comment = sys.argv[1]

filename = 'MEND_namelist.nml'
folder = ''

with open(filename, 'r') as file:
    for line in file:
        line = line.strip()
        if 'Dir_Output  = ' in line:
            folder = line.split('=')[1].strip().replace("'", '').replace('/', '\\')

if folder:
    os.makedirs(folder, exist_ok=True)

logfile = 'log.txt'
with open(logfile, 'r') as log_file, open(os.path.join(folder, logfile), 'w') as output_file:
    for line in log_file:
        output_file.write(line)

with open(os.path.join(folder, logfile), 'a') as output_file:
    output_file.write(folder + '\n')
    output_file.write(now_string + '\n')

shutil.copy(os.path.join(folder, logfile), 'log.txt')
with open(os.path.join(folder, 'comment.txt'), 'w') as comment_file:
    comment_file.write(comment)

shutil.copy(filename, os.path.join(folder, 'MEND_namelist.nml'))
os.system(".\\dist\\Debug\\Cygwin-Windows\\mendokw.exe")

