import tarfile
import shutil
import os, sys, re
import subprocess
from datetime import datetime

def du(path):
    """disk usage in human readable format (e.g. '2,1GB')"""
    return subprocess.check_output(['du','-sh', path]).split()[0].decode('utf-8')

success = 0
failed = 0
number_dir = len(os.listdir('.'))
first = number_dir
last = -1
for directory in os.listdir('.'):
    if 'run-' in directory:
        number = int(directory[4:])
        first = min( number, first )
        last  = max( number, last)

        success += 1
    elif 'fail-' in directory:
        failed += 1

file = open('run-%d/run.log' % first)
for line in file.readlines():
    if re.search(r" CP2K\| source code revision number:.*git\:(.{7})", line):
        git_version = re.search(r" CP2K\| source code revision number:.*git\:(.{7})", line).group(1)
    elif re.search(r".*PROGRAM STARTED AT *(.*)", line):
        start = re.search(r".*PROGRAM STARTED AT *(.*)", line).group(1)

file = open('run-%d/run.log' % last)
for line in file.readlines():
    if re.search(r".*PROGRAM ENDED AT *(.*)", line):
        end = re.search(r".*PROGRAM ENDED AT *(.*)", line).group(1)

date_start = datetime.strptime(start, "%Y-%m-%d %H:%M:%S.%f")
date_end = datetime.strptime(end, "%Y-%m-%d %H:%M:%S.%f")
time_diff = date_end - date_start
hours, remainder = divmod(time_diff.total_seconds(), 3600)
minutes, seconds = divmod(remainder, 60)

name_bucket = os.getcwd().split('/')[-1]

final = "\n\n"
final +=  "INFORMATION FOR BUCKET: %s\n\n" % name_bucket

final += "CP2K version: %s\n" % git_version
final += "Task started at: %s\n" % start
final += "Task ended at: %s\n" % end
final += "The task lasted: %s hours %s minutes %s seconds\n\n" % (hours, minutes, seconds)

final += "Size of the bucket: %s\n\n" % du('.')

final += "Number of sucessful runs: %d\n" % success
final += "Number of failed runs: %d\n" % failed



print final
file = open(name_bucket + '.info', 'w')
file.write(final)
file.close()




