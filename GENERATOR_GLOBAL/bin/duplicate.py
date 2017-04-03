import os, shutil,sys


try:
    task = sys.argv[1]
except:
    print "A TASK TO DUPLICATE IS REQUIRED!"
    raise SystemExit
taskname = task.split('/')[1]


try:
    name = sys.argv[2]
except:
    name = 'TONAME'

if not os.path.exists(name):
    os.mkdir(name)
else:
    print "TONAME ALREADY EXISTS"
    raise SystemExit

dst=name
srclist = ['bin', 'templates', 'scripts','initial']
for src in srclist:
    shutil.copytree(src, dst + '/' + src)

shutil.copytree(task, name + '/' + taskname)