import os, shutil,sys, re


try:
    task = sys.argv[1]
except:
    print "A TASK TO DUPLICATE IS REQUIRED!"
    raise SystemExit


if 'preparation/' in task:
    taskname = task.split('/')[1]
else:
    taskname = task

try:
    name = sys.argv[2]
except:
    name = 'TONAME'


if not os.path.exists(name):
    os.mkdir(name)
else:
    print ("%s ALREADY EXISTS" % name)
    raise SystemExit


dst=name

srclist = ['bin', 'templates', 'scripts']
for src in srclist:
    shutil.copytree(src, dst + '/' + src)

shutil.copytree(task, name + '/' + taskname)

with open('%s/task.py' % task) as file:
    text = file.read()
    try:
        inits = re.findall(r"'FILE_INIT': '([^']*)", text)[0]
        fileinit = True
        print "Hey"
    except:
        fileinit = False
if fileinit:
   os.mkdir('%s/initial' % name)
   os.system('cp -r initial/*-%s %s/initial/' % (inits, name))