import os

#rm run-r-1.out
#rm run-mix-1.ener
#rm input-1.psf

os.mkdir('fail')
for dir in os.listdir('.'):
    if 'run-fssh' in dir:
        os.system('rm %s/*restart*' % dir)
        os.system('rm %s/input-1.psf' % dir)
        os.system('rm %s/run-mix-1.ener' % dir)
        os.system('rm %s/run-r-1.out' % dir)
        failed = True
	try:
	        with open('%s/run.log' % dir) as logfile:
        	    for line in logfile.readlines():
                	if 'PROGRAM ENDED AT' in line:
                    		failed = False
	except:
		pass
        if failed:
            os.system('mv %s fail/%s' % (dir, dir))

