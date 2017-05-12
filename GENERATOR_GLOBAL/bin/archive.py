import tarfile
import shutil
import os
import sys


try:
    nworker = int( sys.argv[1] )
except:
    nworker = nworker = -1



def archive(directory):
    if 'run' in directory:
        numero = directory.split('-')[-1]
        tmp = 'tmp-' + numero
        if os.path.exists(tmp):
            os.system('mv %s %s' % (tmp, directory))
        print 'Creating archive of %s' % directory
        tar = tarfile.open(directory + '.tar.gz', 'w:gz')
        tar.add(directory)
        shutil.rmtree(directory)



if nworker == 0:
    for directory in os.listdir('.'):
        archive(directory)
else:
    from multiprocessing import Pool, cpu_count
    if nworker == -1:
        nworker = cpu_count()
    pool = Pool(nworker)
    pool.map(archive, os.listdir('.'))



name_bucket = os.getcwd().split('/')[-1]
if name_bucket != 'GENERATOR_GLOBAL':
    new_name =  name_bucket + '-ARCHIVED'
    os.chdir('..')
    os.system('mv %s %s' % (name_bucket, new_name))




