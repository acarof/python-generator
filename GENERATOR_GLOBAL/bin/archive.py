import tarfile
import shutil
import os


for directory in os.listdir('.'):
    if 'run' in directory:
        print 'Creating archive of %s' % directory
        tar = tarfile.open(directory + '.tar.gz', 'w:gz')
        tar.add(directory)
        shutil.rmtree(directory)
