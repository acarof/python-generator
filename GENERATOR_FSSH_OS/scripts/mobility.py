# standard modules
import sys, os, time
from multiprocessing import Pool, cpu_count


# custom modudules
from utils_analyse import *

scripts = 'mobility'
keywords = ['TEMPERATURE']
dict_properties = {
    'Block-runs-average' : ['MSD']
}
number_blocks = 1
# FOR HISTO
nbin=100
max_ = 0.004
min_ = -0.004

#FIND THE TITLE
name_bucket = os.getcwd().split('/')[-1]
short_time = time.strftime("%y%m%d%H%M", time.localtime())
title = '%s-%s' % (name_bucket, short_time,)
dirlist = os.listdir('.')
dataname = 'data-%s-%s' % (scripts, title)
if not os.path.isdir(dataname):
    os.mkdir(dataname)



# CREATE LIST OF RUNS
run_dict = {}
os.system('cd ..')
for i, directory in enumerate(dirlist):
    if 'run-fssh' in directory and 'per' not in directory:
        keys = tuple(keywords)
        if run_dict.get(keys ) is None:
            run_dict[keys] = []
        run_dict[keys].append(directory)





# PARALLEL OR SERIAL CALCULATION
def super_analyse(tuple):
    return analyse_properties(tuple, run_dict, dict_properties, number_blocks=number_blocks, nbin=nbin, max_=max_, min_=min_)
try:
    nworker = int( sys.argv[1] )
except:
    nworker = -1
if nworker == -1:
    nworker = cpu_count()
if nworker == 0:
    results = []
    for keys in run_dict.keys():
        results.append( super_analyse(keys) )
else:
    pool = Pool(nworker)
    results = pool.map(super_analyse, run_dict.keys())



# PRINT THE RESULTS
results_dict = {}
with open('%s/List-tuple.dat' % dataname, 'w') as filetuple:
    filetuple.write('%s\n' % ('  '.join(keywords)))
    for tuple, result in zip( run_dict.keys(), results):
        filetuple.write('%s\n' % '  '.join(tuple))
        results_dict[tuple] = result
        properties = results_dict[tuple]
        for property in dict_properties.get('Mean', []):
            create_file(property, title, properties[property + 'info'], 'Mean', dataname, tuple = tuple)
        for property in dict_properties.get('Runs-average', []):
            properties[property] = average_dict(properties[property], properties['Number runs'])
            print_dict( properties[property], property, tuple, dataname  )
        for property in dict_properties.get('Block-runs-average', []):
            new_list = []
            for element in properties[property]:
                new_list.append(average_dict(element, properties['Length-block']))
            print_list_dict(new_list, property, tuple, dataname)
        for property in dict_properties.get('Specific', []):
            create_file(property, title, properties[property + 'info'], 'Spec', dataname, tuple = tuple)
        for property in dict_properties.get('Initial', []):
            create_file(property, title, properties[property + 'info'], 'Initial', dataname, tuple = tuple)
        for property in dict_properties.get('Histogram', []):
            print_histo(property, properties[property], properties[property + 'bins'], dataname, tuple = tuple)
            #print properties[property]
            #print properties[property + 'bins']





print "End of the analysis at %s" %\
    datetime.strftime(datetime.now(), "%Y %m %d %H:%M:%S ")


