# standard modules
import sys, os, time
from multiprocessing import Pool, cpu_count


# custom modudules
from utils_analyse import *


scripts = 'general'
keywords = ['TEMPERATURE']
dict_properties = {
    'Runs-average' :  ['Adiabatic-populations', 'Surface-populations', 'Delta_E', 'Populations', 'MSD', 'IPR'],
    'Block-runs-average' : [ 'MSD'],
    'Specific' : ['FSSH', 'Detailed-FSSH'],
    'Mean' : ['Total-energy', 'Temperature'],
    'Histogram' : ['Off-diagonals', 'Delta_E'],
    'Initial'   : ['Delta_E'],
    'Last' : ['IPR']
}
number_blocks = 5
# FOR HISTO
histo_info = {
    'Off-diagonals' : {
        'nbin' : 100,
        'max_' : 0.01,
        'min_' : -0.01
    },
    'Delta_E' : {
        'nbin': 100,
        'max_': 0.01,
        'min_': -0.01
    }
}

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
    return analyse_properties(tuple, run_dict, dict_properties, number_blocks=number_blocks, histo_info=histo_info)
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
        for method in dict_properties:
            dict_ = properties[method]
            for property in dict_properties[method]:
                if method in ['Mean', 'Specific', 'Initial', 'Last']:
                    create_file(property, title, dict_[property + 'info'], method, dataname, tuple=tuple)
                elif method == 'Runs-average':
                    dict_[property] = average_dict(dict_[property], properties['Number runs'])
                    print_dict(dict_[property], property, tuple, dataname)
                elif method == 'Block-runs-average':
                    new_list = []
                    for element in dict_[property]:
                        new_list.append(average_dict(element, properties['Length-block']))
                    print_list_dict(new_list, property, tuple, dataname)
                elif method == 'Histogram':
                    print_histo(property, dict_[property], dict_[property + 'bins'], dataname, tuple=tuple)
                else:
                    print "Unknown method: %s" % method
                    raise SystemExit



print "End of the analysis at %s" %\
    datetime.strftime(datetime.now(), "%Y %m %d %H:%M:%S ")


