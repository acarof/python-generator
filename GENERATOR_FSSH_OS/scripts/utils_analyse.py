# standard modules

# custom modudules
from utils_scripts import *
from datetime import datetime



def sum_two_dict( dict1, dict2):
    result = {}
    if dict1 is None or dict1 == {}:
        result = dict2
    else:
        for key in dict1:
            value = np.array(dict1[key]) + np.array(dict2[key])
            result[key] = value
    return result


def average_dict(dict1, number):
    result = {}
    for key in dict1:
        result[key] = np.array( dict1[key] ) / number
    return result

def print_list(list_, property, tuple, dest = '.'):
    filename = property + '-' + '-'.join(tuple) + '.dat'
    with open('%s/%s' % (dest, filename), 'w') as file:
        for element in list_:
            line = '%s\n' % element
            file.write(line)


def print_histo(property, values, bins, dest = '.', tuple = ''):
    filename = property + '-' + '-'.join(tuple) + '.dat'
    with open('%s/%s' % (dest, filename), 'w') as file:
        file.write('#  Bins   Values')
        for bin, value in zip( bins, values):
            file.write('%s   %s\n' % (bin, value))


def print_dict( dict_, property, tuple, dest = '.'):
    filename = property + '-' + '-'.join(tuple) + '.dat'
    file = open('%s/%s' % (dest, filename), 'w')
    for time in sorted(dict_):
        line = '%f  %s\n' % (time, '   '.join(map(str, dict_[time])) )
        file.write(line)
    file.close()
    return filename

def print_list_dict( list_, property, tuple, dest = '.'):
    filename = property + '-' + '-'.join(tuple) + '.dat'
    file = open('%s/%s' % (dest, filename), 'w')
    for time in sorted(list_[0]):
        result = []
        for dict_ in list_:
            result.append(dict_[time])
        line = '%f  %s\n' % (time, '   '.join(map(str, result)) )
        file.write(line)
    file.close()
    return filename


def append_two_dict(dict1, dict2):
    result = {}
    if dict1 is None:
        result = dict2
    else:
        for key in dict1:
            value = dict1[key] + dict2[key]
            result[key] = value
    return result


def analyse_properties(tuple, run_dict, dict_properties, number_blocks = 5, nbin = 10, min_ = 0, max_=0):
    real_start_time = datetime.now()
    start_time = datetime.strftime( real_start_time, "%Y %m %d %H:%M:%S ")
    print "One worker for: %s starts at %s" % (tuple, start_time)
    list_dir = run_dict[tuple]
    properties_dict = {}
    index = -1
    properties_dict['Number runs'] = len(list_dir)
    length_block = int(len(list_dir) / number_blocks)
    print "Number runs", len(list_dir)
    print "Length block", length_block
    properties_dict['Length-block'] = length_block
    for directory in list_dir:
        dir = FSSHRun(directory)
        index += 1
        for property in dict_properties.get('Block-runs-average', []):
            prop = dir.extract(property)
            block = int(index/length_block)
            if properties_dict.get(property) is None:
                properties_dict[property] = [{}] * number_blocks
            properties_dict[property][block] = sum_two_dict(properties_dict.get(property)[block], prop)
        for property in dict_properties.get('Time-runs-average', []):
            prop = dir.extract(property)
            #print prop.keys()
            #print list_time
            new_prop = dict( (element, prop[element]) for element in list_time if element in prop)
            #print new_prop
            properties_dict[property] = sum_two_dict(properties_dict.get(property), new_prop)
            properties_dict[property + 'for-index-%s' % index] = properties_dict[property]
        for property in dict_properties.get('Runs-average', []):
            prop = dir.extract(property)
            properties_dict[property] = sum_two_dict(properties_dict.get(property), prop)
        for property in dict_properties.get('Mean', []):
            prop = dir.extract(property)
            list = statistics( prop )
            if properties_dict.get(property) is None:
                properties_dict[property] = []
                properties_dict[property + 'info'] = ''
            properties_dict[property].append(list)
            line = '%s      %s\n' % ( directory, '    '.join(map(str, list )))
            properties_dict[property + 'info'] += line
        for property in dict_properties.get('Specific', []):
            prop = dir.extract(property)
            list = prop
            if properties_dict.get(property) is None:
                properties_dict[property] = []
                properties_dict[property + 'info'] = ''
            line = '%s      %s\n' % (directory, '    '.join(map(str, list)))
            properties_dict[property + 'info'] += line
        for property in dict_properties.get('Initial', []):
            prop = dir.extract(property, init = True)
            line = '%s      %s\n' % (directory, '    '.join(map(str, prop.values()[0])))
            if properties_dict.get(property) is None:
                properties_dict[property] = []
                properties_dict[property + 'info'] = ''
            properties_dict[property + 'info'] += line
        for property in dict_properties.get('Histogram', []):
            prop = dir.extract(property)
            if properties_dict.get(property) is None:
                properties_dict[property] = np.array( [0] * nbin )
                properties_dict[property + 'bins'] = np.linspace(min_, max_, nbin)
            for time in prop:
                properties_dict[property]+= np.array( bin_list(nbin, min_, max_, prop[time]) )

    real_end_time = datetime.now()
    end_time = datetime.strftime( real_end_time, "%Y %m %d %H:%M:%S ")
    print "One worker for: %s ends at %s" % (tuple, end_time)
    time_diff = real_end_time - real_start_time
    hours, remainder = divmod(time_diff.total_seconds(), 3600)
    minutes, seconds = divmod(remainder, 60)
    print "The worker %s lasted: %s hours %s minutes %s seconds" % (tuple, hours, minutes, seconds)
    return properties_dict



def bin_list(nbin, min_, max_, list_):
    histo = [0] * nbin
    step = (max_ - min_) / nbin
    for element in list_:
        length = element - min_
        ind = int(length/step)
        if 0 < ind and ind < nbin:
            histo[ind-1] += 1
        else:
            print "One element not is the list:", element
    return histo




def print_run_dict(run_dict, keywords):
    filename = 'List-run.dat'
    file = open(filename, 'w')
    for tuple in run_dict:
        line = '%s %s\n' % ( tuple, '  '.join(run_dict[tuple]) )
        file.write(line)
    file.close()
    return filename

def create_file(property, title, text,  label='None', dest = '.', tuple = 'None'):
    if tuple is not 'None':
        filename = '%s-%s.dat' % (property, '-'.join(tuple))
    else:
        filename = '%s-%s.dat' % (property, title)
    file = open('%s/%s' % (dest, filename), 'w')
    unit = units.get(property)
    if (label == 'Mean'):
        replace = (unit,) * 4
        header = '# Run  Mean (%s)  Std (%s)  Drift (%s/fs)   QMean (%s) \n' % replace
    elif (label == 'Spec'):
        header = headers.get(property)
    elif (label == 'Initial'):
        header = '# Run   Initial value (%s)\n' % unit
    else:
        header = None
    if header is not None:
        file.write(header)
    file.write(text)
    file.close()
    return  filename


def statistics(prop):
    times = []
    props = []
    for time in sorted(prop):
        times.append(time)
        props.append(prop.get(time))
    mean = np.mean(props)
    std = np.std(props)
    drift = np.polyfit(times, props, 1)[0][0]
    qmean = np.sqrt(np.mean(np.square(props)))
    return [mean, std, drift, qmean]



