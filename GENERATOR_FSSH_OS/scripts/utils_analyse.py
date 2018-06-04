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
    filename = 'Histo-' + property + '-' + '-'.join(tuple) + '.dat'
    with open('%s/%s' % (dest, filename), 'w') as file:
        file.write('#  Bins   Values')
        for bin, value in zip( bins, values):
            file.write('%s   %s\n' % (bin, value))


def print_dict( dict_, property, tuple, dest = '.'):
    filename = property + '-' + '-'.join(tuple) + '.dat'
    file = open('%s/%s' % (dest, filename), 'w')
    for time in sorted(dict_):
        try:
            line = '%f  %s\n' % (time, '   '.join(map(str, dict_[time])) )
        except:
            line = '%f  %s\n' % (time, str(dict_[time]))
        file.write(line)
    file.close()
    return filename

def print_list_dict( list_, property, tuple, dest = '.'):
    filename = 'Block' + '-' + property + '-' + '-'.join(tuple) + '.dat'
    with open('%s/%s' % (dest, filename), 'w') as file:
        for time in sorted(list_[0]):
            result = []
            for index, dict_ in enumerate(list_):
                result = dict_[time]
                line = '%f  %s %s\n' % (time, index, '   '.join(map(str, result)) )
                file.write(line)
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


def initialize_properties_dict(properties_dict, dict_properties, number_blocks, histo_info):
    for method in dict_properties:
        properties_dict[method] = {}
        dict_ = properties_dict[method]
        for property in dict_properties[method]:
            if method == 'Block-runs-average':
                dict_[property] = [{}] * number_blocks
            elif method == 'Runs-average':
                dict_[property] = {}
            elif method in ['Mean', 'Specific', 'Initial', 'Last']:
                dict_[property] = []
                dict_[property + 'info'] = ''
            elif method == 'Histogram':
                nbin = histo_info[property]['nbin']
                min_ = histo_info[property]['min_']
                max_ = histo_info[property]['max_']
                dict_[property] = np.array([0] * nbin)
                dict_[property + 'bins'] = np.linspace(min_, max_, nbin)
            else:
                print "Unknown method: %s" % method
                raise SystemExit

    return properties_dict


def analyse_properties(tuple, run_dict, dict_properties, number_blocks = 5, histo_info = {}, msd_info = 0.0,
                       psf_file = '', coms=[]):
    list_dir = run_dict[tuple]
    properties_dict = {}
    index = -1
    properties_dict['Number runs'] = len(list_dir)
    length_block = int(len(list_dir) / number_blocks)
    print "Number runs", len(list_dir)
    print "Length block", length_block
    properties_dict['Length-block'] = length_block
    properties_dict = initialize_properties_dict(properties_dict, dict_properties, number_blocks, histo_info)
    for directory in list_dir:
        print "Do dir:", directory
        dir = FSSHRun(directory)
        index += 1
        for method in dict_properties:
            dict_ = properties_dict[method]
            for property in dict_properties[method]:
                prop = dir.extract(property, msd_info=msd_info, psf_file=psf_file, coms=coms)
                if method == 'Block-runs-average':
                    block = int(index / length_block)
                    #print block
                    dict_[property][block] = sum_two_dict(dict_.get(property)[block], prop)
                elif method == 'Runs-average':
                    dict_[property] = sum_two_dict(dict_.get(property), prop)
                elif method == 'Mean':
                    list = statistics( prop )
                    dict_[property].append(list)
                    line = '%s      %s\n' % ( directory, '    '.join(map(str, list )))
                    dict_[property + 'info'] += line
                elif method == 'Specific':
                    list = prop
                    line = '%s      %s\n' % (directory, '    '.join(map(str, list)))
                    dict_[property + 'info'] += line
                elif method == 'Initial':
                    #prop = dir.extract(property, init = True)
                    line = '%s      %s\n' % (directory, '    '.join(map(str, prop[min(prop.keys())])))
                    dict_[property + 'info'] += line
                elif method == 'Last':
                    try:
                        line = '%s      %s\n' % (directory, '    '.join(map(str, prop[max(prop.keys())])))
                    except:
                        line = '%s      %s\n' % (directory, str(prop[max(prop.keys())]))
                    dict_[property + 'info'] += line
                elif method == 'Histogram':
                    nbin = histo_info[property]['nbin']
                    min_ = histo_info[property]['min_']
                    max_ = histo_info[property]['max_']
                    for time in prop:
                        dict_[property]+= np.array( bin_list(nbin, min_, max_, prop[time]) )
                else:
                    print "Unknown method: %s" % method
                    raise SystemExit

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
    filename = label + '-' + filename
    with open('%s/%s' % (dest, filename), 'w') as file:
        unit = units.get(property)
        if (label == 'Mean'):
            replace = (unit,) * 4
            header = '# Run  Mean (%s)  Std (%s)  Drift (%s/fs)   QMean (%s) \n' % replace
        elif (label == 'Specific'):
            header = headers.get(property)
        elif (label == 'Initial'):
            header = '# Run   Initial value (%s)\n' % unit
        else:
            header = None
        if header is not None:
            file.write(header)
        file.write(text)
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



