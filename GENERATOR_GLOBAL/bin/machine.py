import os, pwd
from subprocess import Popen, PIPE


path_cp2k_dict = {
    'seacourt.phys.ucl.ac.uk' : {
        'acarof' : '/scratch/acarof/bin/cp2k.sopt'
    },
    'scutum.phys.ucl.ac.uk'   : {
        'acarof' : '/scratch/grudorff/antoine/bin/cp2k.sopt'
    },
    'bulstake.phys.ucl.ac.uk'   : {
        'sgiannini' : '/scratch/sgiannini/bin/TRIVIAL_merged_code/cp2k.sopt',
        'acarof' : '/scratch/acarof/bin/cp2k.sopt'
    },
    'thame.phys.ucl.ac.uk': {
        'acarof' : '/scratch/acarof/bin/cp2k.sopt'
    },
    'shuttle.phys.ucl.ac.uk': {
        'acarof' : '/scratch/acarof/bin/cp2k.sopt'
    }
}

source_cp2k_dict = {
    'seacourt.phys.ucl.ac.uk' : {
        'acarof' : '/scratch/acarof/src/CP2K/flavoured-cptk/cp2k/tools/toolchain/install/setup'
    }
}

def get_cp2k_path():
    hostname = os.uname()[1]
    print "HOSTNAME: %s" % hostname
    username = pwd.getpwuid( os.getuid() )[ 0 ]
    print "USER: %s" % username
    if path_cp2k_dict.get(hostname) is not None:
        if path_cp2k_dict[hostname].get(username) is not None:
            if not os.path.isfile(path_cp2k_dict[hostname][username]):
                raise SystemExit('WARNING: check path for CP2K executable in bin/machine.py')
            else:
                return path_cp2k_dict[hostname][username]
        else:
            print "WARNING: please provide a correct username."
            print "Check bin/machine.py: path_cp2k_dict"
            print "List: %s" % path_cp2k_dict[hostname].keys()
            raise SystemExit
    else:
        print "WARNING: please provide a correct hostname."
        print "Check bin/machine.py: path_cp2k_dict"
        print "List: %s" % path_cp2k_dict.keys()
        raise SystemExit




def source(script, update=1):
    pipe = Popen(". %s; env" % script, stdout=PIPE, shell=True)
    data = pipe.communicate()[0]

    env = dict((line.split("=", 1) for line in data.splitlines()))
    if update:
        os.environ.update(env)

    return env


def shell_source(script):
    """Sometime you want to emulate the action of "source" in bash,
    settings some environment variables. Here is a way to do it."""
    import subprocess, os
    pipe = subprocess.Popen(". %s; env" % script, stdout=subprocess.PIPE, shell=True)
    output = pipe.communicate()[0]
    env = dict((line.split("=", 1) for line in output.splitlines()))
    os.environ.update(env)


def source_cp2k():
    hostname = os.uname()[1]
    username = pwd.getpwuid( os.getuid() )[ 0 ]
    if source_cp2k_dict.get(hostname) is not None:
        if source_cp2k_dict[hostname].get(username) is not None:
            cmd = source_cp2k_dict[hostname][username]
            source( cmd )
            #shell_source(cmd)
            #execfile(cmd)
            return cmd
            print "Source command: %s." % cmd
    else:
        print "No source command to do."
