#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

from __future__ import print_function
import sys
import os
import re
import glob

#
# Note that an issue was encountered processing (with Python3) a Fortran90 source file
# that was not UTF-8 encoded.  Because the purpose of this program in to scan C and
# Fortran source files to generate prerequisites/dependencies and because for these
# languages keywords are ASCII, this program deals with encoding issues using the
# "errors=replace" option of the open function.  Any encoding errors will occur within
# comments and can be ignored.
#
# Note that to support using both Python2 and Python3, the open function has been
# aliased.  Use the input_open function to open input files; use the output_open
# function to open output files.
#
if sys.version_info.major == 2:
    import io
    input_open = io.open
    output_open = open
else:
    input_open = open
    output_open = open


def print_err(*pargs, **kwargs):
    """
    Print error message.

    pargs is the list of positional arguments.

    kwargs is the dictionary of key-word arguments.

    pargs contains the message to print.
    """
    print('[ERR]', *pargs, **kwargs)


def print_dbg(*pargs, **kwargs):
    """
    Print debug message.

    pargs is the list of positional arguments.

    kwargs is the dictionary of key-word arguments.

    pargs contains the message to print.

    This routine first checks the verbose setting before printing the
    message.  Note that cli_args is a global variable.
    """
    if cli_args.VERBOSE >= cli_args.DBG:
        print('[DBG]', *pargs, **kwargs)


def print_warn(*pargs, **kwargs):
    """
    Print warning message.

    pargs is the list of positional arguments.

    kwargs is the dictionary of key-word arguments.

    pargs contains the message to print.

    This routine first checks the verbose setting before printing the
    message.  Note that cli_args is a global variable.
    """
    if cli_args.VERBOSE >= cli_args.WARN:
        print('[WARN]', *pargs, **kwargs)


def print_info(*pargs, **kwargs):
    """
    Print info message.

    pargs is the list of positional arguments.

    kwargs is the dictionary of key-word arguments.

    pargs contains the message to print.

    This routine first checks the verbose setting before printing the
    message.  Note that cli_args is a global variable.
    """
    if cli_args.VERBOSE >= cli_args.INFO:
        print('[INFO]', *pargs, **kwargs)


def print_status(*pargs, **kwargs):
    """
    Print status message.

    pargs is the list of positional arguments.

    kwargs is the dictionary of key-word arguments.

    pargs contains the message to print.

    This routine first checks the verbose setting before printing the
    message.  Note that cli_args is a global variable.
    """
    if cli_args.VERBOSE >= cli_args.STATUS:
        print(*pargs, **kwargs)


def get_cli_exclude_files(stype):
    """
    Return either the set of files to exclude from searching.

    stype specifies which set of exclude files to return.

    Note that cli_args is a global variable.
    """
    if stype == 'fmod':
        return cli_args.exclude_fmod_files
    elif stype == 'finc':
        return cli_args.exclude_finc_files
    else:
        raise Exception('[ERR] Unknown exclude files type {}'.format(stype))


def get_cli_search_dirs():
    """
    Return the list of search directories.

    Note that cli_args is a global variable.
    """
    return cli_args.dirs


def get_suffix(name):
    """
    Return the suffix of name.

    Example: given name = 'file.f90', return '.f90'.
    """
    return os.path.splitext(name)[1]


def basename_without_suffix(name):
    """
    Return name without its suffix.

    Example: given name = 'file.f90', return 'file'.
    """
    return os.path.splitext(os.path.basename(name))[0]


def basename_new_suffix(name, suffix):
    """
    Return name with a new suffix.

    Example: given name = 'file.f90' and suffix = '.X', return 'file.X'.

    Note: suffix must contain the dot.
          suffix should be given as '.X' not just 'X'.
    """
    return os.path.splitext(os.path.basename(name))[0]+suffix


def suffix_check(f):
    """
    Return the type of file f by examining its suffix.
    """
    fortran90 = {'.f90', '.F90'}
    fortran77 = {'.f', '.F'}
    c = {'.c', '.cc'}
    suffix = get_suffix(f)
    if suffix in fortran90:
        return 'fortran90'
    elif suffix in fortran77:
        return 'fortran77'
    elif suffix in c:
        return 'c'
    else:
        raise Exception('[ERR] Unknown suffix {}; skipping.'.format(suffix))


def find_file(fname):
    """
    Return the full file name (path plus name) to file fname.

    fname is the name of the file to find.

    If the file fname is not found, then simply return None.
    """
    for d in get_cli_search_dirs():
        full_filename = os.path.join(d, fname)
        if os.path.exists(full_filename):
            return full_filename
    raise Exception('[WARN] Could not find {}; skipping.'.format(fname))


def add_prerequisite(prereqs, fname, suffix=None):
    """
    Return an updated set containing the prerequisites/dependencies for the
    file being processed.

    prereqs is the set containing the prequisites.

    fname is the filename of the new prerequisite to add to the set.

    suffix is an optional argument used to overwrite the suffix of fname.
    """
    if suffix:
        prereqs.add(basename_new_suffix(fname, suffix))
    else:
        prereqs.add(os.path.basename(fname))
    return prereqs


def find_module_file(desired_mod):
    """
    Return the name of the file containing the module specified by desired_mod.

    desired_mod is tne name of the module to find.

    This function searches all '*.F90' and '*.f90' files looking for the one
    that contains desired_mod.

    If desired_mod is in the exclude_fmod_files set, then simply return None.
    We do not want to add external modules, like ESMF, as a dependency.

    This function first tries to find desired_mod in a file named
    desired_mod.F90 or an upper or lower-case version of that name.
    If that search is unsuccessful, then this routine searches all Fortran90
    files.
    """
    #exclude_fmod_files = {'esmf', 'netcdf', 'grib_api', 'hdf5', 'mpi'}
    exclude_fmod_files = get_cli_exclude_files('fmod')
    if desired_mod.lower() in exclude_fmod_files:
        return None
    check_files = [desired_mod+'.F90', desired_mod+'.f90',
                   desired_mod.lower()+'.F90', desired_mod.lower()+'.f90',
                   desired_mod.upper()+'.F90', desired_mod.upper()+'.f90']
    for d in get_cli_search_dirs():
        for cf in check_files:
            filename = os.path.join(d, cf)
            if os.path.isfile(filename):
                if contains_module_definition(filename, desired_mod):
                    return cf
    for d in get_cli_search_dirs():
        f90s = os.path.join(d, '*.f90')
        F90s = os.path.join(d, '*.F90')
        check_files = glob.glob(f90s)
        check_files += glob.glob(F90s)
        for cf in check_files:
            if contains_module_definition(cf, desired_mod):
                return cf
    print_warn('module {} not found.'.format(desired_mod))
    return None


def contains_module_definition(file_name, desired_mod):
    """
    Return whether file_name contains the module definition for desired_mod.
    """
    desired_mod = desired_mod.lower()
    with input_open(file_name, 'r', errors='replace') as f:
        for line in f:
            result = module_def_pattern.match(line)
            if result:
                use_mod = result.group(1).lower()
                if use_mod == desired_mod:
                    return True
        else:
            return False


def find_use_module_statement(line):
    """
    Determine whether line contains a Fortran USE statement.  If so,
    find the file that contains the definition of the module being USE'd
    and return that file's name.

    line is a line of code from the Fortran file being processed.
    """
    result = use_statement_pattern.match(line)
    if result:
        use_mod = result.group(1)
        return find_module_file(use_mod)
    else:
        return None


def find_f_include_statement(line):
    """
    Determine whether line contains a Fortran INCLUDE statement.  If so,
    return the file name being INCLUDE'd.

    line is a line of code from the Fortran file being processed.

    If the INCLUDE'd file is in the exclude_finc_files set, then simply
    return None.  We do not want to add external include files, like mpif.h,
    as a prerequisite.
    """
    #exclude_finc_files = {'mpif.h', 'netcdf.inc', 'CXMLDEF.FOR'}
    exclude_finc_files = get_cli_exclude_files('finc')
    result = f_include_pattern.match(line)
    if result:
        inc_file = result.group(2)
        if inc_file in exclude_finc_files:
            return None
        else:
            return inc_file
    else:
        return None


def find_c_include_file(desired_file):
    """
    Determine whether the header file desired_file is in the list
    of directories specified by get_cli_search_dirs.

    line is a line of code from the file being processed.

    This function searches all directories specified by get_cli_search_dirs
    to find the file being INCLUDE'd.  If the file cannot be found, then
    it is assumed to be a <system> header file, like <string.h>, and it
    should not be added to the dependency list.
    """
    for d in get_cli_search_dirs():
        filename = os.path.join(d, desired_file)
        if os.path.isfile(filename):
            return desired_file
    print_info('excluding {} from dependency list'.format(desired_file))
    return None


def find_c_include_statement(line):
    """
    Determine whether line contains a C '#include' statement.  If so,
    return the file name being INCLUDE'd.

    line is a line of code from the file being processed.

    This function searches calls find_c_include_file to find the file
    being INCLUDE'd.  If the file cannot be found, then it is assumed
    to be a <system> header file, like <string.h>, and it should not be
    added to the dependency list.
    """
    result = c_include_pattern.match(line)
    if result:
        header_file = result.group(1)
        return find_c_include_file(header_file)
    else:
        return None


def process_fortran90_file(fname, prereqs):
    """
    Return the set of prerequisites/dependencies for file fname.

    fname is the file name of the Fortran90 file being processed.

    prereqs is a set containing the prerequisites/dependencies for file fname.
    """
    try:
        ffname = find_file(fname)
        f = input_open(ffname, 'r', errors='replace')
    except IOError as e:
        print_err("Cannot open file '{}'; skipping.".format(fname))
    except Exception as e:
        print(e)
    else:
        for line in f:
            name = find_use_module_statement(line)
            if name:
                prereqs = add_prerequisite(prereqs, name, '.o')
            name = find_f_include_statement(line)
            if name:
                prereqs = add_prerequisite(prereqs, name)
                print_dbg('recursing')
                prereqs = process_fortran90_file(name, prereqs)
            name = find_c_include_statement(line)
            if name:
                prereqs = add_prerequisite(prereqs, name)
                print_dbg('recursing')
                prereqs = process_fortran90_file(name, prereqs)
        f.close()
    return prereqs


def process_fortran77_file(fname, prereqs):
    """
    Return the set of prerequisites/dependencies for file fname.

    fname is the file name of the Fortran90 file being processed.

    prereqs is a set containing the prerequisites/dependencies for file fname.
    """
    try:
        ffname = find_file(fname)
        f = input_open(ffname, 'r', errors='replace')
    except IOError as e:
        print_err("Cannot open file '{}'; skipping.".format(fname))
    except Exception as e:
        print(e)
    else:
        for line in f:
            name = find_f_include_statement(line)
            if name:
                prereqs = add_prerequisite(prereqs, name)
                prereqs = process_fortran77_file(name, prereqs)
            name = find_c_include_statement(line)
            if name:
                prereqs = add_prerequisite(prereqs, name)
                prereqs = process_fortran77_file(name, prereqs)
        f.close()
    return prereqs


def process_c_file(fname, prereqs):
    """
    Return the set of prerequisites/dependencies for file fname.

    fname is the file name of the Fortran90 file being processed.

    prereqs is a set containing the prerequisites/dependencies for file fname.
    """
    try:
        ffname = find_file(fname)
        f = input_open(ffname, 'r', errors='replace')
    except IOError as e:
        print_err("Cannot open file '{}'; skipping.".format(fname))
    except Exception as e:
        print(e)
    else:
        for line in f:
            name = find_c_include_statement(line)
            if name:
                prereqs = add_prerequisite(prereqs, name)
        f.close()
    return prereqs


def write_prerequisites(depfile, base_filename, target, prereqs):
    """
    Write the set of prerequisites/dependencies for file base_filename
    into the file pointed to be depfile.

    depfile is a file handle to the dependency file to write
    (base_filename.d).

    base_filename is the of the file being processed.

    target is the target of the prerequisites (base_filename.o).

    prereqs is a set containing the prerequisites/dependencies
    for file base_filename.
    """
    depfile.write('{} {} : {}\n'.format(
        target,
        basename_new_suffix(target, '.d'),
        base_filename))
    for p in prereqs:
        depfile.write('{} : {}\n'.format(target, p))


def process_command_line():
    import argparse

    description = "Scan files to generate compile-time dependencies."
    epilog = """
   Notes:
   Files with suffixes '.f90' or '.F90' are considered Fortran90 source files.
   Files with suffixes '.f' or '.F' are considered Fortran77 source files.
   Files with suffixes '.c' or '.cc' are considered C source files.
   Scanning files with any other suffix is considered an error.

   No pre-processing is done to the files before scanning.

   This program assumes that there is only one Fortran "use module_name"
   statement in any given line, and it assumes that the "use" key-word and
   the module name are in the same line.  Similarly for Fortran
   "include file_name" statements.
   """

    parser = argparse.ArgumentParser(description=description, epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose',
                        action='store', dest='verbose', default=None,
                        choices=['status', 'info', 'warn', 'dbg'],
                        help='enable status messages (default: status messages are disabled)')

    parser.add_argument('--exclude-fmod-files',
                        action='store', dest='exclude_fmod_files',
                        type=lambda s: {i.strip() for i in s.split(',')},
                        default=set(),
                        help="comma-separated list of Fortran module files to exclude "
                        "from the dependency scan; e.g. 'esmf,netcdf'")

    parser.add_argument('--exclude-finc-files',
                        action='store', dest='exclude_finc_files',
                        type=lambda s: {i.strip() for i in s.split(',')},
                        default=set(),
                        help="comma-separated list of Fortran include files to exclude "
                        "from the dependency scan; e.g. 'mpif.h,netcdf.inc'")

    parser.add_argument('files',
                        action='store',
                        type=lambda s: {i.strip() for i in s.split(':')},
                        help='colon-separated list of files to scan')

    parser.add_argument('dirs',
                        action='store',
                        type=lambda s: {i.strip() for i in s.split(':')},
                        help='colon-separated list of directories to search '
                        'for dependencies (Fortran modules, include files, etc.)')

    cli_args = parser.parse_args()

    cli_args.STATUS = 1
    cli_args.INFO = 2
    cli_args.WARN = 3
    cli_args.DBG = 4
    if cli_args.verbose == 'status':
        cli_args.VERBOSE = cli_args.STATUS
    elif cli_args.verbose == 'info':
        cli_args.VERBOSE = cli_args.INFO
    elif cli_args.verbose == 'warn':
        cli_args.VERBOSE = cli_args.WARN
    elif cli_args.verbose == 'dbg':
        cli_args.VERBOSE = cli_args.DBG
    else:
        cli_args.VERBOSE = 0
    return cli_args


#
# main
#

# regular expressions for matching module and include statements
module_def_pattern = re.compile(r'\A\s*module\s+(\w+)', re.IGNORECASE)
use_statement_pattern = re.compile(r'\A\s*use\s+(\w+)', re.IGNORECASE)
f_include_pattern = re.compile(r'\A\s*include\s+([\'"])(.*?)\1', re.IGNORECASE)
c_include_pattern = re.compile(
    r'\A\s*#\s*include\s+[\'"<](.*?)[\'">]', re.IGNORECASE)

cli_args = process_command_line()
print_dbg('cli_args: ', cli_args)

files = cli_args.files

for full_filename in files:
    print_status('Generating dependencies for {}'.format(full_filename))
    # derivatives of full_filename
    base_filename = os.path.basename(full_filename)
    depname = basename_new_suffix(base_filename, '.d')
    target = basename_new_suffix(base_filename, '.o')

    # open and initialize the dependency list
    depfile = output_open(depname, 'w')
    prereqs = set()

    try:
        suffix = suffix_check(base_filename)
    except Exception as e:
        print(e)
    else:
        if suffix == 'fortran90':
            prereqs = process_fortran90_file(full_filename, prereqs)
        elif suffix == 'fortran77':
            prereqs = process_fortran77_file(full_filename, prereqs)
        elif suffix == 'c':
            prereqs = process_c_file(full_filename, prereqs)
        # remove circular references; see below.
        prereqs.discard(target)
        write_prerequisites(depfile, base_filename, target, prereqs)
    finally:
        depfile.close()

#
# A circular reference such as file.o depends on file.o can occur when
# file.F90 contains two module definitions, say for module A and module B,
# and module B USEs module A.
#
