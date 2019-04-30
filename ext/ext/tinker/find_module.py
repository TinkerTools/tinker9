#!/usr/bin/env python


'''
e.g.
(1) ./find_module.py list [fortran_files] | bash
'''

import os
import sys

Option = ''
Filenames = []
if len(sys.argv) > 1:
    Option = sys.argv[1]
if len(sys.argv) > 2:
    Filenames = sys.argv[2:]


MODULE_FILES = []
PROGRAM_FILES = []
SUBROUTINE_FILES = []

UNKNOWN_TYPE = 'UNKNOWN_TYPE'
MODULE_TYPE = 'MODULE_TYPE'
PROGRAM_TYPE = 'PROGRAM_TYPE'
SUBROUTINE_TYPE = 'SUBROUTINE_TYPE'


def determine_module_subroutine_program(fortran_filename):
    global MODULE_FILES
    global PROGRAM_FILES
    global SUBROUTINE_FILES

    # 'source/dynamic.f' -> 'dynamic.f'
    dirname = os.path.dirname(fortran_filename)
    base = os.path.basename(fortran_filename)
    stem, ext = os.path.splitext(base)  # 'dynamic.f' -> 'dynamic', '.f'
    content = [line.rstrip().lower() for line in open(fortran_filename)]

    use_list = []
    file_type = UNKNOWN_TYPE

    def type_must_be_unknown(t):
        if t != UNKNOWN_TYPE:
            raise BaseException('Cannot parse file: %s' % fortran_filename)

    for line in content:
        # '      module foo'
        if len(line) > 13:
            if line[0:13] == '      module ':
                type_must_be_unknown(file_type)
                file_type = MODULE_TYPE
                MODULE_FILES.append(stem)
        # '      program foo'
        if len(line) > 14:
            if line[0:14] == '      program ':
                type_must_be_unknown(file_type)
                file_type = PROGRAM_TYPE
                PROGRAM_FILES.append(stem)
        # '      use foo'
        if len(line) > 10:
            if line[0:10] == '      use ':
                temp_list = line.split()
                if len(temp_list) > 1:
                    word = temp_list[1]
                    if word not in use_list:
                        use_list.append(word)
    if file_type == UNKNOWN_TYPE:
        file_type = SUBROUTINE_TYPE
        SUBROUTINE_FILES.append(stem)

    use_list.sort()
    return stem, file_type, use_list, dirname


def option_list(files):
    dirname = ''
    for f in files:
        _, _, _, dirname = determine_module_subroutine_program(f)

    text1 = '''#!/bin/bash
HEADER=tinker.mod.h
DIR=detail

rm -f $HEADER
rm -f $DIR/*.hh

##########

cat << ENDOFFILE >> $HEADER
#ifndef TINKER_MOD_H_
#define TINKER_MOD_H_

ENDOFFILE
'''

    text2 = '''
##########

cat << ENDOFFILE >> $HEADER

#endif
ENDOFFILE
'''
    print(text1)
    for m in MODULE_FILES:
        print('./parse.py fpp %s/%s.f > $DIR/%s.hh' % (dirname, m, m))
    for m in MODULE_FILES:
        print('echo \'#include "detail/%s.hh"\' >> $HEADER' % m)
    print(text2)


if __name__ == '__main__':
    if Option.lower() == 'list':
        option_list(Filenames)
