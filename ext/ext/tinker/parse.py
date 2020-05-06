#!/usr/bin/env python

import sys
import re
import os

_Option = ''
_Filename = ''
_FilenameStem = ''
_RawContents = ''

if len(sys.argv) > 1:
    _Option = sys.argv[1]
if len(sys.argv) > 2:
    _Filename = sys.argv[2]  # /usr/local/foo.f
    _FilenameStem = os.path.basename(_Filename)  # foo.f
    _FilenameStem = os.path.splitext(_FilenameStem)[0]  # foo
    _RawContents = [line.rstrip().lower() for line in open(_Filename)]

_ModuleKeys = []
_Modulus = {}
_CurrentModule = ''

# _SharedHeader = '"util/macro.h"'
# _SharedHeader = '"util_macro.h"'
_SharedHeader = '"macro.h"'

_TINKER_MOD = 'TINKER_MOD'

_KnownTypes = {
    'integer': 'int',
    'integer*8': 'unsigned long long',  # used as memory handle
    'logical': 'int',
    'real*8': 'double',
    'real*4': 'float',
    'character': 'char',
    'character*1': 'char',
}
for i in range(2, 1025):
    _KnownTypes['character*%d' % i] = 'char'

_CXXKeywords = [
    'alignas', 'alignof', 'and', 'and_eq', 'asm',
    'atomic_cancel', 'atomic_commit', 'atomic_noexcept', 'auto',
    'bitand', 'bitor', 'bool', 'break',
    'case', 'catch', 'char', 'char8_t', 'char16_t', 'char32_t',
    'class', 'compl', 'concept', 'const', 'consteval', 'constexpr',
    'const_cast', 'continue', 'co_await', 'co_return', 'co_yield',
    'decltype', 'default', 'delelte', 'do', 'double', 'double_cast',
    'else', 'enum', 'explicit', 'export', 'extern',
    'false', 'float', 'for', 'friend',
    'goto',
    'if', 'import', 'inline', 'int',
    'long',
    'module', 'mutable',
    'namespace', 'new', 'noexcept', 'not', 'not_eq', 'nullptr',
    'operator', 'or', 'or_eq',
    'private', 'protected', 'public',
    'reflexpr', 'register', 'reinterpret_cast', 'requires', 'return',
    'short', 'signed', 'sizeof', 'static', 'static_assert',
    'static_cast', 'struct', 'switch', 'synchronized',
    'template', 'this', 'thread_local', 'throw',
    'true', 'try', 'typedef', 'typeid', 'typename',
    'union', 'unsigned', 'using',
    'virtual', 'void', 'volatile',
    'wchar_t', 'while',
    'xor', 'xor_eq'
]


class entry:
    def __init__(self):
        self.module = '_Uninitialized'
        self.symbol = '_Uninitialized'
        self.value = '_Uninitialized'
        self.type = '_Uninitialized'
        self.is_const = False
        self.dimension = []

    def __str__(self):
        this_type = self.type
        if self.is_const:
            this_type = 'const ' + this_type
        line = '%s::%s = %s of %s %s;' % (
            self.module, self.symbol, self.value, this_type, str(self.dimension))
        return line

    def print_fpp(self, op):
        local_symbol = self.symbol
        line = '#error ERROR'
        if local_symbol in _CXXKeywords:
            local_symbol += '_'
        if self.is_const:
            # const int a = 100; <- integer, parameter :: a=100
            if op == 'header':
                line = 'const %s %s = %s;' % (
                    self.type, local_symbol, self.value)
            elif op == 'extern_c':
                line = ''
            elif op == 'define':
                line = ''
        elif len(self.dimension):
            if ':' not in self.dimension:
                # extern int (&a)[maxatm]; <- integer a(maxatm)
                # extern char (&b)[6][5][4][240]; <- character*240 b(4,5,6)
                if op == 'header':
                    line = 'extern %s (&%s)' % (self.type, local_symbol)
                elif op == 'extern_c':
                    line = 'extern "C" %s %s(%s, %s)' % (
                        self.type, _TINKER_MOD, self.module, self.symbol)
                elif op == 'define':
                    line = '%s (&%s)' % (self.type, local_symbol)
                for x in reversed(self.dimension):
                    line += '[%s]' % x
                if op == 'define':
                    line += ' = %s(%s, %s)' % (_TINKER_MOD,
                                               self.module, self.symbol)
                line += ';'
            elif self.dimension[0] == ':':
                # extern float*& c; <-> real*4, allocatable :: c(:,:)
                if op == 'header':
                    line = 'extern %s*& %s;' % (self.type, local_symbol)
                elif op == 'extern_c':
                    line = 'extern "C" %s* %s(%s, %s);' % (
                        self.type, _TINKER_MOD, self.module, self.symbol)
                elif op == 'define':
                    line = '%s*& %s = %s(%s, %s);' % (self.type,
                                                      local_symbol,
                                                      _TINKER_MOD,
                                                      self.module,
                                                      self.symbol)
            else:
                # extern char (*&d)[6]; <-> character*6, allocatable :: d(:,:)
                known_dimension = [x for x in self.dimension if x != ':']
                if op == 'header':
                    line = 'extern %s (*&%s)' % (self.type, local_symbol)
                elif op == 'extern_c':
                    line = 'extern "C" %s (*%s(%s, %s))' % (
                        self.type, _TINKER_MOD, self.module, self.symbol)
                elif op == 'define':
                    line = '%s (*&%s)' % (self.type, local_symbol)
                for x in reversed(known_dimension):
                    line += '[%s]' % x
                if op == 'define':
                    line += ' = %s(%s, %s)' % (_TINKER_MOD,
                                               self.module, self.symbol)
                line += ';'
        else:
            # extern int& n;
            if op == 'header':
                line = 'extern %s& %s;' % (self.type, local_symbol)
            elif op == 'extern_c':
                line = 'extern "C" %s %s(%s, %s);' % (
                    self.type, _TINKER_MOD, self.module, self.symbol)
            elif op == 'define':
                line = '%s& %s = %s(%s, %s);' % (
                    self.type, local_symbol, _TINKER_MOD, self.module, self.symbol)
        if line != '':
            print(line)


class module:
    def __init__(self):
        self.name = '_Uninitialized'
        self.has_symbol = False
        self.depends_on = []
        self.entries = []

    def find_entry_by_symbol(self, s):
        for e in self.entries:
            if e.symbol == s:
                return True, e
        return False, ''

    def print_fpp(self):
        local_symbol = self.name
        if local_symbol in _CXXKeywords:
            local_symbol += '_'
        if len(self.depends_on):
            for h in self.depends_on:
                print('#include "%s.hh"' % h)
        print('\nnamespace tinker { namespace %s {' % local_symbol)
        if len(self.depends_on):
            for h in self.depends_on:
                print('using namespace %s;' % h)
            print('')
        for x in self.entries:
            x.print_fpp(op='header')
        non_const = 0
        for x in self.entries:
            if not x.is_const:
                non_const += 1
        if non_const:
            print('\n#ifdef TINKER_MOD_CPP_')
            for x in self.entries:
                x.print_fpp(op='extern_c')
            print('')
            for x in self.entries:
                x.print_fpp(op='define')
            print('#endif')
        print('} }')


def parse_ampersand(raw_file):
    '''Handle the & and $ characters in the fixed F77 format.'''
    totl = len(raw_file)
    if totl == 0 or totl == 1:
        return raw_file
    # totl >= 2
    contents = []
    for line in _RawContents:
        if len(line) > 6 and line[5] in ['&', '$']:
            contents[-1] = contents[-1] + line[6:]
        else:
            contents.append(line)
    return contents


def split_by_comma_outside_parentheses(s):
    '''Example:
real*8,dimension(:,:),allocatable -> [real*8, dimension(:,:), allocatable]'''
    return re.split(r',\s*(?![^()]*\))', s)


def expressions_within_parentheses(s):
    '''Examples:
a -> a
b6(0:a,0:a) -> 0:a,0:a
b5(a,a),b6(0:a,0:a) -> a,a'''
    if '(' not in s:
        return s
    else:
        return re.search(r'\((.*?)\)', s).group(1)


def fortran_double_precision_expression(s):
    result = re.search(r'[+-]?(\d+(\.\d*)?|\.\d+)([dD][+-]?\d+)?', s)
    if result:
        s2 = result.group(0)
        s2 = s2.replace('d', 'e')
        return s2.replace('D', 'e')
    else:
        return s


def convert_range_to_length(e):
    '''Examples:
: -> :
100 -> 100
0:100 -> 101
a:100 -> 101-a
0:a -> a+1
100:100+a -> 100+a-99
a:2*a -> 1+2*a-a'''
    totlen = e
    if ':' == e:  # :
        pass
    elif ':' in e and len(e) > 1:
        r = e.split(':')
        back = r[1]
        iback = False
        nback = -1
        try:
            nback = int(back)
            iback = True
        except ValueError:
            pass
        begin = r[0]
        ibegin = False
        nbegin = -1
        try:
            nbegin = int(begin)
            ibegin = True
        except ValueError:
            pass
        if iback and ibegin:  # 0:100
            totlen = '%d' % (nback + 1 - nbegin)
        elif iback and not ibegin:  # a:100
            totlen = '%d-%s' % (nback + 1, begin)
        elif not iback and ibegin:  # 0:a, 100:100+a
            if nbegin == 0:
                totlen = '%s+1' % back
            else:
                totlen = '%s-%d' % (back, nbegin - 1)
        else:  # a:2*a
            totlen = '1+%s-%s' % (back, begin)
    else:  # 100
        pass
    return totlen


def parse_line(raw_line):
    global _ModuleKeys
    global _Modulus
    global _CurrentModule

    line = raw_line.lstrip()
    if (len(line) == 0):
        return

    is_comment = (line[0] == '!' or line[0] == '#' or raw_line[0] != ' ')
    word = line.split()[0]
    word1 = word.split(',')[0]
    keys = _KnownTypes.keys()
    if is_comment:
        pass
    elif 'implicit ' in line and ' none' in line:
        pass
    elif line == 'save':
        pass
    elif 'module ' in line:
        _CurrentModule = line.split()[-1]
        m = module()
        m.name = _CurrentModule
        _ModuleKeys.append(_CurrentModule)
        _Modulus[_CurrentModule] = m
    elif line == 'end':
        _CurrentModule = ''
    elif line[0:4] == 'use ':
        m = _Modulus[_CurrentModule]
        m.depends_on.append(line.split()[-1])
    elif word in keys or word1 in keys or word in ['parameter']:
        parse_line_1(line)
    else:
        print("#error I don't understand this line: " + line)
    return


def parse_line_1(line):
    global _Modulus
    global _CurrentModule

    m = _Modulus[_CurrentModule]

    part1 = ''
    part2 = ''
    if '::' in line:
        y = line.split('::')
        part1 = y[0]
        part2 = y[1]
    else:
        y = line.split()
        part1 = y[0]
        part2 = ''.join(y[1:])
    # clear whitespace
    type_aspects = ''.join(part1.split(' '))
    symbol_aspects = ''.join(part2.split(' '))

    # parse type_aspects
    type_array = split_by_comma_outside_parentheses(type_aspects)
    this_type = ''
    for x in type_array:
        if x in _KnownTypes.keys():
            this_type = _KnownTypes[x]  # e.g. character*7 -> char
            break

    this_is_const = False
    if 'parameter' in type_array:
        this_is_const = True

    type_dimension = []
    for x in type_array:
        if 'character*' in x:  # e.g. character*7
            y = x.split('*')
            yy = int(y[1])
            type_dimension.append(y[1])
    for x in type_array:  # e.g. dimension(2), dimension(:,:)
        if 'dimension(' in x:
            ranges = expressions_within_parentheses(x).split(',')
            for r in ranges:
                totlen = convert_range_to_length(r)
                type_dimension.append(totlen)

    # parse symbol_aspects
    symbol_array = split_by_comma_outside_parentheses(symbol_aspects)
    for sa in symbol_array:
        e = entry()
        is_an_update = False
        this_dimension = [d for d in type_dimension]
        if '=' in sa:
            s = expressions_within_parentheses(sa)
            tmp = s.split('=')
            this_symbol = tmp[0]
            this_value = tmp[-1]
            this_value = fortran_double_precision_expression(this_value)
            assert(this_is_const)
            if this_type != '':  # integer, parameter :: x=1000, y=a
                e.type = this_type
            else:  # parameter (a=1000)
                is_an_update = True
                found, eref = m.find_entry_by_symbol(this_symbol)
                if found:
                    e = eref
            e.value = this_value
            e.is_const = this_is_const
        else:
            # integer a
            # integer b(3,3)
            this_symbol = sa.split('(')[0]
            if '(' in sa:
                ranges = expressions_within_parentheses(sa).split(',')
                for r in ranges:
                    totlen = convert_range_to_length(r)
                    this_dimension.append(totlen)
        if not is_an_update:
            e.module = _CurrentModule
            e.symbol = this_symbol
            e.type = this_type
            e.is_const = this_is_const
            e.dimension = this_dimension
            m.entries.append(e)


def print_file_contents():
    '''Convert the Fortran source code (not the comments) to lower case.'''
    raw = [line.rstrip() for line in open(_Filename)]
    for line in raw:
        comment_line = False
        if (len(line)):
            first_nonwhite = len(line) - len(line.lstrip())
            c = line[first_nonwhite]
            if (c == '!' or c == '#' or line[0] != ' '):
                comment_line = True
        if comment_line:
            print(line)
        else:
            print(line.lower())


def print_fpp():
    # guard = 'TINKER_MOD_%s_HH_' % _FilenameStem.upper()
    # print('#ifndef %s' % guard)
    # print('#define %s' % guard)
    print('#pragma once')
    print('\n#include %s' % _SharedHeader)
    for k in _ModuleKeys:
        m = _Modulus[k]
        m.print_fpp()
    # print('\n#endif')
    return


foo0_f = r'''      module foo0
      implicit none
      integer, parameter :: a=100
      save
      end
'''

foo_f = r'''c     ######################
c     ##  COPYRIGHT INFO  ##
c     ######################
c
c
      module TYPENAME
      integer a0
      parameter (a0=10)
      logical a1,a2(a0),a3(10)
      character*240 a4,a5(4,5,6)
      real*4, allocatable :: a6(:,:)
      character*7 j1(a0:2*a0),j2(a0:30),
     &            j3(0:a0),j4(0:100),j5(100,100+a0)
      save
      end
c
c
c
      module STATIC_CAST
      use foo0
      implicit none
      integer, parameter :: class=100,auto=a
      real*8, parameter :: namespace=3.1415d0
      real*8, dimension(3,10) :: b1,b2(4,5)
      real*8, dimension(:,:), allocatable :: b3
      real*4 b5(a,a), b6(0:a,0:a)
      save
      end
'''

##########################################################################

if __name__ == '__main__':
    for raw_line in parse_ampersand(_RawContents):
        parse_line(raw_line)
    _Option = _Option.lower()
    if _Option == 'fpp':
        print_fpp()
    elif _Option == 'gen_test':
        print(foo0_f)
        print(foo_f)
