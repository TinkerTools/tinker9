#!/usr/bin/env python

'''Run this Python script inside a git-enabled Tinker9 directory. \
Internet connection is mandatory. \
GitHub username and authentication may be required. \
Do not abuse this script because GitHub.com limits the rate of API requests.'''

import datetime
import json
import requests
from requests.auth import HTTPBasicAuth
import subprocess
import sys

tinker9_url = 'https://github.com/TinkerTools/tinker9'
tinker9_apiurl = 'https://api.github.com/repos/tinkertools/tinker9'
iso_8601_format = '%Y-%m-%dT%H:%M:%SZ'
git_log_format = '''TZ=UTC git log -1 --date='format-local:%Y-%m-%dT%H:%M:%SZ' \
--pretty=format:'{%n\
"commit": "%H",%n\
"abbreviated_commit": "%h",%n\
"tree": "%T",%n\
"author": {"name": "%aN","email": "%aE","date": "%ad"},%n\
"commiter": {"name": "%cN","email": "%cE","date": "%cd"}%n}'\
'''

def print_timedelta(delta:datetime.timedelta):
   sec = delta.seconds
   h = sec // 3600
   m = (sec - 3600*h) // 60
   s = sec - 3600*h - 60*m
   msg = '{}:{}:{}'.format(h, m ,s)
   if delta.days == 1:
      msg = '1 day {}'.format(msg)
   elif delta.days > 1:
      msg = '{} days {}'.format(delta.days, msg)
   return msg

def get_github_head_json():
   try:
      r = requests.get(tinker9_apiurl+'/commits/HEAD', auth=HTTPBasicAuth('user', 'pass'))
      return r.json()
   except:
      print('Cannot access GitHub API. Aborting...')
      sys.exit(1)

def get_libtinker_json(version:str='HEAD'):
   try:
      r = requests.get(tinker9_apiurl+'/contents/tinker?ref={}'.format(version), auth=HTTPBasicAuth('user', 'pass'))
      return r.json()
   except:
      print('Cannot access GitHub API. Aborting...')
      sys.exit(1)

def exec_command(cmd:str):
   sp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
   out = sp.communicate()[0]
   if sp.returncode != 0:
      print('Command -- {} -- returned error code {}.'.format(cmd, sp.returncode))
      sys.exit(sp.returncode)
   else:
      return out.decode('utf-8').strip()

def check_tinker9():
   print('Checking Tinker9...')
   # GitHub Tinker9
   d0 = get_github_head_json()
   datestr = d0['commit']['author']['date']
   t0 = datetime.datetime.strptime(datestr, iso_8601_format)
   # GitHub libtinker
   d0 = get_libtinker_json()
   sha8 = d0['sha']

   # local version
   out = exec_command(git_log_format)
   j = json.loads(out)
   sha1 = d0['sha']
   datestr = j['author']['date']
   t1 = datetime.datetime.strptime(datestr, iso_8601_format)
   if t1 == t0:
      print('The current version {} is up-to-date.'.format(sha1))
   elif t1 < t0:
      print('The current version {} is {} behind {}.'.format(sha1, print_timedelta(t0-t1), tinker9_url))
      # local libtinker
      print('Checking libtinker...')
      sha9 = exec_command(r"git submodule status | grep 'tinker (' | awk '{print $1}'")
      if sha9 == '':
         print('Submodule tinker is not initialized. Exit.')
      elif sha8 != sha9:
         print('libtinker needs recompilation after update.')
      else:
         print('libtinker is up-to-date.')
   else:
      print('The current version {} is {} ahead of {}.'.format(sha1, print_timedelta(t1-t0), tinker9_url))

if __name__ == '__main__':
   check_tinker9()
