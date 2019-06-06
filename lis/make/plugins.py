#!/usr/bin/env python
from __future__ import print_function
import StringIO
import ConfigParser

'''
   This program processes the default configuration file, default.cfg,
   and the user-modifiable file, user.cfg, to enable/disable the various
   plugins components within LIS.

   This program creates two files, Filepath and LIS_plugins.h,
   which are used by the Makefile and preprocessor to compile the
   LIS source code.

   By default, all optional plugins are enabled and all restricted and
   old/unsupported plugins are disabled.  Note that the restricted plugins
   are not distributed with the public release of the LIS source code.
'''


warning_msg1 = '''\
#
#           !!! WARNING !!!
# This is a generated file.  Do not manually edit.
# To enable/disable optional components, please create a user.cfg file.
# See the LIS Users' Guide for more informaion.
#

'''


warning_msg2 = '''\
#if 0

          !!! WARNING !!!
This is a generated file.  Do not manually edit.
To enable/disable optional components, please create a user.cfg file.
See the LIS Users' Guide for more informaion.

#endif
'''


#
# Define the required components
#
required_filepath='''\
core
plugins
interp
'''


#
# Process sections and enable/disable table
#
config = ConfigParser.SafeConfigParser()
config.read('default.cfg')

try:
   for line in open('user.cfg', 'r'):
      if '#' in line:
         line = line[0:line.find('#')]
      if line.strip():
         sect, enabled = (l.strip() for l in line.split(':'))
         config.set(sect, 'enabled', enabled)
except IOError:
   print('No user.cfg file found.  Using defaults.')

# Create virtual sections
config.add_section('virtual_irrigation')
config.set('virtual_irrigation', 'enabled', 'false')
config.set('virtual_irrigation', 'virtual', 'true')
if ( config.getboolean('Sprinkler', 'enabled') or
     config.getboolean('Flood', 'enabled') or
     config.getboolean('Drip', 'enabled') ):
      config.set('virtual_irrigation', 'enabled', 'true')

config.add_section('virtual_routing')
config.set('virtual_routing', 'enabled', 'false')
config.set('virtual_routing', 'virtual', 'true')
if ( config.getboolean('NLDAS router', 'enabled') or
     config.getboolean('HYMAP router', 'enabled') ):
      config.set('virtual_routing', 'enabled', 'true')

config.add_section('virtual_da')
config.set('virtual_da', 'enabled', 'false')
config.set('virtual_da', 'virtual', 'true')
if ( config.getboolean('Direct insertion', 'enabled') or
     config.getboolean('EnKF', 'enabled') or
     config.getboolean('EnKS', 'enabled') ):
      config.set('virtual_da', 'enabled', 'true')

config.add_section('virtual_optue')
config.set('virtual_optue', 'enabled', 'false')
config.set('virtual_optue', 'virtual', 'true')
if ( config.getboolean('OPTUE ES', 'enabled') or
     config.getboolean('OPTUE LM', 'enabled') or
     config.getboolean('OPTUE GA', 'enabled') or
     config.getboolean('OPTUE SCEUA', 'enabled') or
     config.getboolean('OPTUE MCSIM', 'enabled') or
     config.getboolean('OPTUE RWMCMC', 'enabled') or
     config.getboolean('OPTUE DEMC', 'enabled') or
     config.getboolean('OPTUE DEMCz', 'enabled') ):
      config.set('virtual_optue', 'enabled', 'true')

config.add_section('virtual_da_obs_snodep')
config.set('virtual_da_obs_snodep', 'enabled', 'false')
config.set('virtual_da_obs_snodep', 'virtual', 'true')
if ( config.getboolean('virtual_da', 'enabled') and
     config.getboolean('DA OBS SNODEP', 'enabled') ):
      config.set('virtual_da_obs_snodep', 'enabled', 'true')

config.add_section('virtual_da_obs_ldtsi')
config.set('virtual_da_obs_ldtsi', 'enabled', 'false')
config.set('virtual_da_obs_ldtsi', 'virtual', 'true')
if ( config.getboolean('virtual_da', 'enabled') and
     config.getboolean('DA OBS LDTSI', 'enabled') ):
      config.set('virtual_da_obs_ldtsi', 'enabled', 'true')

#
# Write Filepath and LIS_plugins.h
#
filepath = open('Filepath', 'w')
defineh = open('LIS_plugins.h', 'w')

filepath.write(warning_msg1)
defineh.write(warning_msg2)

fpath = ' ../{0}'
dentry = '#{0} {1}\n'

filepath.write('dirs := .')
paths = StringIO.StringIO(required_filepath)
for line in paths:
   filepath.write(fpath.format(line.strip()))

for sect in config.sections():
   if not config.has_option(sect, 'virtual'):
      if config.getboolean(sect, 'enabled'):
         defineh.write(dentry.format('define', config.get(sect, 'macro')))
         for p in config.get(sect, 'path').split(','):
            filepath.write(fpath.format(p.strip()))
         if config.has_option(sect, 'dependent_comps'):
            dcomps = config.get(sect, 'dependent_comps')
            if dcomps:
               for d in dcomps.split(','):
                  d = d.strip()
                  if config.getboolean(d, 'enabled'):
                     d_path = d+' path'
                     for p in config.get(sect, d_path).split(','):
                        filepath.write(fpath.format(p.strip()))
      else:
         defineh.write(dentry.format('undef', config.get(sect, 'macro')))

filepath.close()
defineh.close()
