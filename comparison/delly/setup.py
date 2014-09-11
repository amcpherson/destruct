import os
import sys

conf_filename = 'setup.config'

if os.path.exists(conf_filename):
    print 'Configuration file ' + conf_filename + ' exists'
    sys.exit()

with open(conf_filename, 'w') as conf_file:
    conf_file.write('# Comparison configuration file\n\n')
    conf_file.write('# Root directory for comparison data\n')
    conf_file.write('comparison_directory = \n')

print 'Edit configuration file ' + conf_filename + ' to complete the setup.'
