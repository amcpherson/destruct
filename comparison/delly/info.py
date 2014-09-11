import os

config = {}
execfile('setup.config', config) 

comparison_directory = config['comparison_directory']

data_directory = os.path.join(comparison_directory, 'data')
install_directory = os.path.join(comparison_directory, 'install')
results_directory = os.path.join(comparison_directory, 'results')

