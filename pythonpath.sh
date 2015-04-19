#
# Adds the required directories to the python path
#

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export PYTHONPATH=$DIR/pypeliner:$DIR/pygenes:$PYTHONPATH