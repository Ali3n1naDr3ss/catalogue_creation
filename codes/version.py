


import numpy
print("version: ",numpy.__version__)
print("file: ", numpy.__file__)

export PYTHONPATH=/home/ullyott/.local/lib/python3.11/site-packages:$PYTHONPATH
python3 -c "import numpy; print(numpy.__version__)"
