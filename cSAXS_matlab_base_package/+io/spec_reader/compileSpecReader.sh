module add gcc/8.2.0
module add hdf5_serial/1.10.5
module add openmpi/3.1.2

export LIBRARY_PATH=$(pwd)/../../cxs_shared_lib/:$LIBRARY_PATH
export CPATH=$(pwd)/../../cxs_shared_lib/cxsLib/includes/:$(pwd)/../../cxs_shared_lib/tinyxml/includes/:$(pwd)/../../cxs_shared_lib/exprtk/includes/:$(pwd)/../../cxs_shared_lib/json/includes/:$CPATH

make init
make cleanAll
make
