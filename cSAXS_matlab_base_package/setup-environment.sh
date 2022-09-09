host_full=`hostname`
host=`echo $host_full | awk -F'-' '{print $1}'`
host_node=`echo $host_full | awk -F'-' '{print $2}'`

module clear

run_environment=1

if [ $host == 'x12sa' ]; then
	if [ ! $host_node == 'gpu' ] && [ ! $host_node == 'cn' ] && [ ! $host_node == 'cons' ] && [ ! $host_node == 'comp' ];then
		echo 'The cSAXS environment is not available for' $host_full
		run_environment=0
	fi
fi

if [ $run_environment == 1 ]; then

ptycho=1
tomo=1
SAXS=1

if [ "$#" -gt 0 ]; then
    ptycho=0;
    tomo=0;
    SAXS=0;
    for i in "$@"; do 
            if [ $i == "tomo" ]; then
                tomo=1
            fi
            if [ $i == "ptycho" ]; then
                ptycho=1
            fi
	    if [ $i == "SAXS" ]; then
                SAXS=1
	    fi
    done
fi

source /opt/psi/config/profile.bash

if [ "$tomo" == 1 ] || [ "$ptycho" == 1 ]; then
module add cuda/10.0.130
fi

module add gcc/8.2.0  # needed for openmpi/3.0.1 and hdf5_serial/1.10.1 , mexcuda in matlab needs  GCC 4.8.5 !!
# force matlab to use our shared libraries
export LD_PRELOAD=/opt/psi/Programming/gcc/8.2.0/lib64/libstdc++.so.6:$LD_PRELOAD

module add hdf5_serial/1.10.5 # needed for data storage and loading 

if [ "$ptycho" == 1 ]; then
export MKLROOT=/opt/psi/Programming/intel/18.2/mkl
LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH 
module add openmpi/3.1.2

#export LD_LIBRARY_PATH=$pa1/cxs_shared_lib/:$LD_LIBRARY_PATH
if [[ -z "${DISPLAY// }" ]]; then
  echo 'No display? Okay...'
  alias matlab="matlab -nodesktop -nosplash"
fi
fi

module add matlab/2018a
fi

