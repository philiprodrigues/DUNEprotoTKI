echo calling DUNEprotoTKI/setup.sh 

export DUNEPROTOTKI=$(pwd)
export LOCALBIN=${DUNEPROTOTKI}/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${DUNEPROTOTKI}/style
export PATH=$PATH:${LOCALBIN}


