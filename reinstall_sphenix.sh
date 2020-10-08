#!/bin/bash

source /opt/sphenix/core/root/bin/thisroot.sh &&
source /opt/sphenix/core/bin/sphenix_setup.sh -n new &&
export SPHENIX=$PWD &&
export MYINSTALL="$SPHENIX"/install &&
mkdir -p $MYINSTALL &&
mkdir -p $MYINSTALL/include &&
mkdir -p $MYINSTALL/lib &&
mkdir -p $MYINSTALL/bin &&
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL &&

cd $PWD/treeMaker/
source autogen.sh --prefix=$MYINSTALL &&
make
make install
