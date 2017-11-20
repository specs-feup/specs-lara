#!/bin/bash


# make sure to have margot
MARGOT_FOLDER_NAME=$PWD/autotuner
BUILD_FOLDER_NAME=$MARGOT_FOLDER_NAME/build
INSTALL_FOLDER_NAME=$MARGOT_FOLDER_NAME/install
$PWD/bootstrap.sh


# get a copy of the high level interface (removing the previous one, if any)
SOURCE_MARGOT_HEEL=$MARGOT_FOLDER_NAME/margot_heel/margot_heel_if
DEST_MARGOT_HEEL=$PWD/high_level_interface
rm -rf $DEST_MARGOT_HEEL
cp -r $SOURCE_MARGOT_HEEL $DEST_MARGOT_HEEL


# remove the example configuration files
rm $DEST_MARGOT_HEEL/config/*.conf


# copy the cconfiguration file (is the same for autotuning, but w/out operating points)
cp $PWD/config/adaptation.conf $DEST_MARGOT_HEEL/config
if ! [ -f $PWD/config/oplist.conf ]; then
	cp $PWD/config/oplist.conf $DEST_MARGOT_HEEL/config
fi


# build the high level interface of margot
PREVIOUS_WORKING_DIR=$PWD
mkdir $DEST_MARGOT_HEEL/build
cd $DEST_MARGOT_HEEL/build
cmake -DCMAKE_BUILD_TYPE=Release -DMARGOT_CONF_FILE=adaptation.conf .. || exit -1
make || exit -1
cd $PREVIOUS_WORKING_DIR


# build the actual application
PREVIOUS_WORKING_DIR=$PWD
mkdir -p $PWD/build
cd $PWD/build
cmake -DCMAKE_BUILD_TYPE=Release .. || exit -1
make || exit -1
cd $PREVIOUS_WORKING_DIR
