#!/bin/bash

# clone the autotuner framework
MARGOT_FOLDER_NAME=$PWD/autotuner
if ! [ -d $MARGOT_FOLDER_NAME ]; then
  echo "[INFO] clone mARGOt"
  git clone https://gitlab.com/margot_project/core.git $MARGOT_FOLDER_NAME
fi


# build the framework
BUILD_FOLDER_NAME=$MARGOT_FOLDER_NAME/build
INSTALL_FOLDER_NAME=$MARGOT_FOLDER_NAME/install
if ! [ -d $INSTALL_FOLDER_NAME ]; then
  echo "[INFO] build mARGOt"
  mkdir $BUILD_FOLDER_NAME
  PREVIOUS_WORKING_DIR=$PWD
  cd $BUILD_FOLDER_NAME
  cmake -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL_FOLDER_NAME -DCMAKE_BUILD_TYPE=Release .. || exit -1
  make install
  cd $PREVIOUS_WORKING_DIR
fi
