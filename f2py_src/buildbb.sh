#!/bin/bash

module purge; module load bluebear
module load bear-apps/2018a
module load imkl 2019.5.281-gompi-2019b
module load Python 3.7.4-GCCcore-8.3.0
module load matplotlib 3.1.1-foss-2019b-Python-3.7.4

rm -R /rds/homes/b/bxp912/apps/f2py_lib/eampa_lib
mkdir -p /rds/homes/b/bxp912/apps/f2py_lib/eampa_lib

top=$("pwd")
echo $top
for d in */ ; do
  echo "####################################"
  echo $d
  echo "####################################"
  cd $d
  if test -f "build.sh"; then
    ./build.sh
    cp *.so /rds/homes/b/bxp912/apps/f2py_lib/eampa_lib
  fi
  cd $top
done

