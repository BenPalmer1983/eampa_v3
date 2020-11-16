#!/bin/bash

mkdir -p ../../f2py_lib/eampa_lib/
top=$("pwd")
echo $top
for d in */ ; do
  echo "####################################"
  echo $d
  echo "####################################"
  cd $d
  if test -f "build.sh"; then
    ./build.sh
    cp *.so ../../f2py_lib/eampa_lib/
  fi
  cd $top
done
