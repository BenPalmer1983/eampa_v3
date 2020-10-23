#!/bin/bash
python3 pack/pack.py eampa.py

libdir=$("pwd")"/f2py_lib/eampa_lib"
mkdir -p $libdir
echo $libdir

top=$("pwd")
cd f2py_src
srcdir=$("pwd")
for d in */ ; do
  echo "####################################"
  echo $d
  echo "####################################"
  cd $d
  #if test -f "build.sh"; then
  #  ./build.sh
  #  cp *.so $libdir
  #fi
  cd $srcdir
done

cd $top

export_line="export PYTHONPATH=\$PYTHONPATH:\""$libdir"\""
profile_file=$HOME"/.bash_profile"
touch $profile_file
if grep -q "$export_line" "$profile_file";
then
  echo $profileFile" already has the bin directory path"
else
  echo "Adding "$export_line" to "$profile_file
  echo $export_line >> $profile_file
fi













