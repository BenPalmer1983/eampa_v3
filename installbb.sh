#!/bin/bash
module purge; module load bluebear
module load bear-apps/2018a
module load imkl 2019.5.281-gompi-2019b
module load Python 3.7.4-GCCcore-8.3.0
module load matplotlib 3.1.1-foss-2019b-Python-3.7.4

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
  if test -f "build.sh"; then
    ./build.sh
    cp *.so $libdir
  fi
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













