#!/bin/bash
echo "Python code packaged"
python3 pack/pack.py eampa.py

bindir=$HOME"/pybin"
echo "Make bin dir: "$bindir
mkdir -p $bindir

libdir=$HOME"/pylib"
echo "Make lib dir: "$libdir
mkdir -p $libdir


eampalibdir=$HOME"/pylib/eampa_lib"
echo "Make eampa lib dir: "$eampalibdir
mkdir -p $eampalibdir


top=$("pwd")
cd f2py_src
srcdir=$("pwd")
d="bp"
cd $d
if test -f "build.sh"; then
  ./build.sh
  cp *.so $eampalibdir
fi
cd $srcdir


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
