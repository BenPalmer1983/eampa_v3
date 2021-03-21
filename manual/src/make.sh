#!/bin/bash
#################################################################################
# Bash Script
#
#################################################################################
# Load arguments to array
args=("$@")
# Change to directory script is being run from:
# cd ${0%/*}
# Set useful variables
scriptName="Bash Template"
workingDir=$(pwd)
currentUser=$(whoami)
printHeading=0
# Print
if [ $printHeading == 1 ]; then
  echo "Script:       "$scriptName
  echo "Working dir:  "$workingDir
  echo "User:         "$currentUser
fi
#################################################################################
#
# Cleanup
#########################
rm *.aux
rm *.xml
rm *.toc
rm *blx.bib
rm *.blg
rm *.log
#
# Make plots
#########################
for directory in $workingDir/*/ ; do
  for filename in $directory*.plot; do
    if [ -f $filename ]; then
      echo $filename
      cd $directory
      gnuplot $filename
      cd $workingDir
    fi
  done
  for subdirectory in $directory*/ ; do
    for filename in $subdirectory*.plot; do
      if [ -f $filename ]; then
        echo $filename
        cd $subdirectory
        gnuplot $filename
        cd $workingDir
      fi
    done
  done
done
# Build first time
pdflatex main.tex
# Make bibliography
bibtex8 main
# Build with bibliography
pdflatex main.tex
cp main.pdf ../../manual.pdf

#################################################################################
