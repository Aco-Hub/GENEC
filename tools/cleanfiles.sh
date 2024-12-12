#!/bin/bash
# usage ~/path/to/cleanfiles P1p1z14S0

DIR_SOURCE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" > /dev/null && pwd )"
ScriptBZ="python "$DIR_SOURCE"/../../UtilsEvol/BurningZonesCalc.py"
catag="python "$DIR_SOURCE"/catag.py"
wg_reductor="python "$DIR_SOURCE"/wg_reductor.py"

if [ `uname -s` = "SunOS" ];
then
  tarprog="/opt/csw/bin/gtar"
fi
if [ `uname -s` = "Linux" ];
then
  tarprog="tar"
fi
if [ `uname -s` = "Darwin" ];
then
  tarprog="tar"
fi

gzipflag=''
BZ=false
red=false
cat_flags='-f'

while getopts ":bfhpruv" optin; do
   case $optin in
      b  ) BZ=true
           ;;
      f  ) gzipflag="-f"
           ;;
      h | --help )
           echo
           echo "-----------------------------------------------"
           echo "               CLEANFILES HELP"
           echo
           echo "This script deals with the reduction of files in the calculation "
           echo "directory after a run. You must at least enter the star name."
           echo
           echo "USAGE"
           echo "       cleanfiles [-b -f -p -r -u] StarName"
           echo
           echo "OPTIONS"
           echo "       -b: compute burning zones"
           echo "       -f: forced gzipping"
           echo "       -p: keeps preMS lines in reduced file"
           echo "       -r: creates reduced file star.dat"
           echo "       -u: uncleanfiles version"
           echo "-----------------------------------------------"
           echo
           exit 0
           ;;
      p  ) preMS=true
           cat_flags=$cat_flags"p"
           ;;
      r  ) red=true
           cat_flags=$cat_flags"r"
           echo "cat_flags: "$cat_flags
           ;;
      u  ) unclean=true
           ;;
      v | --version )
           echo
           echo "-----------------------------------------------"
           echo "        Script CLEANFILES 3.2 (May 2016)"
           echo
           echo "based on the initial script written by Raphael Hirschi"
           echo "modified and extended by Cyril Georgy & Sylvia Ekstrom"
           echo
           echo "Complaints or suggestions to be adressed to"
           echo "           sylvia.ekstrom@unige.ch"
           echo "-----------------------------------------------"
           echo
           exit 0
           ;;
      \? ) echo "usage : cleanfiles [-b -f -p -r -u] StarName"
           echo
           echo "for a help on this script, type:"
           echo "        cleanfiles -h"
           exit 0
   esac
done

# Decalage des arguments :
shift $(($OPTIND - 1))

if [ -z "$1" ]
then
   echo
   echo "You must at least enter the star name."
   echo
   echo "usage : cleanfiles [-b -f -p -r -u] P007z14S4"
   echo
   echo "for a help on this script:"
   echo "        cleanfiles -h"
   exit
fi

echo "file names begin with" $1

if [ $unclean ];
then
  ls $1.w* *files.tgz
  chmod 644 $1.w* afiles.tgz gfiles.tgz
  echo "untarring afiles ..."
  $tarprog xfz afiles.tgz
  ls $1.a*
  echo "untarring gfiles ..."
  $tarprog xfz gfiles.tgz
  ls $1.g*
  echo "untarring StrucData files ..."
  $tarprog xfz $1_StrucData.tgz
  ls $1_StrucData_*
  if [ -s xfiles.tgz ];
  then
    echo "untarring xfiles ..."
    $tarprog xfz xfiles.tgz
    ls $1.x*
  fi
  if [ -s yfiles.tgz ];
  then
    echo "untarring yfiles ..."
    $tarprog xfz yfiles.tgz
    ls $1.y*
  fi
  if [ -s zfiles.tgz ];
  then
    echo "untarring zfiles ..."
    $tarprog xfz zfiles.tgz
    ls $1.z*
  fi

  echo "unzipping .ws --> "$1".s0000001"
  gunzip $1.ws.gz
  mv $1.ws $1.s0000001

  echo "files "$1".w[a,g] are moved to P"$1".w[a,g]"
  mv $1.wa P$1.wa
  mv $1.wg P$1.wg

  echo "removing tgz files ..."
  rm *.tgz
else
  echo "other files:"
  ls -dp [^P]*

  if [ -s .PlotData_$1 ] && [ $BZ != "true" ];
  then
    BZ=true
  fi

  if [ $BZ ];
  then
    echo "finding burning zones: now processing..."
    $ScriptBZ
  fi

  echo "finding and removing empty files: now processing..."
  find . -size 0 -exec rm {} \;

  echo "cat "$1".a0* > "$1".wa"
  find . -name $1".a0*" -print0 | sort -z | xargs -0 cat > $1.wa


  echo "tar&zipping .a files: tar -cvzf afiles.tar.gz "$1".a0* now processing..."

  find . -name $1".a0*" -print | $tarprog cfz afiles.tgz --files-from -

  echo "looking for ZAMS model..."
  testzams=$(ls $1.s0000[0-8]*)
  if [ -n "$testzams" ];
  then
    zamsmod=$(grep -l 'isol=0' $1.s0000[0-8]* |head -1)
	echo "ZAMS model: " $zamsmod
    modanf=$(egrep modanf $zamsmod | tr -d " ")
    modanf=${modanf##*modanf=}
    modanf=$(printf '%05d' $modanf)
    bzams=$1".b"${modanf}
    echo "zipping ZAMS .b files: $bzams"
    gzip $gzipflag $bzams
  fi

  blast=$(find . -name $1".b[0-9]*" -print | sort | tail -2)
  echo "last .b: " $blast
  #blast=$(ls $1".b[0-9]*" | tail -2)
  echo "zipping last .b files: $blast"
  gzip $gzipflag $blast

  echo "removing unnecessary .b files: rm $1.b*[^0,z] now processing..."
  find . -name $1".b[0-9]*[^0,l,n,z]" -exec rm {} \;

  echo "zipping .b files: gzip $1.b*[^z] now processing..."
  gzip $gzipflag $1.b[0-9]*[^z]

  echo "cat $1.g0* > $1.wg"
  echo "catag launched with flags "$cat_flags
  if [ $red ];
  then
    echo "making reduced file $1.dat"
  fi
  $catag $cat_flags $1

  echo "tar&zipping .g files: tar cvzf gfiles.tar.gz "$1".g0* now processing..."
  find . -name $1".g0*" -print | $tarprog cfz gfiles.tgz --files-from -

  echo "removing unnecessary .l files: rm "$1".l*[^00]1.gz now processing..."
  find . -name $1".l*[^0]1" -exec rm {} \;
  find . -name $1".l*[^0]01" -exec rm {} \;
  gzip $gzipflag $1.l*[^z]

  echo "cat "$1".s0* and zipping " $1".ws"
  find . -name $1".s0*" -print0 | sort -z | xargs -0 cat > $1.ws
  gzip $1.ws

  echo "zipping .v files: gzip "$1".v0* now processing..."
  find . -name $1".v0*" -exec gzip $gzipflag {} \;

  testx=$(find . -name $1".x0*" -print | head -1)
  if [ -n "$testx" ];
  then
    echo "tar&zipping .x files: tar -cvzf xfiles.tar.gz "$1".x0* now processing..."
    find . -name $1".x0*[^z]" -print | $tarprog cfz xfiles.tgz --files-from -
  fi

  testy=$(find . -name $1".y0*" -print | head -1)
  if [ -n "$testy" ];
  then
    echo "tar&zipping .y files: tar -cvzf yfiles.tar.gz "$1".y0* now processing..."
    find . -name $1".y0*[^z]" -print | $tarprog cvfz yfiles.tgz --files-from -
  fi

  testz=$(find . -name $1".z0*" -print | head -1)
  if [ -n "$testz" ];
  then
    if [ -s $testz ];
    then
      echo "tar&zipping .z files: tar -cvzf zfiles.tar.gz "$1".z0* now processing..."
      find . -name $1".z0*" -print | $tarprog cfz zfiles.tgz --files-from -
      find . -name $1".z0*" -exec rm {} \;
    fi
  fi

  teststruc=$(find . -name $1"_StrucData_0*" -print | head -1)
  if [ -n "$teststruc" ];
  then
    if [ -s $teststruc ];
    then
      echo "tar&zipping entire structure files: tar -cfz "$1"_StrucData.tgz now processing..."
      find . -name $1"_StrucData_0*" -print | $tarprog cfz $1"_StrucData.tgz" --files-from -
    fi
  fi

  chmod 444 $1.w* afiles.tgz gfiles.tgz

  echo "Removing files..."
  if [ -s afiles.tgz ];
  then
    find . -name $1".a0*" -exec rm {} \;
  fi
  if [ -s gfiles.tgz ];
  then
    find . -name $1".g0*" -exec rm {} \;
  fi
  if [ -s $1.ws.gz ];
  then
    find . -name $1".s0*" -exec rm {} \;
  fi
  if [ -s xfiles.tgz ];
  then
    find . -name $1".x0*[^z]" -exec rm {} \;
  fi
  if [ -s yfiles.tgz ];
  then
    find . -name $1".y0*[^z]" -exec rm {} \;
  fi
  if [ -s $1"_StrucData.tgz" ];
  then
    find . -name $1"_StrucData_[0-9]*.dat" -exec rm {} \;
  fi

fi

echo "end"
echo -e "\a"