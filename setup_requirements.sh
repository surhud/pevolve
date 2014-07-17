#!/bin/sh

rm -rf build/ install/ src/*.py src/*wrap*

# First test if there is a gsl installation
gslex=`gsl-config --libs` 

# Exit shell script on error from here on
set -e

# Check for gsl installation
if [ -z "$gslex" ]
then
  # Install GSL
  echo "++++++++++++++++++++++++++++++++++++++++++++++"
  echo "GSL is not installed. Trying to install it:"
  echo "++++++++++++++++++++++++++++++++++++++++++++++"
  wget https://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
  tar -zxvf gsl-latest.tar.gz
  cd gsl-*
  ./configure --prefix=`pwd`/../install
  make -j8 
  make install
  cd ../install
  idir=`pwd`
  # Setup flags for GSL
  incflag="-I$idir/include"
  libflag="-L$idir/lib -lgsl -lgslcblas"
  lddir=$idir/lib
  cd ../
  echo "# Please add $lddir to your path:\n LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$lddir" > instructions.txt
  echo "# Please add $idir/bin to your path:\n PATH=\$PATH:$idir/bin" >> instructions.txt
else
  echo "GSL is installed"
  # Setup flags for GSL
  libflag=`gsl-config --libs`
  incflag=`gsl-config --cflags` 
fi

echo "LDFLAGS is: " $libflag
echo "CPPFLAGS is: " $incflag

set +e

# Now test for swig

swigex=`swig -help`

set -e
if [ -z "$swigex" ]
then
  # Install swig
  echo "++++++++++++++++++++++++++++++++++++++++++++++"
  echo "swig is not installed. Trying to install it:"
  echo "++++++++++++++++++++++++++++++++++++++++++++++"
  wget http://prdownloads.sourceforge.net/swig/swig-3.0.2.tar.gz
  tar -zxvf swig-*.tar.gz
  cd swig-*
  ./configure --prefix=`pwd`/../install
  make -j8
  make install
  cd ../install
  idir=`pwd`
  echo "# Please add $idir/bin to your path: \n PATH=\$PATH:$idir/bin" >> instructions.txt
else
  echo "swig is installed"
fi

set +e

# Now mandc code
mandcexec=`PATH=$PATH:. mandc.x`
# Check for gsl installation
if [ -z "$mandcexec" ]
then
  # Install mandc.x
  echo "++++++++++++++++++++++++++++++++++++++++++++++++"
  echo "mandc is not installed. Trying to install it:"
  echo "++++++++++++++++++++++++++++++++++++++++++++++++"
  set -e
  wget 202.127.29.4/dhzhao/mandc_code/mandc-1.03main.tar.gz
  tar -zxvf mandc-1.03main.tar.gz
  cd mandc-1.03main
  set +e
  
  success=1
  for comp in ifort f77 gfortran f95 g77
      do
          hash $comp 2>&-
          success=$?
          if [[ $success -ne 0 ]]; then
              echo "I could not find $comp";
              continue;
          else
              echo "Found $comp"
              echo "Executing patch makefile ../patches/patch.$comp"
              patch makefile ../patches/patch.$comp
              break;
          fi
      done
  
  if [[ $success -eq 1 ]]; then
      echo "Sorry could not find any of the common fortran compilers:"
      echo "f77 g77 gfortran f95 or ifort"
      echo "If you have installed these please add the path to the"
      echo "compilers to your system variable PATH by executing"
      echo "export PATH=\$PATH:path_to_compiler"
      echo "and try again. Exiting now!"
  else
      echo Found compiler $comp
      set -e
      make
      cd ../install/bin
      ln -sf ../../mandc-1.03main/mandc.x .
      echo "The code mandc.x is ready to use!"
      set +e
  fi
  echo "# Please add $idir/bin to your path: \n PATH=\$PATH:$idir/bin" >> instructions.txt
else
  echo "mandc.x is installed"
fi

# To compile
PATH=$PATH:./install/bin LDFLAGS="$LDFLAGS $libflag" CPPFLAGS="$CPPFLAGS $incflag" python setup.py install --prefix=`pwd`/install
# Ensure all files are installed correctly after recompile
PATH=$PATH:./install/bin LDFLAGS="$LDFLAGS $libflag" CPPFLAGS="$CPPFLAGS $incflag" python setup.py install --record "installlog.txt" --prefix=`pwd`/install

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "# Add the following to your PYTHONPATH variable in your bashrc:"
echo "# Add the following to your bashrc:" >> instructions.txt
cat installlog.txt | grep cosmology.pyc | sed s/cosmology.pyc//
newppath=`cat installlog.txt | grep cosmology.pyc | sed s/cosmology.pyc//`
echo PYTHONPATH=$PYTHONPATH:$newppath >>instructions.txt
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

echo "Trying to install documentation"
sed s:inspath:`cat installlog.txt | grep cosmology.pyc | sed s/cosmology.pyc//  `: doc/source/conf.py.orig> doc/source/conf.py

#"Check for Sphinx"
sphinxex=`sphinx-build --version` 
set -e
if [ -z "$sphinxex" ]
then
  # Install Sphinx
  echo "++++++++++++++++++++++++++++++++++++++++++++++++"
  echo "sphinx is not installed. Trying to install it:"
  echo "No worries if this step does not succeed."
  echo "++++++++++++++++++++++++++++++++++++++++++++++++"
  pip install --install-option="--prefix=`pwd`/install" sphinx
else
  echo "sphinx is installed"
fi

cd doc
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$lddir sphinx-build -b html ./source ./build
