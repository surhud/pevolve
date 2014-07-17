Installation
============

The pseudo-evolution package depends upon the following
    1. `swig <http://www.swig.org>`_
    2. `GNU scientific library <http://www.gnu.org/software/gsl>`_
    3. `The public mandc code by Donghai Zhao <http://202.127.29.4/dhzhao/mandc.html>`_
    4. Any of the following fortran compilers: (gfortran, g77, f77, ifort)

Optional
    - `Sphinx for documentation <http://sphinx-doc.org>`_

The installer can install the first three dependencies if the computer has access to the internet and at least python (and its header files) is installed. On Ubuntu install most of the dependencies using:

.. sourcecode:: bash

    $ sudo apt-get install python python-dev python-sphinx swig libgsl0ldbl libgsl0-dev gfortran

On Fedora linux:

.. sourcecode:: bash

    $ sudo yum install python python-devel python-sphinx swig gsl gcc-gfortran

On a mac, if you have homebrew,

.. sourcecode:: bash

    $ brew install python swig gsl gfortran

Step by step instructions
=========================

Download source code
    - `Source code for pevolve <http://github.com/surhud/pevolve>`_

.. sourcecode:: bash

    $ git clone http://github.com/surhud/pevolve
    $ cd pevolve
    $ chmod a+x setup_requirements.sh
    $ ./setup_requirements.sh

If all goes well, then pevolve should be installed in the install directory
inside pevolve. Now add the directory where pevolve was installed to you
PYTHONPATH variable and the current directory to the PATH variable.

.. sourcecode:: bash

    $ export PYTHONPATH=$PYTHONPATH:`pwd`/install/lib/python2.7/site-packages
    $ export PATH=$PATH:`pwd`

Please modify the python path in the code above depending upon your system
install. Check installlog.txt to see the path to pevolve for example. 

If the installer had to install gsl on your system then you need to add the path
where gsl was installed to your LD_LIBRARY_PATH environment variable. 

.. sourcecode:: bash

    $ # Only if gsl had to be installed by the installer
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/gsl-1.16/install/lib/

The code should now be setup for use! Please see the cosmology documentation for
further examples.
