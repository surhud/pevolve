Installation
============

The pseudo-evolution package depends upon the following
    1. `swig <http://www.swig.org>`_
    2. `GNU scientific library <http://www.gnu.org/software/gsl>`_
    3. `The public mandc code by Donghai Zhao <http://202.127.29.4/dhzhao/mandc.html>`_

Optional
    - `Sphinx for documentation <http://sphinx-doc.org>`_

The installer can install item number 2 and 3 if the computer has access to the internet.

Step by step instructions
=========================

Download source code
    - `Source code for pevolve <http://github.com/surhud/pevolve>`_

.. sourcecode:: bash

    $ git clone http://github.com/surhud/pevolve
    $ cd pevolve
    $ sh setup_requirements.sh

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
