.. KAI documentation master file, created by
   sphinx-quickstart on Thu Feb 17 12:07:34 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Keck AO Imaging (KAI) data reduction pipeline
=================================================

The Keck AO Imaging (KAI) data reduction pipeline is a tool to reduce imaging observations taken with the `NIRC2 <https://github.com/Keck-DataReductionPipelines/KAI/blob/dev/kai/reduce/TheReductionGuide.ipynb>`_ and `OSIRIS <https://www2.keck.hawaii.edu/inst/osiris/>`_ near-infrared imagers at the W. M. Keck Observatory.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   autoapi/index

Download
--------

The latest release can be downloaded from the github repository `here
<https://github.com/Keck-DataReductionPipelines/KAI>`_.

Installation
------------

1. Create a separate `conda <https://docs.conda.io/en/latest/miniconda.html>`_ environment to run KAI. The pipeline uses IRAF/PyRAF, and we recommend using the ```environment_iraf27.yml`` <environment_iraf27.yml>`_ file in this repository to `create a conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_ with the dependencies correctly installed. The environment file will create a new conda environment called ``iraf27``, and must be activated before running KAI using::

    conda activate iraf27
   
**Note**: The IRAF / PyRAF dependency currently requires Python 2.7 and operating systems that support 32-bit software.

2. Clone this git repository. For example::

    cd ~/software/KAI
    git clone git@github.com:Keck-DataReductionPipelines/KAI.git
   

3. Install KAI by going to your cloned repository and running the `setup.py` script. For example::

    conda activate iraf27
    cd ~/software/KAI/
    python setup.py install

4. Test your installation by importing KAI in python. For example::

    from kai.reduce import data

After installation, try running the `reduction tutorial <kai/TheReductionGuide.ipynb>`_ to get up to speed with KAI.


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
