.. KAI documentation master file, created by
   sphinx-quickstart on Thu Feb 17 12:07:34 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Keck AO Imaging (KAI) data reduction pipeline
=================================================

The Keck AO Imaging (KAI) data reduction pipeline is a tool to reduce imaging observations taken with the `NIRC2 <https://www2.keck.hawaii.edu/inst/nirc2/>`_ and `OSIRIS <https://www2.keck.hawaii.edu/inst/osiris/>`_ near-infrared imagers at the W. M. Keck Observatory.

Download
--------

The latest release can be downloaded from the github repository `here
<https://github.com/Keck-DataReductionPipelines/KAI>`_.

Installation
------------

1. Create a separate `conda <https://docs.conda.io/en/latest/miniconda.html>`_ environment to run KAI. The pipeline uses IRAF/PyRAF, and we recommend using the ``environment_iraf27.yml`` file in this repository (available `here <https://github.com/Keck-DataReductionPipelines/KAI/blob/dev/environment_iraf27.yml>`_) to `create a conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_ with the necessary dependencies correctly installed::
    
    conda env create -f environment_iraf27.yml


The environment file will create a new conda environment called ``iraf27``, and must be activated before running KAI using::

    conda activate iraf27
   
**Note**: KAI's IRAF / PyRAF dependency currently requires Python 2.7 and operating systems that support 32-bit software.

2. Clone this git repository. For example::

    cd ~/software/KAI
    git clone git@github.com:Keck-DataReductionPipelines/KAI.git
   

3. Install KAI by going to your cloned repository and running the `setup.py` script. For example::

    conda activate iraf27
    cd ~/software/KAI/
    python setup.py install

4. Test your installation by importing KAI in python. For example::

    from kai.reduce import data

After installation, try running the `reduction tutorial <https://github.com/Keck-DataReductionPipelines/KAI/blob/dev/kai/TheReductionGuide.ipynb>`_ to get up to speed with KAI.

Example Reduction Template Scripts
----------------------------------

The reduction template scripts included in this repository provide a complete run-through of the reduction procedure for imaging data: creating darks and flats, reducing skies, cleaning science images, and combining multiple clean science images into a combo science image.

* `Reduction template script for NIRC2 imaging data <https://github.com/Keck-DataReductionPipelines/KAI/blob/dev/kai/reduce/reduce_template.py>`_
* `Reduction template script for OSIRIS imaging data <https://github.com/Keck-DataReductionPipelines/KAI/blob/dev/kai/reduce/reduce_template_osiris.py>`_

Example Dataset
---------------

An example dataset with scripts can be found at this `Google Drive link <https://drive.google.com/drive/folders/1FpTN3wiG4U826H328JIJcPLbScNCTRQW?usp=sharing>`_. This is a great place to start to test the pipeline.


Index
-----

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   autoapi/index


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
