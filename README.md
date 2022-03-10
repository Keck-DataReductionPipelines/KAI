# KAI

The Keck AO Imaging (KAI) data reduction pipeline is a tool to reduce imaging observations taken with the [NIRC2](https://www2.keck.hawaii.edu/inst/nirc2/) and [OSIRIS](https://www2.keck.hawaii.edu/inst/osiris/) near-infrared imagers at the W. M. Keck Observatory.

Installation instructions are below, while more detailed API documentation is available [here](https://keck-datareductionpipelines.github.io/KAI/).

[![Powered by Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org)

## Installation

1. Create a separate [conda](https://docs.conda.io/en/latest/miniconda.html) environment to run KAI. The pipeline uses IRAF/PyRAF, and we recommend using the [`environment_iraf27.yml`](environment_iraf27.yml) file in this repository to [create a conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) with the necessary dependencies correctly installed. The environment file will create a new conda environment called `iraf27`, and must be activated before running KAI using

   ```bash
   conda activate iraf27
   ```

   **Note**: KAI's IRAF / PyRAF dependency currently requires Python 2.7 and operating systems that support 32-bit software.

2. Clone this git repository. For example:

   ```bash
   cd ~/software/KAI
   git clone git@github.com:Keck-DataReductionPipelines/KAI.git
   ```

3. Install KAI by going to your cloned repository and running the `setup.py` script. For example:

   ```bash
   conda activate iraf27
   cd ~/software/KAI/
   python setup.py install
   ```

4. Test your installation by importing KAI in python. For example:

   ```python
   from kai.reduce import data
   ```

After installation, try running the [reduction tutorial](kai/TheReductionGuide.ipynb) to get up to speed with KAI.

## Example Reduction Template Scripts

The reduction template scripts included in this repository provide a complete run-through of the reduction procedure for imaging data: creating darks and flats, reducing skies, cleaning science images, and combining multiple clean science images into a combo science image.
* [Reduction template script for NIRC2 imaging data](kai/reduce/reduce_template.py)
* [Reduction template script for OSIRIS imaging data](kai/reduce/reduce_template_osiris.py)

## Example Dataset

An example dataset with scripts can be found at this [Google Drive link](https://drive.google.com/drive/folders/1FpTN3wiG4U826H328JIJcPLbScNCTRQW?usp=sharing). This is a great place to start to test the pipeline.



Contributing
------------

We love contributions! KAI is open source, built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not ready to be an open source contributor; that your skills aren't nearly good enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at all, you can contribute code to open source. Contributing to open source projects is a fantastic way to advance one's coding skills. Writing perfect code isn't the measure of a good developer (that would disqualify all of us!); it's trying to create something, making mistakes, and learning from those mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can help out by writing documentation, tests, or even giving feedback about the project (and yes - that includes giving feedback about the contribution process). Some of these contributions may be the most valuable to the project as a whole, because you're coming to the project with fresh eyes, so you can see the errors and assumptions that seasoned contributors have glossed over.

Note: This disclaimer was originally written by [Adrienne Lowe](https://github.com/adriennefriend) for a [PyCon talk](https://www.youtube.com/watch?v=6Uj746j9Heo), and was adapted by KAI based on its use in the README file for the [MetPy project](https://github.com/Unidata/MetPy).

License
-------

This project is Copyright (c) J.R. Lu, A. K. Gautam, T. Do and licensed under the terms of the BSD 3-Clause license. This package is based upon the [Astropy package template](https://github.com/astropy/package-template) which is licensed under the BSD 3-clause license. See the licenses folder for more information.
