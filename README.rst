KAI
---

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge


License
-------

This project is Copyright (c) J.R. Lu, A. Gautam, T. Do and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause license. See the licenses folder for
more information.

Example Dataset
---------------
An example dataset with scripts can be found at this 
`Google Drive link <https://drive.google.com/drive/folders/1FpTN3wiG4U826H328JIJcPLbScNCTRQW?usp=sharing>`_. 
This is a great place to start to test the pipeline.

Installation Guide
------------------

Make sure the python paths point to the appropriate code directories. For example, this can be done by using setup.py:

``python setup.py install``

Or manually putting in the path to the repository in your ``.bash_profile``, such as:

``export PYTHONPATH=$PYTHONPATH:path/to/KAI``

It is best to create a separate conda environment. Currently, the pipeline uses IRAF/PyRAF, and it is helpful to create an environment as instructed here on the `STScI astroconda website <https://astroconda.readthedocs.io/en/latest/installation.html>`_. Note that this pipelines requires operating systems that support 32-bit software.

Contributing
------------

We love contributions! KAI is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
KAI based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.
