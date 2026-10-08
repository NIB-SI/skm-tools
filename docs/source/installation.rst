============
Installation
============

Install from GitHub
===================

.. code-block:: bash

   pip install "skm-tools @ git+https://github.com/NIB-SI/skm-tools.git"

Or download and install:

.. code-block:: bash

   git clone https://github.com/NIB-SI/skm-tools.git
   cd skm-tools
   pip install .

skm-tools needs Python 3.10 or later. The required dependencies are
`networkx <https://networkx.org>`_ (3 or later) and `pandas <https://pandas.pydata.org>`_
(2 or later).


Optional dependencies
=====================

.. code-block:: bash

   pip install "skm-tools[<optional-dependencies>] @ git+https://github.com/NIB-SI/skm-tools.git"

Where ``<optional-dependencies>`` is a comma-separated list of:

* ``cytoscape``: Cytoscape automation (:mod:`skm_tools.cytoscape_utils`), requires
  `py4cytoscape <https://py4cytoscape.readthedocs.io>`_ and a running
  `Cytoscape <https://cytoscape.org>`_ (3.10 or later); node charts also need the
  `enhancedGraphics <https://apps.cytoscape.org/apps/enhancedgraphics>`_ app
* ``pdf``: cropping and combining PDFs, e.g. exported networks (:mod:`skm_tools.pdf_utils`),
  requires pypdf, pdfCropMargins and reportlab
* ``tutorials``: to run the :doc:`tutorials <tutorials/tutorial-pss>` (matplotlib,
  py4cytoscape and a Jupyter kernel)
* ``docs``: to build this documentation
* ``test``: to run the tests (pytest)

For example:

.. code-block:: bash

   pip install "skm-tools[cytoscape,pdf] @ git+https://github.com/NIB-SI/skm-tools.git"

Graphviz layouts additionally need `pygraphviz <https://pygraphviz.github.io>`_.


.. _logging:

Messages
========

The functions report what they did (e.g. how many nodes a filter removed, a download, or
nodes that weren't found) with Python's standard :mod:`logging`, under the ``skm_tools``
logger. Python shows only warnings by default; to see all messages, e.g. at the top of a
notebook:

.. code-block:: python

   import logging
   logging.basicConfig(level="INFO")

Or only those of skm-tools: ``logging.getLogger("skm_tools").setLevel("INFO")`` (with a
handler, e.g. from ``logging.basicConfig()``).


Tests
=====

.. code-block:: bash

   git clone https://github.com/NIB-SI/skm-tools.git
   cd skm-tools
   pip install ".[test]"
   pytest

See https://github.com/NIB-SI/skm-tools/tree/main/tests. The Cytoscape tests mock
py4cytoscape, so Cytoscape doesn't need to run.

The tutorials double as integration tests: the CI runs them without their Cytoscape cells.
To run them completely, with a running Cytoscape, use ``tutorials/run-tutorials.sh``.


Building the documentation
==========================

.. code-block:: bash

   pip install ".[docs]"
   cd docs
   make clean
   make html

The pages are written to ``docs/build/html``.
