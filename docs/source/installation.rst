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

The required dependencies are `networkx <https://networkx.org>`_ and
`pandas <https://pandas.pydata.org>`_.


Optional dependencies
=====================

.. code-block:: bash

   pip install "skm-tools[<optional-dependencies>] @ git+https://github.com/NIB-SI/skm-tools.git"

Where ``<optional-dependencies>`` is a comma-separated list of:

* ``cytoscape``: Cytoscape automation (:mod:`skm_tools.cytoscape_utils`), requires
  `py4cytoscape <https://py4cytoscape.readthedocs.io>`_ and a running
  `Cytoscape <https://cytoscape.org>`_
* ``pdf``: batch export of Cytoscape networks to PDF (:mod:`skm_tools.cytoscape_pdf_utils`),
  requires pypdf, pdfCropMargins and reportlab
* ``tutorials``: to run the :doc:`tutorials <tutorials/tutorial-pss>` (matplotlib,
  py4cytoscape and a Jupyter kernel)
* ``docs``: to build this documentation
* ``test``: to run the tests (pytest)

For example:

.. code-block:: bash

   pip install "skm-tools[cytoscape,pdf] @ git+https://github.com/NIB-SI/skm-tools.git"

Graphviz layouts additionally need `pygraphviz <https://pygraphviz.github.io>`_.


Tests
=====

.. code-block:: bash

   git clone https://github.com/NIB-SI/skm-tools.git
   cd skm-tools
   pip install ".[test]"
   pytest

See https://github.com/NIB-SI/skm-tools/tree/main/tests. The Cytoscape tests mock
py4cytoscape, so Cytoscape doesn't need to run.


Building the documentation
==========================

.. code-block:: bash

   pip install ".[docs]"
   cd docs
   make clean
   make html

The pages are written to ``docs/build/html``.
