ionize
=====
[![Code Climate](https://codeclimate.com/github/lewisamarshall/ionize/badges/gpa.svg)](https://codeclimate.com/github/lewisamarshall/ionize) [![Build Status](https://travis-ci.org/lewisamarshall/ionize.svg?branch=master)](https://travis-ci.org/lewisamarshall/ionize)  [![Coverage Status](https://coveralls.io/repos/lewisamarshall/ionize/badge.svg?branch=master&service=github)](https://coveralls.io/github/lewisamarshall/ionize?branch=master) [![Docs Status](https://readthedocs.org/projects/ionize/badge/?version=latest)](https://ionize.readthedocs.org)

**ionize** is a Python package for calculating the properties of ions and
electrolytes. Originally designed with electrophoresis in mind,
these calculations can help any time the pH or electrical properties of an
aqueous solution impact system performance.

The goal of **ionize** is to make full simulations of electrolyte properties
easy, ubiquitous, and accurate. The **ionize** models consider the impact of pH,
ionic strength, temperature, and the interactions between different ions.
Considering all of these impacts makes ionize accurate over the widest
available range of operating conditions. In addition, **ionize** includes
warnings when experimental conditions are outside the range of model
assumptions.

Installation
------------
One-line install using [pip](https://pypi.python.org/pypi/pip):

    pip install ionize

Tutorial
--------
Want to use **ionize**? Read the [tutorial][tutorial].

Examples
--------
Want to see some examples of **ionize** in action? Take a look at the
[examples][examples].

[peakmaster]: http://web.natur.cuni.cz/gas/ "Peakmaster"
[Spresso]: http://stanfordspresso.blogspot.com/ "Spresso"
[STEEP]: http://microfluidics.stanford.edu/download/ "STEEP"
[tutorial]: ./tutorial.ipynb  "ionize Tutorial"
[examples]: ./examples.ipynb  "ionize Examples"
