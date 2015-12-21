"""
    ionize
    ~~~~~~

    Ionize is a Python package for calculating the properties of ions and
    electrolytes. Originally designed with electrophoresis in mind,
    these calculations can help any time the pH or electrical properties of an
    aqueous solution impact system performance.

    :copyright: (c) 2015 by Lewis A. Marshall.

"""
from .__version__ import __version__
from .Solvent import Aqueous
from .Ion import Ion
from .PolyIon import NucleicAcid, Peptide
from .IonComplex import IonComplex, Protein
from .Solution import Solution
from .deserialize import deserialize
from .Database import Database

# TODO: include pitzer model for high ionic strength NaCl activity
