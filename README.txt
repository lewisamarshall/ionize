======
ionize
======

A unified Python package for calculating buffer properties.

**ionize** calculates the properties of individual ionic species in
aqueous solution, as well as aqueous solutions containing arbitrary
sets of ions.

Components
=========

**ionize** is composed of three main components:

The Ion Class
-------------
The core of **ionize** is the **Ion** class, which  represents a single ionic
species. An ion contains a name, a set of ionization states, and an optional
temperature parameter. Each ionization state contains a charge, a pKa, and
an absolute mobility. An ionization  state may also include values for
&Delta;H and &Delta;Cp of ionization to improve the accuracy of temperature
correction.

The Solution Class
------------------
The **Solution** class is used to represent an aqueous solution containing any
number of ionic species. A **Solution** contains a list of **Ion** objects, and
a second list containing the concentrations of each species. **Solution** can
also take an optional temperature parameter. **Solution** solves for pH,
iteratively accounting for the ionic strength. The ionic strength and pH are
used to calculate the properties of the **Ions**, and bulk properties of the
solution such as conductivity.

The ionize Database
-------------------
**ionize** also contains its own database, containing the combined entries of
the Spresso and STEEP databases. This database can be accessed through the
**load_ion**() and **search_ion**() functions.
