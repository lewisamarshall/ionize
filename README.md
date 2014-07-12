AspPy
=====

A python port of Asp, the Accessible Solution Parameters electrolyte simulation tool.

Ion
----
The Ion class represents represents ionic species. Initialize an ion with
Ion('name', [z1,...], [pKa1, ...], [µ1, ...]). Ions can also be initialized
using the load_ion function.

Solution
--------
Class Solution represents an aqueous solution containing one or more ions. This
class automatically solves for the pH of the solution , iteratively accounting
for the ionic strength. The ionic strength and pH are used to calculate the
properties of the ions, and bulk properties of the solution such as
conductivity.

load_ion
--------
Initialize an ion using load_ion('name'). The database
is stored in ions_shelve.db, based on the Peakmaster database, with additions from literature.
