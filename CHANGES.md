Ionize Changelog
================

0.8.0
-----
Fixes to diffusivity for partially ionized species.
Using load_ion with an invalid key will now raise an error instead of returning None.

0.7.0
-----
Added get_concentration as a method for Solutions.

0.6.0
-----
Added diffusivity as a method for Ions.

0.5.3
-----
Fixed a bug in which taps had an extra value for dCp.

0.5.2
-----
Fixed a problem in which pKa adjustment for temperature was not being performed.

0.5.1
-----
Corrected a problem in which the nucleic acid mobility was inaccurate.

0.5.0
-----
Added Debye length as a property of Solutions.
Added Solution.set_T().

0.4.0
-----
Added support to initialize solutions with ion names instead of ions.
