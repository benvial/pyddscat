*********************
:mod:`pyddscat.targets`
*********************

.. automodule:: pyddscat.targets

Inheritance Diagram
===================

.. inheritance-diagram:: pyddscat.targets
    :parts: 1

Class Summary
=============
.. currentmodule:: pyddscat.targets

Base Classes
------------
Otherclasses derive from these. Generally not intended to be used on their own.

.. autosummary::
    Target
    Target_Builtin
    Periodic

Isolated Builtin Targets
------------------------
Classes for defining isolated finite targets based on target definitions
built into DDSCAT.

.. autosummary::
    Target_Builtin
    RCTGLPRSM
    CYLNDRCAP
    ELLIPSOID
    CYLINDER
    Sphere
    Cube

Arbitrarily Shaped Targets
--------------------------
Classes for defining isolated finite targets of arbitrary shape.

.. autosummary::
    FROM_FILE
    FRMFILPBC
    Iso_FROM_FILE
    Ellipsoid_FF
    Helix
    SpheresHelix
    Conical_Helix

Periodic Targets
----------------
Classes for defining targets with 1D and 2D semi-infinite periodicity.

.. autosummary::
    Periodic
    FRMFILPBC
    RCTGL_PBC
    CYLNDRPBC


Class Definitions
=================

.. automodule:: pyddscat.targets
   :members:
   :undoc-members:
   :show-inheritance:
   :member-order: bysource