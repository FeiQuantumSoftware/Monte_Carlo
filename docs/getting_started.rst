Getting Started
===============

This page details how to get started with monte_carlo. monte_carlo is a package 
which was developed for the MolSSI Bes t Practices workshop.

Installation 
------------
To install molecool, you will need an environment with the following packages:

* Python 3.9
* NumPy 

Once you have these packages installed, you can install molecool in the same environment using 
:: 

    pip install -e .


Usage 
------
Once installed, you can use the package. This example shows how to initialize a spinlist.
::

    import monte_carlo

    example_spinlist = molecool.SpinConfig(8)
    example_spinlist.init_decimal_input(10) = [0, 0, 0, 0, 1, 0, 1, 0]

