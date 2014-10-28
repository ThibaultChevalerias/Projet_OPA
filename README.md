Projet_OPA
==========

(to be completed)

This code is written in the framework of an engineering school project of two third-year students of École Centrale Paris. 

It aims at modelling an ensemble of spins in a material as a cubic or rectangular cuboid 3-dimensional set of spins.
We currently use the Ising model, thus the spins have values + 1 or -1. A Heisenberg model might be used later.
The values (+ 1 or - 1) of each spin is stored in a Tensor, which is an object created to emulate a tensor of order 3 using a C++ vector variable.

A Monte-Carlo simulation is implemented, in which we choose to flip a random spin or not at each step of modelisation, given the change in energy it induces. Several Monte-Carlo steps (usually arond 1 million) are performed to have a good equilibrium.
This simulation is then runned starting from a given temperature, and decreasing it by little steps. We get an evolution of interesting physical variables as a function of temperature, and we aim at finding a transition temperature between a disordered and an ordered phase, and identify these phases.
This study is runned for several values of the parameter J2 (exchange constant), in order to find the full phase diagram in the (-J2/J1; kT/J0) space.
