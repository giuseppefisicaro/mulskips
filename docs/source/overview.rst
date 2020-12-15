The MulSKIPS code
=================

`MulSKIPS` is a Kinetic Monte Carlo super-Lattice code, designed to study with an atomic resolution
the growth kinetics of elements, alloys and compounds characterized by the sp3 bond symmetry.
Formalization and implementation details of the code are discussed in Ref. [1_].
The code is open source and it is distributed according to the GNU public license.
MulSKIPS is available on GitLab_, from where it can be downloaded as a tar file.

.. _GitLab: https://github.com/giuseppefisicaro/mulskips

Deposition and evaporation of the substrate atoms are the active Monte Carlo events,
driving the stochastic evolution. In MulSKIPS, a dense super-lattice correctly
accommodates the original lattice of the ideal crystal along with a large class
of defective configurations [2_]. This feature makes the code unique in the range
of lattice Kinetic Monte Carlo codes currently available for sp3 materials.
Indeed, the code is able to simulate the evolution of both point and extended defects,
like stacking faults of different symmetries, antiphase boundaries and grain boundaries.
Moreover, MulSKIPS can simulate the morphological evolution during the growth process,
e.g. the epitaxial growth or etching of flat, structured, or patterned substrates,
as well as nanoparticles of various shapes.
In the case of surfaces, periodic boundary conditions are applied in the planes
orthogonal to the growth direction.

MulSKIPS functionalities
========================

Following we list all MulSKIPS functionalities:

* Epitaxial growth of nanoparticles, flat and patterned substrates, micro- nano-electronic devices by:

  * physical vapor deposition [released];
  * chemical vapor deposition [under development];

* Melting laser annealing processes [under development].


References:

.. [1] A. L. Magna et al., Simulation of the Growth Kinetics in Group IV Compound Semiconductors, physica status solidi (a) vol. 216, no. 10, p. 1800597, 2019, doi: 10.1002/pssa.201800597

.. [2] G. Fisicaro et al., Genesis and Evolution of Extended Defects: The Rrole of Evolving Interface Instabilities in Cubic SiC, Applied Physics Reviews vol. 7, no. 2, p. 021402, Apr. 2020, doi: 10.1063/1.5132300
