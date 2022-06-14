# mvf

<img width="400" alt="visualisation of a multivector field" src="https://user-images.githubusercontent.com/95714236/173676786-4778f3a2-0cea-42cd-bada-533adc37a394.png">
  
A small project built for DSC291-E00 Topological Data Analysis.  Implements part of "Persistence of the Conley Index in Combinatorial Dynamical Systems" by Tamal K. Dey, Marian Mrozek, and Ryan Slechta (https://arxiv.org/abs/2003.05579).

A multivector field is a discrete representation of a dynamical system that represents it as a partitioning of an input mesh (or simplicial complex) into "multivectors".  Using these multivectors, one can describe the topological properties of the dynamical system using the theory of simplicial homology.

In "Persistence of the Conley Index in Combinatorial Dynamical Systems", Dey et al. use the framework of persistent homology to describe changes in the topology of a dynamical system as the system undergoes changes over time or with respect to some parameter.  Specifically, the work captures changes in the _Conley index_, which represents the dynamical behavior of the invariant sets of a system using relative homology.

Due to time constraints, this repository is incomplete and does not work for vector fields which change over time, but still manages to correctly compute the Conley index for fixed, noise-free vector fields.

Repository also contains code for creating and visualising meshes in the plane and multivector fields on those meshes.  Implements the halfedge mesh data structure (common in graphics and elsewhere), and can instantiate from .obj-like representations (i.e. lists of vertices and faces).

Please read the .pdf in the repository for a more detailed (albeit less formal or detailed as the cited paper) explanation of multivector fields and Conley index persistence, as well as a deeper look at my visualisations and results.

Uses Dionysus 2 for zigzag homology, and NumPy / Matplotlib for all other numerical routines and visualisation.
