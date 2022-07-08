<img src="https://raw.githubusercontent.com/ChaosGroup/microsurface-transformations/master/teaser.jpg" alt="Microsurface transformations teaser">

# Microsurface Transformations

Code accompanying the paper ["Microsurface Transformations"](https://docs.chaos.com/display/RESEARCH/Microsurface+Transformations) by Asen Atanasov, Vladimir Koylazov, Rossen Dimov and Alexander Wilkie from EGSR 2022.

## Derivation

In the folder *derivation* we provide the [Mathematica](https://www.wolfram.com/mathematica/) notebook for the derivation of the Jacobian that ensures normalization of linearly transformed microfacet distributions.

## Numerical Validation

In *microsurface_transformations.h* we provide the implementations of a variety of microsurfaces: [GTR](https://disneyanimation.com/publications/physically-based-shading-at-disney/), [Anisotropic GGX and Beckmann](https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.html), [STD](https://hal.archives-ouvertes.fr/hal-01535614), [Phong](https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.html), [Sheen](http://www.aconty.com/pdf/s2017_pbs_imageworks_sheen.pdf) and [Discrete GGX](https://www.cs.cornell.edu/projects/stochastic-sg14/). In *test.cpp* we instantiate all of these microsurfaces with random parameters and verify that they fullfill normalization and shadowing constraints. Then these instances are transformed using random tangential transformations and the constraints are verified once again.