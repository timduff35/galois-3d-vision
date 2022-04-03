# galois-3d-vision
Code and examples accompanying the paper "[Galois/monodromy groups for decomposing minimal problems in 3D reconstruction](https://arxiv.org/pdf/2105.04460.pdf)."

# Summary of files and subdirectories

| File/Directory | Description |
| ---------------| ----------- |
| ``groups-relative-pose/``| Contains [GAP](https://www.gap-system.org/) code for analyzing the 4 Galois/monodromy groups appearing in Result 4.1|
| ``P3P.gp``| [GAP](https://www.gap-system.org/) code output by ``example-tour.m2`` giving the Galois/monodromy group of Grunert's equations for P3P|
| ``groups-absolute-pose/``| contains [GAP](https://www.gap-system.org/) code for analyzing the 4 Galois/monodromy groups appearing in Result 3.1|
| ``common.m2`` | Various utility functions in [Macaulay2](http://www2.macaulay2.com/Macaulay2/)|
|``example-21.m2``| [GAP](https://www.gap-system.org/) code output by ``example-tour.m2`` giving the Galois/monodromy group of a three-fold cyclic branched cover|
|``example-tour.m2``| Illustrates examples & computations appearing throughout the paper using [Macaulay2](http://www2.macaulay2.com/Macaulay2/)|
|``gal-5pt-20.m2``| [GAP](https://www.gap-system.org/) code output by ``example-tour.m2`` giving the Galois/monodromy group of the five-point problem with translations and depths in a common projective space|
|``gal-5pt-40.m2``| [GAP](https://www.gap-system.org/) code output by ``example-tour.m2`` giving the Galois/monodromy group of the five-point problem with unit-length translation|
|``plmp-problem-builder.m2``| [Macaulay2](http://www2.macaulay2.com/Macaulay2/) code called by ``result-41.m2`` for particular problems. Much of this code is derived from a [previous project](https://github.com/timduff35/PLMP).|
|``result-31.m2``| [Macaulay2](http://www2.macaulay2.com/Macaulay2/) code computing permutations in Galois/monodromy groups appearing in Result 3.1|
|``result-41.m2``| [Macaulay2](http://www2.macaulay2.com/Macaulay2/) code computing permutations in Galois/monodromy groups appearing in Result 4.1|
| ``P3P.gp``| [GAP](https://www.gap-system.org/) code output by ``example-tour.m2`` giving the Galois/monodromy group of the "sparsification" of Grunert's equations (ie. ignoring the special structure in their coefficients.)|
