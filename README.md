Description
--

This repository contains Magma code for reconstructing plane quartic curves from their Dixmier--Ohno invariants and for calculating isomorphisms between plane quartic curves. The latter code is written in terrible style for now, and it only works over the rationals and finite fields. In some non-generic cases, we use an approach due to Andrew Sutherland that involves Groebner bases.

Prerequisites
--

An installation of Magma.

The algorithms use the Magma algorithm `MinimizeReducePlaneQuartic`, due to Elsenhans and Stoll, to simplify any output over the rationals. In order for this to work properly, one needs a bug fix of the file `magma/package/Geometry/SrfDP/`, which can be found in `magma/BugFix.m`. Be warned that this new file is attached by default and supersedes its version in Magma.

Installation
--

You can enable the functionality of this package in Magma by attaching the `quartic/magma/spec` file with `AttachSpec`. To make this independent of the directory in which you find yourself, and to active this on startup by default, you may want to indicate the relative path in your `~/.magmarc` file, by adding the line
```
AttachSpec("~/Programs/quartic/magma/spec");
```

Usage
--

Examples are given in the directory [`examples`](examples) (a full list of intrinsics is [here](intrinsics.md)).

Verbose comments are enabled by
```
SetVerbose("Quartic", n);
```
where `n` is either `1` or `2`. A higher value gives more comments.

Credits
--

This package uses code from the Magma package [`echidna`](http://iml.univ-mrs.fr/~kohel/alg/index.html) by David Kohel for an implementation of the Dixmier--Ohno invariants.

Citing this code
--

Please cite the following preprint if this code has been helpful in your research:

Reynald Lercier, Christophe Ritzenthaler, and Jeroen Sijsling
*Reconstructing plane quartics from their invariants*
[Discrete Comput. Geom. 63 (2020), no. 1, 73–113](https://doi.org/10.1007/s00454-018-0047-4)
