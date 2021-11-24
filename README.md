Description
--

This repository contains Magma code for reconstructing plane quartic curves from their Dixmier--Ohno invariants and for calculating isomorphisms between plane quartic curves. The latter code only works over the rationals and finite fields. In some non-generic cases, we use an approach due to Andrew Sutherland that involves Groebner bases.

This package is now included by default in Magma-v2.26.

Prerequisites
--

An installation of Magma.

The algorithms use the Magma algorithm `MinimizeReducePlaneQuartic`, due to Elsenhans and Stoll, to simplify any output over the rationals. In order for this to work properly, one needs a bug fix of a function at the end of the file `magma/package/Geometry/SrfDP/minred.m`, which can be found in `magma/BugFix.m`. Be warned that this new file is attached by default and supersedes its version in Magma.

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

Reynald Lercier, Jeroen Sijsling, Christophe Ritzenthaler,
*Functionalities for genus 2 and 3 curves*,
[arXiv:2102.04372](https://arxiv.org/abs/2102.04372)
