Exported intrinsics
--

### Calculate invariants

```
intrinsic DixmierOhnoInvariants(f::RngMPolElt, p::RngIntElt :
    normalize := false,
    PrimaryOnly := false,
    IntegralNormalization := false,
    degmax := Infinity(), degmin := 1,
    PolynomialOnly := false) -> SeqEnum, SeqEnum
intrinsic DixmierOhnoInvariants(f::RngMPolElt :
    normalize := false,
    IntegralNormalization := false,
    PrimaryOnly := false, degmax := 10^6, degmin := 1,
    PolynomialOnly:=true) -> SeqEnum, SeqEnum
intrinsic DixmierOhnoInvariants(C::Crv :
    normalize := false,
    IntegralNormalization := false,
    PrimaryOnly := false, degmax := 10^6, degmin := 1,
    PolynomialOnly:=true) -> SeqEnum, SeqEnum

intrinsic QuarticDiscriminant(f::RngMPolElt :
    IntegralNormalization := false) -> Any
intrinsic DiscriminantFromDixmierOhnoInvariants(DO::SeqEnum) -> .

intrinsic DixmierOhnoInvariantsEqual(DO1::SeqEnum, DO2::SeqEnum) -> BoolElt

intrinsic DixmierOhnoAlgebraicRelations(DOinv::SeqEnum) -> SeqEnum
```

### Reconstruct curve from invariants

```
intrinsic PlaneQuarticFromDixmierOhnoInvariants(DO::SeqEnum :
    exact := false, minimize := true, descent := true, search_point := true) -> Crv, SeqEnum
intrinsic TernaryQuarticFromDixmierOhnoInvariants(DO::SeqEnum :
    exact := false, minimize := true, descent := true, search_point := true) -> RngMPolElt, SeqEnum
```

### Compute twists

```
intrinsic Twists(C::Crv :
    AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm
```

### Isomorphisms

```
intrinsic IsIsomorphicQuartic(X1::CrvPln, X2::CrvPln :
    geometric := false) -> .
intrinsic IsIsomorphicQuartic(f1::RngMPolElt, f2::RngMPolElt :
    geometric := false) -> .

intrinsic AutomorphismGroupQuartic(X::CrvPln :
    geometric := false) -> .
intrinsic AutomorphismGroupQuartic(f::RngMPolElt :
    geometric := false) -> .

intrinsic GeometricAutomorphismGroup(C::Crv) -> GrpPerm
```
