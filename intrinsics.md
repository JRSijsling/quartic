Exported intrinsics
--

### Invariants

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

intrinsic DiscriminantOfTernaryQuartic(f::RngMPolElt :
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

### Isomorphisms

```
intrinsic IsIsomorphicPlaneQuartics(X1::CrvPln, X2::CrvPln :
    geometric := false) -> BoolElt, SeqEnum
intrinsic IsomorphismsOfPlaneQuartics(X1::CrvPln, X2::CrvPln :
    geometric := false) -> SeqEnum

intrinsic IsIsomorphicTernaryQuartics(f1::RngMPolElt, f2::RngMPolElt :
    geometric := false) -> BoolElt, SeqEnum
intrinsic IsomorphismsOfTernaryQuartics(f1::RngMPolElt, f2::RngMPolElt :
    geometric := false) -> SeqEnum


intrinsic AutomorphismsOfPlaneQuartic(X::CrvPln :
    geometric := false) -> SeqEnum
intrinsic AutomorphismGroupOfPlaneQuartic(X::CrvPln :
    geometric := false, explicit := false) ->  GrpPerm, Map
intrinsic AutomorphismGroupOfPlaneQuartic(X::CrvPln, Autos::SeqEnum :
    explicit := false) ->  GrpPerm, Map

intrinsic AutomorphismsOfTernaryQuartic(f::RngMPolElt :
    geometric := false) -> SeqEnum
intrinsic AutomorphismGroupOfTernaryQuartic(f::RngMPolElt :
    geometric := false, explicit := false) ->  GrpPerm, Map
intrinsic AutomorphismGroupOfTernaryQuartic(f::RngMPolElt, Autos::SeqEnum :
    explicit := false) -> GrpPerm, Map

intrinsic GeometricAutomorphismGroup(C::Crv) -> GrpPerm
```

### Twists

```
intrinsic TwistsOfPlaneQuartic(C::Crv :
    AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm
intrinsic TwistsOfPlaneQuartic(C::Crv, Autos::SeqEnum  :
    AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm

intrinsic Twists(C::Crv :
    AutomorphismGroup := false) -> SeqEnum, GrpPerm
```
