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

### Isomorphisms

```
intrinsic IsIsomorphicQuartic(X1::CrvPln, X2::CrvPln :
    geometric := false) -> BoolElt, SeqEnum
intrinsic QuarticIsomorphisms(X1::CrvPln, X2::CrvPln :
    geometric := false) -> SeqEnum

intrinsic IsIsomorphicQuartic(f1::RngMPolElt, f2::RngMPolElt :
    geometric := false) -> BoolElt, SeqEnum
intrinsic QuarticIsomorphisms(f1::RngMPolElt, f2::RngMPolElt :
    geometric := false) -> SeqEnum


intrinsic QuarticAutomorphisms(X::CrvPln :
    geometric := false) -> SeqEnum
intrinsic AutomorphismGroupQuartic(X::CrvPln :
    geometric := false, explicit := false) ->  GrpPerm, Map
intrinsic AutomorphismGroupQuartic(X::CrvPln, Autos::SeqEnum :
    explicit := false) ->  GrpPerm, Map

intrinsic QuarticAutomorphisms(f::RngMPolElt :
    geometric := false) -> SeqEnum
intrinsic AutomorphismGroupQuartic(f::RngMPolElt :
    geometric := false, explicit := false) ->  GrpPerm, Map
intrinsic AutomorphismGroupQuartic(f::RngMPolElt, Autos::SeqEnum :
    explicit := false) -> GrpPerm, Map

intrinsic GeometricAutomorphismGroup(C::Crv) -> GrpPerm
```

### Twists

```
intrinsic QuarticTwists(C::Crv, Autos::SeqEnum  :
    AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm
intrinsic QuarticTwists(C::Crv :
    AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm

intrinsic Twists(C::Crv :
    AutomorphismGroup := false) -> SeqEnum, GrpPerm
intrinsic Twists(C::Crv, Autos::SeqEnum  :
    AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm
```
