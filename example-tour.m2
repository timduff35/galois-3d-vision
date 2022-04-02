-*
Example 1.1
*-

-*
Example 2.4
*-

-*
Examples 2.7 and 2.13
*-

-*
Example 2.11
*-

-*
Example 2.12
*-


-*
Section 3.1
*-

-*
Section 3.2
*-

-*
Section 4.1
A dominant, rational map from (P^2)^20  - -> Grass(P^3, P^8)
*-

restart
needs "common.m2"
--needsPackage "MonodromySolver"
xs = for i from 1 to 5 list gateMatrix{{declareVariable x_(i,1),declareVariable x_(i,2),inputGate 1}}
ys = for i from 1 to 5 list gateMatrix{{declareVariable y_(i,1),declareVariable y_(i,2),inputGate 1}}
L = fold(apply(xs,ys,(x,y) -> gateMatrix{flatten entries(transpose x*y)}),(a,b)->a||b)
(L1, L2) = (L_{0..4}, L_{5..8})
stiefelL = (inverse L1) * L2;
gateOutput = transpose gateMatrix{flatten entries(stiefelL)};
gateInput = gateMatrix fold(apply(xs, x -> x_{0,1}), (a,b) -> a|b)| fold(apply(ys, y -> y_{0,1}), (a,b) -> a|b),
stiefelGateSystem = gateSystem(gateInput, gateOutput)
xy0 = random(QQ^20,QQ^1) -- try replacing with your favorite point correspondences (x1,x2,..,x5,y1,..,y5)
printWidth = 10000
J = evaluateJacobian(stiefelGateSystem, point xy0)
rank J

-*
Section 4.2
*-

-*
First, we verify the degree of the problem using Groebner bases over a finite field.
Here, he variable "SAT" enforces the constraint that depths are nonzero.
*-
restart
FF = ZZ/nextPrime 2022
RNG = FF[r_(1,1)..r_(3,3), a_1..a_4, b_1..b_4, t_1..t_3, SAT]
R = transpose genericMatrix(RNG,3,3)
tt = matrix{{t_1},{t_2},{t_3}}
xs = for i from 1 to 4 list random(FF^3, FF^1)
ys = for i from 1 to 4 list random(FF^3, FF^1)
-- ideal of point correspondence equations
Icorr = ideal for i from 1 to 4 list (
    a_i * R * xs#(i-1) + tt - b_i * sub(ys#(i-1),RNG)
    )
Irot = ideal(R*transpose R - id_(FF^3), det R - 1)
worldPoints = fold(for i from 1 to 4 list (a_i*sub(xs#(i-1),RNG))||matrix{{1}}, (a,b) -> a|b)
Iplanar = ideal det worldPoints
Isat = ideal((SAT * product for i from 1 to 4 list a_i*b_i) - 1)
I = Icorr + Iplanar + Irot + Isat + ideal(t_3-1)-- We may verify the degree by working in the affine chart on depths/translation space where t3 = 1.
-- next line timed at ~.25s on MBP
elapsedTime G = groebnerBasis(I, Strategy => "F4");
inI = ideal leadTerm G;
dim inI, degree inI
-*
 We now verify our formula for the deck transformation Psi2 as follows:
   * express the plane normal n as a function of input data and unknown quantities
   * check that application of Psi2 preserves homogeneous equations defining I in the function field of X.
   (Here, we switch to homogeneous coordinates on depths/translation space.)
*-
nd = lift((gens ker sub(transpose worldPoints, RNG/Iplanar))_{0}, RNG)
assert((transpose nd * worldPoints) % Iplanar == 0)
FRNG = frac(RNG)
nd = promote(nd, FRNG)
(n,d) = (nd^{0,1,2}, -nd_(3,0))
ntn = (transpose n *n)_(0,0)
nnt = n * transpose n
tt = promote(tt, FRNG)
R = promote(R, FRNG)
Psi2R = R * ((2/ntn)*nnt - id_(FRNG^3))
Psi2t = -tt - (2*d/ntn)*R*n
psi2 = map(FRNG, FRNG, (flatten entries Psi2R)|{a_1,a_2,a_3,a_4,-b_1,-b_2,-b_3,-b_4} | (flatten entries Psi2t)| {SAT});
assert(gens FRNG == psi2 \ psi2 \ (gens FRNG)) -- Psi2 has order 2

-- next two lines verify that Psi2 preserves planarity of the scene
Iplanar = promote(Iplanar, FRNG)
psi2 Iplanar 
-- next two lines verify that Psi2 preserves point correspondence constraints.
Icorr = promote(Icorr, FRNG)
apply(numerator \ (psi2 Icorr)_*, e -> e % (ideal G))
-- next two lines verify that Psi2 preserves SO(3) constraints.
Irot = promote(Irot, FRNG)
apply(numerator \ (psi2 Irot)_*, e -> e % (ideal G))

-- Now we verify that the composition of Psi1 and Psi2 preserves the homography matrix H
normSquared = x -> (assert(numcols x == 1); (transpose x * x)_(0,0))
Psi1a = for i from 1 to 4 list -a_i * (normSquared tt)/(normSquared(b_i*promote(ys#(i-1),FRNG))-normSquared(a_i*promote(xs#(i-1),FRNG)))
Psi1b = for i from 1 to 4 list  b_i * (normSquared tt)/(normSquared(b_i*promote(ys#(i-1),FRNG))-normSquared(a_i*promote(xs#(i-1),FRNG)))
Psi1R = ((2/normSquared(tt)) * tt * transpose tt - id_(FRNG^3)) * R
Psi1tt = tt
psi1 = map(FRNG, FRNG, (flatten entries Psi1R)| Psi1a | Psi1b | (flatten entries Psi1tt)| {SAT});
H = R + (1/d)*tt*transpose n
H12 = psi1 psi2 H;
U = unique(denominator \ (flatten entries H))
U12 = unique(denominator \ (flatten entries H12))
assert(#U == 1 and #U12 == 1)
-- Timings recorded for running the next line were on the order of 2min
elapsedTime HminusH12Polynomial = ((first U)*matrix apply(entries H12, r -> numerator \ r) -  (first U12)*matrix apply(entries H, r -> numerator \ r));
assert(HminusH12Polynomial % (ideal G) == 0)

-*
Section 4.3
*-

