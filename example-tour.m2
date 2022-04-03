-*
Example 1.1
*-

-*
Example 2.4
*-
restart
needsPackage "MonodromySolver"
unknownMatrix = gateMatrix{{declareVariable x}}
parameterMatrix = gateMatrix{{declareVariable p}}
equationMatrix = gateMatrix{{x^3-p}}
G = gateSystem(parameterMatrix, unknownMatrix, equationMatrix) -- of class "GateSystem"
setRandomSeed 2022
monodromyGroup(G, "msOptions" => {NumberOfNodes => 5}, FileName => "example-21.gp")

-*
Examples 2.8 and 2.14
*-
restart
FF = ZZ/nextPrime 2022
RNG = FF[a_0..a_3,b_0..b_3]
E = matrix{ -- coordinates of the map X -> E
    {
        -- first row
        a_0*b_0 - a_1*b_1 - a_2*b_2 + a_3*b_3, 
        a_0*b_1 + a_1*b_0 + a_2*b_3 + a_3*b_2,
        a_0*b_2 + a_2*b_0 - a_1*b_3 - a_3*b_1
        },
    {
        -- second row
        a_0*b_1 + a_1*b_0 - a_2*b_3 - a_3*b_2, 
        -a_0*b_0 + a_1*b_1 - a_2*b_2 + a_3*b_3,
        a_1*b_2 + a_2*b_1 + a_0*b_3 + a_3*b_0
        },
    {
        -- third row
        a_0*b_2 + a_2*b_0 + a_1*b_3 + a_3*b_1, 
        a_1*b_2 + a_2*b_1 - a_0*b_3 - a_3*b_0,
        -a_0*b_0 - a_1*b_1 + a_2*b_2 + a_3*b_3
        }

}
IX = ideal sum for i from 0 to 3 list a_i * b_i
-- Check that the matrix E satisfies Demazure's trace constraints:
(E * transpose E * E - (1/2) * trace(E * transpose E) * E) % IX
-- Check that swapping a & b has the same effect as the twisted pair:
phi = map(RNG, RNG, {b_0,b_1,b_2,b_3,a_0,a_1,a_2,a_3})
E-phi E
-- Verify (in the affine chart a_0 = b_0 = 1 on P^3 x P^3) that we get 20 = deg(X/Z) solutions for (a,b)
I = IX + ideal apply(5, i -> random(FF^1,FF^3) * E * random(FF^3, FF^1));
elapsedTime I = saturate(I, ideal(a_0..a_3)*ideal(b_0..b_3)) + ideal(a_0-1,b_0-1);
dim I, degree I

-*
Section 3.1
*-
-- Galois group of equations (16) is the full wreath product...
restart
needsPackage "MonodromySolver"
unknownMatrix = gateMatrix{toList vars(x_1,x_2,x_3)}
parameterMatrix = gateMatrix{toList vars(A,B,C,D,E,F,G,H,I,J,K,L)}
sparseEquationMatrix = transpose gateMatrix{{A*x_1^2+B*x_2^2+C*x_1*x_2+D,E*x_1^2+F*x_3^2+G*x_1*x_3+H,I*x_2^2+J*x_3^2+K*x_2*x_3+L}}
sparseG = gateSystem(parameterMatrix, unknownMatrix, sparseEquationMatrix)
setRandomSeed 2022
monodromyGroup(sparseG, "msOptions" => {NumberOfNodes => 5}, FileName => "sparsified-P3P.gp")

-- ...but the Galois group of P3P (with the same monomial support) is that wreath product's intersection with A8!
P3PEquationMatrix = transpose gateMatrix{{A*(x_1^2+x_2^2)+C*x_1*x_2+D,E*(x_1^2+x_3^2)+G*x_1*x_3+H,I*(x_2^2+x_3^2)+K*x_2*x_3+L}}
P3PG = gateSystem(parameterMatrix, unknownMatrix, P3PEquationMatrix)
setRandomSeed 2022
monodromyGroup(P3PG, "msOptions" => {NumberOfNodes => 5}, FileName => "P3P.gp")

-- Computing the discriminant of "Grunert's equations" for the P3P absolute pose problem.
-- Remark: dij below represent squared distances: however, reparametrizing in terms of distances does not change the Galois group.
restart
tuples = {(1,2), (1,3), (2,3)}
FF = frac(QQ[flatten \\ tuples/(t-> {d_t, c_t})])
R = FF[x_1,x_2,x_3]
I = ideal apply(tuples, t -> (
	(i, j) := (first t, last t);
    	x_i^2 + x_j^2 - c_(i,j)*x_i*x_j - d_(i,j)
	)
    )
elapsedTime G = groebnerBasis(I, Strategy => "F4");
needsPackage "FGLM"
-- The next line was timed at ~ 115s on a 2021 Macbook pro w/ 8G RAM
elapsedTime lexG = fglm(forceGB G, FF[gens R,MonomialOrder => Lex]);
netList flatten entries gens lexG
-- Extract monomials and coefficients from the minimal polynomial:
(m,c) = coefficients (gens lexG)_(0,0)
-- Since the constant term of the minimal polynomial is a square, its discriminant is a square:
factor c_(numrows c -1, 0)

-*
Section 3.2: Investigating symmetries of absolute pose problems with mixed point & line features.
We use the constraints proposed in the paper:
  "Pose Estimation using Both Points and Lines for Geo-Localization" (Ramalingam, Bouaziz, Sturm, ICRA 2011.)
*-
-- 2 points and 1 line (Equations 18 in our paper, Equations 4 and 5 in RBS '11)
restart
FF = ZZ/nextPrime 2022
RNG = FF[b_1,b_2,a_1,a_2,X_2,X_3,X_4,Y_3,Y_4,Z_4,R_(1,1), R_(2,1), R_(3,1), R_(2,2), R_(2,3), T_1,T_2,T_3]
XX = transpose matrix{{R_(1,1), R_(2,1), R_(3,1), R_(2,2), R_(2,3), T_1,T_2,T_3}}
B = matrix{{0},{-b_1},{0},{-b_2},{0},{0}}
A = matrix{
    {       0,       0,       0,   0,   0, -b_1, a_1,   0},
    {       0,       0,       0,   0,   0,    0,  -1, b_1},
    {-b_2*X_2, a_2*X_2,       0,   0,   0, -b_2, a_2,   0},
    {       0,    -X_2, b_2*X_2,   0,   0,    0,  -1, b_2},
    {       0,     X_3,       0, Y_3,   0,    0,   1,   0},
    {       0,     X_4,       0, Y_4, Z_4,    0,   1,   0}
    }
I = ideal(A*XX-B, R_(1,1)^2 + R_(2,1)^2 + R_(3,1)^2 - 1, R_(2,1)^2 + R_(2,2)^2 + R_(2,3)^2 - 1)
-- The deck transformation (R,t) -> (-R, t-2e3) is easily verified
eq1 = A*XX-B
deckXX = -XX - matrix apply(numrows XX, i -> {if i < numrows XX -1 then 0 else 2})
eq2 = A*deckXX-B
assert(eq1 + eq2 == 0)
-- Verify the degree of the problem on random input
specializedRing = FF[flatten entries XX]
specializationMap = map(specializedRing, RNG, apply(numgens RNG - numrows XX, i -> random FF) | (flatten entries sub(XX, specializedRing)))
specializedI = specializationMap I
dim specializedI, degree specializedI

-- 1 point and 2 lines (Equations 18 in our paper, Equations 7 and 8 in RBS '11.)
restart
FF = ZZ/nextPrime 2022
RNG = FF[b_1,b_2,b_4,a_1,a_2,a_4,X_1..X_4,Y_1..Y_4,Z_1..Z_4,R_(1,1),R_(1,2),R_(1,3),R_(2,1),R_(2,2),R_(2,3),T_1,T_2,T_3]
XX = transpose matrix{{R_(1,1),R_(1,2),R_(1,3),R_(2,1),R_(2,2),R_(2,3),T_1,T_2,T_3}}
B = matrix{{0},{-b_1},{0},{0},{0},{0}}
A = transpose matrix{
    {   0,   0,   0,   0, -b_4*X_3, -b_4*X_4},
    {   0,   0,   0,   0, -b_4*Y_3, -b_4*Y_4},
    {   0,   0,   0,   0, -b_4*Z_3, -b_4*Z_4},
    {   0,   0, X_1, X_2,  a_4*X_3,  a_4*X_4},
    {   0,   0, Y_1, Y_2,  a_4*Y_3,  a_4*Y_4},
    {   0,   0, Z_1, Z_2,  a_4*Z_3,  a_4*Z_4},
    {-b_1,   0,   0,   0,     -b_4,     -b_4},
    { a_1,  -1,   1,   1,      a_4,      a_4},
    {   0, b_1,   0,   0,        0,        0}
    }
I = ideal(A*XX-B, R_(1,1)^2 + R_(1,2)^2 + R_(1,3)^2 - 1, R_(2,1)^2 + R_(2,2)^2 + R_(2,3)^2 - 1, R_(1,1)*R_(2,1) + R_(1,2) * R_(2,2) + R_(1,3) * R_(2,3))
-- The deck transformation (R,t) -> (-R, t-2e3) is easily verified
eq1 = A*XX-B
deckXX = -XX - matrix apply(numrows XX, i -> {if i < numrows XX -1 then 0 else 2})
eq2 = A*deckXX-B
assert(eq1 + eq2 == 0)
-- Verify the degree of the problem on random input:
specializedRing = FF[flatten entries XX]
specializationMap = map(specializedRing, RNG, apply(numgens RNG - numrows XX, i -> random FF) | (flatten entries sub(XX, specializedRing)))
specializedI = specializationMap I
dim specializedI, degree specializedI

-*
Section 4.1
A dominant, rational map from (P^2)^20  - -> Grass(P^3, P^8) appearing in Proposition 4.1.
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
assert(rank J == 20)

-*
Section 4.2
*-

-*
First, we verify the degree of the problem using Groebner bases over a finite field.
Here, the variable "SAT" enforces the constraint that depths are nonzero.
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
Psi2 = map(FRNG, FRNG, (flatten entries Psi2R)|{a_1,a_2,a_3,a_4,-b_1,-b_2,-b_3,-b_4} | (flatten entries Psi2t)| {SAT});
assert(gens FRNG == Psi2 \ Psi2 \ (gens FRNG)) -- Psi2 has order 2

-- next two lines verify that Psi2 preserves planarity of the scene
Iplanar = promote(Iplanar, FRNG)
Psi2 Iplanar 
-- next two lines verify that Psi2 preserves point correspondence constraints.
Icorr = promote(Icorr, FRNG)
apply(numerator \ (Psi2 Icorr)_*, e -> e % (ideal G))
-- next two lines verify that Psi2 preserves SO(3) constraints.
Irot = promote(Irot, FRNG)
apply(numerator \ (Psi2 Irot)_*, e -> e % (ideal G))

-- Now we verify that the composition of Psi1 and Psi2 preserves the homography matrix H
load "common.m2"
Psi1a = for i from 1 to 4 list -a_i * (normSquared tt)/(normSquared(b_i*promote(ys#(i-1),FRNG))-normSquared(a_i*promote(xs#(i-1),FRNG)))
Psi1b = for i from 1 to 4 list  b_i * (normSquared tt)/(normSquared(b_i*promote(ys#(i-1),FRNG))-normSquared(a_i*promote(xs#(i-1),FRNG)))
Psi1R = ((2/normSquared(tt)) * tt * transpose tt - id_(FRNG^3)) * R
Psi1tt = tt
Psi1 = map(FRNG, FRNG, (flatten entries Psi1R)| Psi1a | Psi1b | (flatten entries Psi1tt)| {SAT});
H = R + (1/d)*tt*transpose n
H12 = Psi1 Psi2 H;
U = unique(denominator \ (flatten entries H))
U12 = unique(denominator \ (flatten entries H12))
assert(#U == 1 and #U12 == 1)
elapsedTime HminusH12Polynomial = ((first U)*matrix apply(entries H12, r -> numerator \ r) -  (first U12)*matrix apply(entries H, r -> numerator \ r));
assert(HminusH12Polynomial % (ideal G) == 0)

-- It is far easier to check that Psi2(H) = -H
H2 = Psi2 H
assert(H + H2 == 0)

-- Finally, we verify dimension and degree of the ideal in Equation (29)
R = FF[s,S_(1,2),S_(1,3),S_(2,1)..S_(3,3)]
Smatrix = matrix for i from 1 to 3 list for j from 1 to 3 list if (i==1 and j==1) then 1 else S_(i,j)
Icorr = ideal apply(cpMatrix \ ys, xs, (yx, x) -> yx * Smatrix * x )
I = Icorr + ideal(det(transpose Smatrix * Smatrix - s*id_(R^3)))
(dim I, degree I)

-*
Section 4.3
*-
restart
FF = ZZ/nextPrime 2022
RNG = FF[H_(1,1,1)..H_(2,3,3)]
H1 = transpose genericMatrix(RNG,3,3)
H2 = transpose genericMatrix(RNG,H_(2,1,1),3,3)
W1 = transpose H1 * H1 - id_(RNG^3)
W2 = transpose H2 * H2 - id_(RNG^3)
-- I31 = <equations (31)>. 
I31 = ideal(det W1, det W2, det(W1+W2), det(W1-W2));
-- We may verify that the non-reduced affine scheme Spec(R/I31) has the expected degree 6^4=1296 using Groebner bases:
IGenericLinear = ideal apply(14,i->random(1,RNG)-1)
G = groebnerBasis(I31 + IGenericLinear, Strategy=>"F4");
inG = ideal leadTerm G;
dim inG, degree inG
-- Enforcing additional constraints in Equation (32) removes non-reduced solutions, verifying deg(H2) = 336.
Iresultant = (
    RNG12 = FF[gens RNG,n1,n2];
    (W1, W2) = (sub(W1, RNG12), sub(W2, RNG12));
    ideal apply({
	resultant(W1_(2,2) * n1^2 - 2*W1_(0,2) * n1 + W1_(0,0), W2_(2,2) * n1^2 - 2*W2_(0,2) * n1 + W2_(0,0), n1),
	resultant(W1_(2,2) * n2^2 - 2*W1_(1,2) * n2 + W1_(1,1), W2_(2,2) * n2^2 - 2*W2_(1,2) * n2 + W2_(1,1), n2)
	}, r -> sub(r, RNG))
    );
G = groebnerBasis(I31 + IGenericLinear + Iresultant, Strategy=>"F4");
inG = ideal leadTerm G;
dim inG, degree inG

-- Point correspondence constraints, however, are not generic linear equations:
load "common.m2"
xs = (xs := for i from 1 to 3 list random(FF^3, FF^1); xs = xs | {xs#0 + (random FF) * xs#1})
ys = (ys := for i from 1 to 3 list random(FF^3, FF^1); ys = ys | {ys#0 + (random FF) * ys#1})
zs = (zs := for i from 1 to 3 list random(FF^3, FF^1); zs = zs | {zs#0 + (random FF) * zs#1})
Icor = ideal(
    apply(cpMatrix \ ys, xs, (yx, x) -> yx * H1 * x ) | 
    apply(cpMatrix \ zs, xs, (zx, x) -> zx * H2 * x )
    )
-- if we only enforce constraints in Eq 31, the minimal problem has 320 "solutions"
G = groebnerBasis(I31 + Icor, Strategy=>"F4");
inG = ideal leadTerm G;
dim inG, degree inG
-- if we only enforce constraints in Eqs 31 and 32, the minimal problem has 320 "solutions"
G = groebnerBasis(I31 + Iresultant + Icor, Strategy=>"F4");
inG = ideal leadTerm G;
dim inG, degree inG
-- The last constraint we need is det H1 != 0: in fact, we should saturate by the square of det H1 to obtain Result 4.3
Isat = (ideal G):(det H1)^2;
dim Isat, degree Isat -- Degree 64 matches the degree computed in rotation/translations, justifying our claim that the associated branched covers are birationally equivalent.
