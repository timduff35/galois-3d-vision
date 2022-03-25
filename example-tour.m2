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
restart
FF = QQ
RNG = FF[r_(1,1)..r_(3,3), a_1..a_4, b_1..b_4, t_1, t_2, SAT]
R = genericMatrix(RNG,3,3)
t = matrix{{t_1},{t_2},{1}} -- work in a chart on PP^10
xs = for i from 1 to 4 list random(FF^3, FF^1)
ys = for i from 1 to 4 list random(FF^3, FF^1)
-- ideal of point correspondence equations
Icorr = ideal for i from 1 to 4 list (
    a_i * R * xs#(i-1) + t - b_i * sub(ys#(i-1),RNG)
    )
Irot = ideal(R*transpose R - id_(FF^3), det R - 1)
Iplanar = ideal det fold(for i from 1 to 4 list (a_i*sub(xs#(i-1),RNG))||matrix{{1}}, (a,b) -> a|b)
Isat = ideal(SAT * product for i from 1 to 4 list a_i*b_i - 1)
I = Icorr + Iplanar + Irot + Isat
-- next line timed at 7198s on MBP
elapsedTime G = groebnerBasis(I, Strategy => "F4");
inI = ideal leadTerm G;
dim inI, degree inI


-*
Section 4.3
*-

