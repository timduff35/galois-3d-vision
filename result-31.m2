-*
To run this script, start a Macaulay2 session in the directory containing it and run lines 48--56 (immediately following the "end" statement.)
*-
needs "common.m2"
-- options for path-tracking and monodromy
msOptions={Verbose=>true,NumberOfNodes=>3}

noLines = (nlines == 0)
noPoints = (npoints == 0)

-- setting up absolute pose equations
Rot = matrix for i from 1 to 3 list for j from 1 to 3 list declareVariable r_(i,j)
Cam=Rot|transpose matrix{for i from 1 to 3 list declareVariable t_i}
pts3d = for i from 1 to npoints list transpose matrix{for j from 1 to 4 list declareVariable P_(i,j)}
pts2d = for i from 1 to npoints list transpose matrix{for j from 1 to 3 list declareVariable Q_(i,j)}
lns3d = for i from 1 to nlines list matrix{for j from 1 to 4 list declareVariable L1_(i,j)} || matrix{for j from 1 to 4 list declareVariable L2_(i,j)}
lns2d = for i from 1 to nlines list matrix{for j from 1 to 3 list declareVariable M_(i,j)}
FF=CC
ptParam = transpose(rfold pts3d || rfold pts2d)
lnParam = if noLines then random(FF^0, FF^0) else cfold apply(lns3d,l->matrix{flatten entries l}) | cfold lns2d
Param = if noPoints then lnParam else if nlines == 0 then ptParam else ptParam | lnParam
X = gateMatrix{flatten entries Cam}
eqs=flatten entries(Rot*transpose Rot - id_(CC^3))|
   flatten apply(pts3d,pts2d,(p,q) -> allMinors(Cam*p|q,2,Laplace=>true)) |
   flatten apply(lns3d,lns2d,(l,m) -> allMinors(l||m*Cam,3,Laplace=>true))
G=gateSystem(Param,X,transpose gateMatrix{eqs})

--sampling the incidence variety
R0=(-1)*randomOn 3
t0 = random(CC^3,CC^1)
Cam0 = R0|t0
pts3d0 = for i from 1 to npoints list random(CC^4,CC^1)
pts2d0 = for i from 1 to npoints list Cam0* pts3d0#(i-1)
lns3d0 = for i from 1 to nlines list random(CC^2,CC^4)
lns2d0 = apply(lns3d0, l -> (
        p1p2 := numericalKernel(l,Tolerance=>1e-5);
        transpose numericalKernel(transpose(Cam0*p1p2),Tolerance=>1e-5)
        )
    )
lnParam0 = if noLines then random(FF^0, FF^0) else cfold apply(lns3d0,l->matrix{flatten entries l}) | cfold lns2d0
ptParam0 = transpose(rfold pts3d0 || rfold pts2d0)
Param0 = point if noPoints then lnParam0 else if nlines == 0 then ptParam0 else ptParam0 | lnParam0
X0= point matrix{flatten entries Cam0}
evalX0=evaluate(G,Param0,X0)
resid=norm evalX0
assert areEqual(0,resid)
Gsq = squareDown(Param0,X0,numVariables G,G)
end--
restart
needsPackage "MonodromySolver"
for i from 0 to 3 do (
    npoints = i;
    nlines = 3-npoints;
    setRandomSeed 2022;
    load "result-31.m2";
    monodromyGroup(Gsq, Param0, {X0}, "msOptions"=>msOptions, FileName=>"./groups-absolute-pose/npts-"|toString(npoints)|"nlns"|toString(nlines)|".gp");
    )
