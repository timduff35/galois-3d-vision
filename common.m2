-- IMPORTS 
needsPackage "MonodromySolver"

-- FUNCTIONS

size GateMatrix := M -> (numrows M, numcols M)
size Matrix := M -> (numrows M, numcols M)

-*
-- evaluate a gateMatrix G at a matrix x
-- don't use this a lot...
evaluate (GateMatrix, Matrix) := (G, x) -> (
    M := mutableMatrix(FF,numrows G,numcols G);
    E' := makeEvaluator(G,matrix{cameraVars|dataParams});
    evaluate(E',mutableMatrix(x),M);
    matrix M
    )
*-

--random diagonal matrix
randDiag = n -> diagonalMatrix for i from 1 to n list random CC

dehomogenize = method(Options=>{})
dehomogenize (Matrix, ZZ) := o -> (v, n) -> (
    --assumes column vector
    (1/v_(n, 0))*v^(toList splice(0..n-1,n+1..numrows v-1))
    )
dehomogenize Matrix := o -> v -> dehomogenize(v, numrows v -1)

summary = L -> (
    n := #L;
    H := sort L;
    Q1 := (1/2) * (H#(floor((n-1)/3))+H#(ceiling((n-1)/3)));
    med := (1/2) * (H#(floor((n-1)/2))+H#(ceiling((n-1)/2)));
    Q3 := (1/2) * (H#(floor(2*(n-1)/3))+H#(ceiling(2*(n-1)/3)));
    mean := (sum L)/n;
    var := sum(L/(x-> (x - mean)^2))/(n-1);
    << "Min: " << toString(min L) << endl;
    << "1Q: " << toString(Q1) << endl;
    << "Med: " << toString(med) << endl;
    << "Avg: " << toString(sub(mean,RR)) << endl;
    << "3Q: " << toString(Q3) << endl;
    << "Max: " << toString(max L) << endl;
    << "Std Dev: " << toString(sqrt(var)) << endl;
    )    

-- random element in the kernel of M
randKernel = method(Options=>{Tolerance=>1e-4})
randKernel (Matrix, InexactFieldFamily) := o -> (M, FF) -> (
    K := numericalKernel(M, o.Tolerance);
    K*random(FF^(numcols K), FF^1)
    )
randKernel Matrix := o -> M -> randKernel(M, CC)

reshapeCol = p -> if (numrows p == 1) then transpose p else p
reshapeRow = p -> if (numcols p == 1) then transpose p else p

-- RANDOMIZATION FOR PARAMETER POINT p
-- assumes random CC has unit modulus!
gammify = method()
gammify Point := p -> gammify reshapeCol matrix p
gammify Matrix := p -> (
    gammas := for i from 1 to m*#indepLines list random CC;
    indDiag := flatten(gammas/(g->{g,g,g}));
    -- abstract to arbitrary diagram, number of cameras
    depWInd := depLines/(l->first select(1,last D,i->member(l,i)));
    mfoldIntersections := flatten(depWInd/(l->for i from 0 to m-1 list l/(x->m*x+i)));
    -- next line assumes indepenedent lines come first!
    depDiag := flatten(
	mfoldIntersections/(ind -> {
	    conjugate gammas#(ind#0), 
	    conjugate gammas#(ind#1)
	    }
	)
    );
    tChartdiag := toList((3*(m-1)+1):random(CC)); -- t chart gamma
    qChartDiags := flatten for i from 0 to m-2 list toList(5:random(CC));
    p' := diagonalMatrix(indDiag|depDiag|tChartdiag|qChartDiags)*p;
    p'
    )


-- fold along rows
rfold = L -> if (#L ==0) then random(FF^0,FF^0) else fold(L, (a,b) -> a||b)

-- fold along cols
cfold = L -> fold(L, (a,b) -> a|b)


-- write starting parameters and solutions to file
writeStartSys = method(Options=>{Filename=>"startSys"})
writeStartSys (Matrix, List) := o -> (M, sols) -> writeStartSys(point M, sols, o)
writeStartSys (Point, List) := o -> (p, sols) -> (
   assert(instance(o.Filename,String));
   f := openOut o.Filename;
   f << "Parameter values: " << endl;
   f << toExternalString p << endl;
   f << "Solutions : " << endl;
   for s in sols do f << toExternalString s << endl;
   close f;
   )

readStartSys = filename -> (
    l := separate("\n", get filename);
    p0 := value l#1;
    sols := for i from 3 to #l-2 list value l#i;
    (transpose matrix p0, sols/matrix/transpose)
    )

-- for testing the contents of a start system file
startSysTester = (p,sols) -> (
    p0 := (transpose matrix V.BasePoint);
    p1 := random(CC^(#dataParams),CC^1);
    P01 = p0||p1;
    Pspec01 := specialize(PH,P0);
    target01 := trackHomotopy(Pspec01, sols);
    Pspec10 := (gammify p1)|(gammify p0);
    trackHomotopy(Pspec10, target01)
    )
    
-- not printing to high precision -- deprecated?
sol2String = p -> replace("\\{|\\}","",toString p.Coordinates)

-- produces gates for "small" determinants"
det2 = M -> M_(0,0)*M_(1,1)-M_(1,0)*M_(0,1)
det3 = M -> M_(0,0)*det2(M_{1,2}^{1,2})-M_(0,1)*det2(M_{0,2}^{1,2})+M_(0,2)*det2(M_{0,1}^{1,2})
det4 = M -> M_(0,0)*det3(M_{1,2,3}^{1,2,3})-M_(0,1)*det3(M_{0,2,3}^{1,2,3})+M_(0,2)*det3(M_{0,1,3}^{1,2,3})-M_(0,3)*det3(M_{0,1,2}^{1,2,3})

laplaceDet = M -> (
    (m, n) := size M;
    if (m=!=n) then error("not square matrix")
    else if (m>5) then error("no Laplace for matrices larger than 4x4")
    else if (m==2) then det2 M
    else if (m==3) then det3 M
    else -* m==4 *- det4 M
    )

adj2 = X -> ( -- adjugate of a 2x2 GateMatrix
         assert(numrows X == 2 and numcols X == 2);
         matrix{{X_(1,1),-X_(0,1)},{-X_(1,0),X_(0,0)}}
         )

inv2 = X -> (inputGate 1)/(det2 X) * adj2 X
inverse GateMatrix := M -> (
    n := numrows M;
    assert(n == numcols M);
    if (n==1) then matrix{{(inputGate 1)/M_(0,0)}}
    else if (n==2) then inv2 M
    else ( -- use Schur complement
	k := floor(n/2);
	A := M_{0..k-1}^{0..k-1};
	B := M_{k..n-1}^{0..k-1};
	C := M_{0..k-1}^{k..n-1};
	D := M_{k..n-1}^{k..n-1};
	Dinv := inverse D;
	DinvC := Dinv * C;
	BDinv := B * Dinv;
	M'over'D := A - BDinv * C;
	M'over'Dinv := inverse M'over'D;
	(
	    (M'over'Dinv         | -M'over'Dinv * BDinv) ||
	    (-DinvC * M'over'Dinv | Dinv + DinvC*M'over'Dinv*BDinv)
	    )
	)
    )

-- jacobian of GateMatrix wrt. a list of inputGates
--jacobian (GateMatrix, List) := (F,inGates) -> fold(apply(inGates,g->diff(g,F)),(a,b)->a|b)

-- get rotation matrix from cayley parameters
cay2R = method(Options=>{Normalized=>false})
cay2R (Thing,Thing,Thing) := o -> (X,Y,Z) -> (
    if instance(X, RR) then x := X_CC else x = X;
    if instance(Y, RR) then y := Y_CC else y = Y;
    if instance(Z, RR) then z := Z_CC else z = Z;
    M := matrix{
    {1+x*x-(y*y+z*z), 2*(x*y-z), 2*(x*z+y)},
    {2*(x*y+z), 1+y^2-(x*x+z*z), 2*(y*z-x)},
    {2*(x*z-y), 2*(y*z+x), 1 +z*z -(x*x+y*y)}
	};
    if o.Normalized then (1/(1+x^2+y^2+z^2)) * M else M
    )
cay2R List := o -> L -> cay2R(L#0, L#1, L#2, o)

-- get Cayley parameters from rotation matrix
R2Cay = method(Options=>{UnNormalize=>false})
R2Cay Matrix := o -> R -> (
    assert(numcols R == 3);
    assert(numrows R == 3);
    S := (R-id_(CC^3))*(R+id_(CC^3))^-1;
    (S_(2,1), S_(0,2), S_(1,0))
    )

-*/// TEST
restart
needs "common.m2"
(x, y, z) = (random RR, random RR, random RR)
R = cay2R(x, y, z)
(x',y',z') = R2Cay R
R = cay2R(x', y', z')
R2Cay R
///*-

-- get rotation matrix from quaternion parameters
Q2R = method(Options=>{Normalized=>false, FF=>CC})
Q2R (Thing,Thing,Thing, Thing) := o -> (W, X,Y,Z) -> (
    if instance(W, RR) then w := W_CC else w = W;
    if instance(X, RR) then x := X_CC else x = X;
    if instance(Y, RR) then y := Y_CC else y = Y;
    if instance(Z, RR) then z := Z_CC else z = Z;
    M := matrix{
    {w*w+x*x-(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)},
    {2*(x*y+w*z), w^2+y^2-(x*x+z*z), 2*(y*z-w*x)},
    {2*(x*z-w*y), 2*(y*z+w*x), w^2 +z*z -(x*x+y*y)}
	};
    if o.Normalized then (1/(w^2+x^2+y^2+z^2)) * M else M
    )
Q2R List := o -> L -> Q2R(L#0, L#1, L#2, L#3, o)

-- get Cayley parameters from rotation matrix
R2Q = method(Options=>{UnNormalize=>false,FF=>CC})
R2Q Matrix := o -> R -> (
    assert(numcols R == 3);
    assert(numrows R == 3);
    c := (R_(2,1) - R_(1,2));
    b := (R_(0,2) - R_(2,0));
    a := (R_(1,0) - R_(0,1));
    w := (1/2)*sqrt(R_(0,0)+R_(1,1)+R_(2,2)+1);
    x := 1/(4*w) * c;
    y := 1/(4*w) * b;
    z := 1/(4*w) * a;
--    << w^2+x^2+y^2+z^2 << endl;
    (w, x, y, z)
    )

-*/// TEST
R=CC[W]
netList solveSystem {W^4-W^2+1/16}
clean T
T=QQ[a..d]
R=Q2R gens T
S = (R-id_(((QQ)^3)))*adjugate(R+id_((QQ)^3));
S
((first x)/(first L))*L
1/sqrt(sum(x/(y->y^2)))*x
L
///*-


-- cross product of col vectors -- takes Matrice or GateMatrix pair
crossProduct = (y,q) -> matrix{{y_(1,0)*q_(2,0)-y_(2,0)*q_(1,0)},{y_(2,0)*q_(0,0)-y_(0,0)*q_(2,0)},{y_(0,0)*q_(1,0)-y_(1,0)*q_(0,0)}}

cmat = method(Options=>{Normalized=>true})
cmat List := o -> L -> cmat(L#0,L#1,L#2,o)
cmat (Thing,Thing,Thing) := o -> (a,b,c) -> (
    tx := matrix{{0,-c,b},{c,0,-a},{-b,a,0}};
    if o.Normalized and not areEqual(L#2,0.0) then tx = (-1/L#2) * tx;
    tx
    )

--
randomLineThroughPoints = (P, FF) -> ( 
    m := numrows P; -- m = dim + 1
    n := numcols P; -- n = number of points
    assert(m>=3 and m<=4);
    K := numericalKernel(transpose P,1e-6);
    --assert(numcols K == m-n); -- true if points are distinct
    transpose(K * random(FF^(numcols K),FF^(m-2)))
    )

-- constructs problem data given a PL diagram D (complete visibility)
-- returns (camera parameters, lines, camera matrices
-- ASSUMES: "intersections" are sorted
fabricatePair = (D, FF, nvars) -> (    
    (nLines,nGhosts,intersections) := D;
    depPoints := set {};
    scan(nLines,l->(
	    ptsOnl := positions(last D,i->member(l,i));
	    depPoints = depPoints + set drop(ptsOnl,2);
	    )
	);    
    pointsOnLineIndices := apply(nLines+nGhosts, l->positions(intersections,i->member(l,i)));
    worldPointsFF := random(FF^4,FF^0);
    scan(#intersections,i->(
	    pointi := if member(i,depPoints) then (
		li := first select(1,pointsOnLineIndices,l->member(i,l));
		<< li << endl;
		a := random FF;
		a*worldPointsFF_{li#0}+(1-a)*worldPointsFF_{li#1}
		) else random(FF^4,FF^1);
	    worldPointsFF = worldPointsFF | pointi
	    )
	);
    worldPoints := sub(worldPointsFF,CC);    
    helperPoints := apply(pointsOnLineIndices, pp->sub(random(FF^4,FF^(2-min(2,#pp))),CC));
    -- future (line below): may be interesting to sample space of variables differently depending on the field we fabricate data over
    sampleCameraParameters := for i from 1 to nvars list sub(random FF,CC);
    subTable := apply(sampleCameraParameters, cameraVars, (a,b) -> b=>inputGate a);
    sampleC := apply(C,cam -> (
	    M := mutableMatrix(CC, 3, 4);
	    evaluate(cam, mutableMatrix{sampleCameraParameters}, M);
	    matrix M
	    )
	    );
    (
	sampleCameraParameters,
	apply(sampleC, cam->(
	    	P := cam * worldPoints;
	    	L := matrix apply(nLines+nGhosts, l->(
		      Hl := helperPoints#l;
		      Pl := P_(pointsOnLineIndices#l);
		      Pl = Pl|cam*Hl;
	    	      line := randomLineThroughPoints(Pl, FF);
		      {(1/norm(2,line))*line} -- normalization that seemed to benefit the chicago solver
		       ));
	       (P,L)
	       )),
       sampleC
	) 
    )

-- functions which fabricate input of the form (parameter, solution)
-- notation for supp materials is (y,c)---previous notation is (p,x)
encodey = (P,L,projs,FF) -> (
    c := transpose matrix{P};
    allLines := L/last;
    yIndLineMatrix := cfold(
	    allLines/(m->m^(toList indepLines))
	    );
    yInd := matrix(yIndLineMatrix,product toList size yIndLineMatrix,1);
    yDep := if (#depLines == 0) then random(CC^0,CC^1) else 
    rfold(
	depLines/(l-> (
	    	lineInds := first select(1,last D,i->member(l,i));
		triplet := take(lineInds, 2) | {l};
	    	rfold(allLines/(m -> (
	    		n :=numericalKernel(transpose m^triplet, kTol);
	    		(1/n_(2,0))*n^{0,1}
			))
	    	    )
		)
    	    )
	);
    ytChart := sub(randKernel(transpose(c^{4*(m-1)..(4*(m-1)+3*(m-1)-1)}||matrix{{1_CC}}), FF),CC);
    yqChart := sub(rfold(
	for i from 0 to m-2 list randKernel(transpose(c^{4*i..4*i+3}||matrix{{1_CC}}), FF)
	),CC);
    yInd||yDep||ytChart||yqChart
    )

encodeyc = (P, L, projs,FF) -> (
    c := transpose matrix{P};
    y := encodey(P, L, projs,FF);
    (point y, point c)    
    )    

fabricateyc = FF -> (
    (P, L, projs) := fabricatePair(D, FF, nvars); -- these variable names are confusing
    encodeyc(P, L, projs,FF)
    )

-- convenience functions for minors
minors (GateMatrix, ZZ, Sequence, Boolean) := o ->  (M, k, S, laplace) -> (
    (Sm, Sn) := (first S, last S);
    (m,n) := (numrows M, numcols M);
    assert(k<=min(m,n));
    assert(all(Sm,s->#s==k));
    assert(all(Sn,s->#s==k));
    flatten apply(Sm,sm->apply(Sn, sn -> 
	    if (laplace) then laplaceDet submatrix(M,sm,sn)
	    else det submatrix(M,sm,sn)
	    ))
    )

allMinors = method(Options=>{Laplace=>false})
allMinors (GateMatrix, ZZ) := o -> (M, k) -> (
    (m, n ) := (numrows M, numcols M);
    s := (subsets(0..m-1,k),subsets(0..n-1,k));
    minors(M, k, s, o.Laplace)
    )

maxMinors = method(Options=>{Laplace=>false})
maxMinors GateMatrix := o -> M -> allMinors(M,min(numrows M, numcols M), Laplace=>o.Laplace)

-- this seems to work
complexQR = M -> (
    A := mutableMatrix M;
    k := ring A;
    Q := mutableMatrix(k,0,0,Dense=>true);
    R := mutableMatrix(k,0,0,Dense=>true);
    rawQR(raw A, raw Q, raw R, true);
    assert(areEqual(Q*R,A)); -- idk if it will work every time!
    (matrix Q,matrix R)
    )

leverageScores = M -> (
    Q = first complexQR M;
    rsort apply(numrows Q,i->(norm(2,Q^{i}),i))
    )



-- indexes subsystem giving Jacobian rank
rowSelector = method(Options=>{Threshold=>1e-4})
rowSelector Matrix := o -> J0 -> (
    (m , n) := size J0;
    J0' := J0^{0}; -- should technically check this is not identically zero
    i := 1;
    inds := new MutableList from {0};
    k := numericalRank J0';
    while ( (i<m) and (k < n)) do (
	if (numericalRank(J0'||J0^{i},Threshold=>o.Threshold) > k) then (
	    inds#k = i;
	    k = k+1;
	    J0' = J0'||J0^{i};
	    );
	i = i + 1;
	);
    if k < n then  << "WARNING: Jacobian has rank" << toString k << endl;
    toList inds
    )

leverageScoreRowSelector = J0 -> (
    sortedRows := (leverageScores J0)/last;
    r := rowSelector J0^(sortedRows);
    sort(r/(i->sortedRows#i))
    )

log10 = x -> log(x)/log(10)

argCC = z -> atan((imaginaryPart z)/(realPart z))

-- complex number whose real and imag parts are standard normal
gaussCC = () -> (
    (u1,u2):=(random RR,random RR);
    sqrt(-2*log(u1))*cos(2*pi*u2)+ii*sqrt(-2*log(u1))*sin(2*pi*u2)
    )

-- random sample drawn from normal distriution N(mu, var^2)
rNorm = (mu,var) -> mu+var*(realPart gaussCC())_CC

-- random sample from (n-1)-sphere with radius r
sphere = (n,r) -> (
    l:=apply(n,i->rNorm(0,1));
    matrix{r/norm(2,l)*l}
    )

-- assumes "u" of unit length
householder=method()
householder (Matrix,ZZ) := (u,n) -> (
    if (numrows u > 1) then error("householder takes a row vector");
    R:=ring u;
    k:=numcols u;
    id_(R^(n-k))++(id_(R^k)-2*(transpose u)*u)
    )
householder (List,ZZ) := (u,n) -> householder(matrix {u},n)

randomOn = n -> diagonalMatrix(toList((n-1):1_RR)|{(-1)^(random 2)}) * fold(reverse apply(2..n,i->householder(sphere(i,1),n)),(a,b)->a*b)

randomCameraNormalized = () -> (
    R := randomOn 3;
    t := matrix{{random CC},{random CC},{random CC}};
--    t := transpose matrix{sphere(3,1)};
    tnorm := (1 / t_(2,0))* t;
    (R|tnorm)
    )

randomCameraNormalizedCayley = () -> (
    R := cay2R(random CC, random CC, random CC,Normalized=>true);
    t := matrix{{random CC},{random CC},{random CC}};
--    t := transpose matrix{sphere(3,1)};
    tnorm := (1 / t_(2,0))* t;
    (R|tnorm)
    )


randomCamera = () -> (
    R := randomOn 3;
    t := transpose matrix{sphere(3,1)};
    (R|t)
    )

ranks = method(Options=>{})
ranks (Matrix, Matrix) := o -> (x, p) -> (
    if (numcols x > 1) then x = matrix(x, numcols x,1);
    if (numcols p > 1) then p = matrix(p, numcols p,1);
    a := PE/( m -> (
	    evaluate(first m, mutableMatrix(x||p), last m);
	    numericalRank matrix last m
	    )
	    );
    b := LE/( m -> (
	    evaluate(first m, mutableMatrix(x||p), last m);
	    numericalRank matrix last m
	    )
	    );
   (a, b)
   )
ranks (Point,Point) := o -> (x,p) -> ranks(matrix x, matrix p)

rankCheck = method(Options=>{Hard=>true})
rankCheck (Matrix, Matrix) := o -> (x, p) -> (
   (a, b) := ranks(x, p);
   if (o.Hard) then (all(a,x->x==3) and all(b,x->x==2))
     else (all(a,x->x<=3) and all(b,x->x<=2))
   )

cpMatrix = t -> matrix{{0,-t_(2,0),t_(1,0)},{t_(2,0),0,-t_(0,0)},{-t_(1,0),t_(0,0),0}}

essential = (R,t) -> R * cpMatrix t

pCompose = method()
pCompose (MutableHashTable, MutableHashTable) := (H1, H2) -> (
    new MutableHashTable from apply(keys H2,k-> if H1#?(H2#k) then k=> H1#(H2#k))
    )
///TEST
H1 = new MutableHashTable from {0=>1, 1=>2, 2=>0}
H2 = new MutableHashTable from {0=>2, 1=>1, 2=>0}
H1H2=pCompose(H1,H2)
H2H1=pCompose(H2,H1)
assert(H1H2#1 == 2)
assert(H2H1#1 == 0)
///

inverse MutableHashTable := H -> new MutableHashTable from apply(keys H, values H, (k,v) -> v=>k)

-*
G=V.Graph
unvisited = set toList G.Vertices
unvisited = unvisited - {V}
spanningEdges = set{}
cycleMakingEdges = set{}
while #unvisited > 0 do (
    e := first keys unvisited;
    
    
apply(verts, v -> 
*-

writePermutations = method(Options=>{})
writePermutations (List, String) := o -> (L, filename) -> (
    goodPerms := select(L,p->all(p,pi->instance(pi,ZZ)));
    perms := goodPerms/(P->P/(i->i+1)); -- increment letters by 1 for GAP
    file := openOut (currentFileDirectory | filename);
    for i from 0 to #perms-1 do file << "p" << i << ":= PermList(" << toString(new Array from perms#i) << ");" << endl;
    file << "G:=Group(";
    for i from 0 to #perms-2 do file << "p" << i << ", ";
    file << "p" << #perms-1 << ");";
    close file;
    )
-- todo: preselect a non-bottleneck vertex
writePermutations (HomotopyNode, ZZ, String) := o -> (V, rc, filename) -> (
    if rc != length V.PartialSols then (f:=openOut filename; f << "failed"; close f) else (
        idPerm := new MutableHashTable from for i from 0 to rc-1 list i => i;
        -- step 0: extract subgraph of complete correspondences
        G := V.Graph;
        -- EG = edges with complete correspondence
        EG := reverse toList select(G.Edges, e -> (
                ve1 := delete(,unique values e.Correspondence12); 
                ve2 := delete(,values e.Correspondence21); 
                rc == length ve1 and rc == length ve2 
                )
            );
        -- VG = connected component of V in subgraph of EG
        VG := set(apply(EG, e -> e.Node1) | apply(EG, e -> e.Node2));
        neighbor = (v, e) -> if v === e.Node1 then e.Node2 else if v === e.Node2 then e.Node1 else error "edge not incident at vertex";
        -- STEP 1: BUILD SPANNING TREE on (VG, EG)
        (uncoveredV, uncoveredE) := (VG-set{V},set EG);
        T := new MutableHashTable from {};
        while (#uncoveredV > 0) do (
            -- select an uncovered vertex adjacent to a covered vertex
            vList := select(1, keys uncoveredV, v -> any(v.Edges, e -> (u := neighbor(v, e); not member(u, uncoveredV) and member(e, uncoveredE))));
            if # vList == 1 then (
                v := first vList;
                -- select an edge w/ complete correspondence
                eList := select(1, reverse v.Edges, e -> (u := neighbor(v, e); not member(u, uncoveredV) and member(e, uncoveredE)));
                if #eList == 1 then (
                    T#v = first eList; 
                    uncoveredE = uncoveredE - set{e};
                    );
                );
            uncoveredV = uncoveredV - set{v};
            );
        -- STEP 2: extract permutations from cycle basis
        perms := values \ apply(keys uncoveredE, e -> (
                (u, v) := (e.Node1, e.Node2);
                uPath := idPerm;
                while T#?u do (
                    ei = T#u;
                    uPath = if u === ei.Node1 then pCompose(ei.Correspondence12, uPath) else pCompose(ei.Correspondence21, uPath);
                    u = neighbor(u, T#u);
                    );
                vPath := idPerm;
                while T#?v do (
                    ei = T#v;
                    vPath = if v === ei.Node1 then pCompose(ei.Correspondence12, vPath) else pCompose(ei.Correspondence21, vPath);
                    v = neighbor(v, T#v);
                    );
                pCompose(vPath, pCompose(e.Correspondence12, inverse uPath))
                )
            );
        writePermutations(perms,filename);
        );
    )

-- "join" of two GateSystems (take all functions from both)
GateSystem || GateSystem := (P, Q) -> (
    allVars := unique( (flatten entries vars P) | (flatten entries vars Q) );
    allParams := unique( (flatten entries parameters P) | (flatten entries parameters Q) );
    gateSystem(
	gateMatrix{allParams},
	gateMatrix{allVars},
    	(gateMatrix P)||(gateMatrix Q)
	)
    )

-- sum of two GateSystems
GateSystem + GateSystem := (P, Q) -> (
    if (numFunctions P =!= numFunctions Q) then error "can only add GateSystems of the same shape";
    H := P || Q;
    gateSystem(parameters H, vars H, gateMatrix P + gateMatrix Q)
    )

-- take some functions from the GateSystem
GateSystem ^ List := (P, inds) -> gateSystem(parameters P, vars P, (gateMatrix P)^inds)

evaluateJacobian (Point, Point, GateSystem) := (y0, c0, F) -> (
    J := diff(vars F, gateMatrix F);
    (M, N) := (numrows J, numcols J);
    JGS := gateSystem(parameters F, vars F, transpose matrix{flatten entries J});
    matrix(transpose evaluate(JGS, y0, c0),M,N)
    )

-- helpers for rowSelector
orthoProjectQR = (M,L) -> (
    (Q,R)=complexQR L;
--    Q := (SVD L)#1;
    Q*conjugate transpose Q *M
    )

-- orthonormal basis for col(L) using SVD
ONB = L -> (
    (S,U,Vt) := SVD L;
    r := # select(S,s->not areEqual(s,0));
    U_{0..r-1}
    )

-- component of col(L) that is perpendicular to M
perp = method(Options=>{UseSVD=>true})
perp (Matrix, Matrix) := o -> (M, L) -> if areEqual(norm L, 0) then M else (
    Lortho := if o.UseSVD then ONB L else first complexQR L; -- QR seems buggy
    Lperp := M-Lortho*conjugate transpose Lortho * M;
    if o.UseSVD then ONB Lperp else first complexQR Lperp
    )

-- extract a LocalRegularSequence from F near (parameter, solution) = (y0, c0)
-- of length <= n
rowSelector = method(Options=>{BlockSize=>1,UseSVD=>true,Verbose=>false})
rowSelector (Point, Point, ZZ, GateSystem) := o -> (y0, c0, n, F) -> (
    blockSize := o.BlockSize;
    numBlocks = ceiling((numFunctions F)/blockSize);
    numIters=0;
    L := matrix{for i from 1 to 14 list 0_CC}; -- initial "basis" for row space
    r := 0;
    goodRows := {};
    diffIndices := {};
    while (r < n and numIters < numBlocks) do (
    	diffIndices = for j from numIters*blockSize to min((numIters+1)*blockSize,numFunctions F)-1 list j;
	if o.Verbose then << "processing rows " << first diffIndices << " thru " << last diffIndices << endl;
    	newRows := evaluateJacobian(y0,c0,F^diffIndices);
    	for j from 0 to numrows newRows - 1 do (
	    tmp := transpose perp(transpose newRows^{j}, transpose L);
	    if not areEqual(0, norm tmp) then (
		if o.Verbose then << "added row " << blockSize*numIters+j << endl;
	    	if areEqual(norm L^{0}, 0) then L = tmp else L = L || tmp;
	    	goodRows = append(goodRows, blockSize*numIters+j);
		);
    	    );
    	r = numericalRank L;
    	numIters = numIters+1;
	);
    if o.Verbose then << "the rows selected are " << goodRows << endl;
    goodRows
    )

squareDown = method(Options=>{BlockSize=>1, Verbose=>false})
squareDown (Point, Point, ZZ, GateSystem) := o -> (y0, c0, n, F) -> F^(rowSelector(y0, c0, n, F, BlockSize => o.BlockSize, Verbose=>o.Verbose))
