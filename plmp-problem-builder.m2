-- IMPORTS
needs "common.m2"

-- GLOBALS
FF = CC
kTol=1e-4

nLines = D#0 + D#1
depLines = flatten((last D)/(i -> drop(i,2)))
indepLines = sort toList(set(0..nLines-1)-set(depLines))


-*
-- SET TRACKER OPTIONS HERE
-- null indicates default value 
scan({CorrectorTolerance=>null,
	EndZoneFactor=>null,
	InfinityThreshold => null, 
	maxCorrSteps => null, 
	NoOutput => null,
	numberSuccessesBeforeIncrease => null,
	Precision => null,
	Predictor => null,
	stepIncreaseFactor => null,
	tStep => null,
	tStepMin => null
	}, 
    opt -> setDefault(opt))
--view current tracker settings
netList apply({CorrectorTolerance,
	EndZoneFactor,
	InfinityThreshold,
	maxCorrSteps,
	NoOutput,
	numberSuccessesBeforeIncrease,
	Precision,
	Predictor,
	stepIncreaseFactor,
	tStep,
	tStepMin
	}, 
    opt -> {opt, getDefault(opt)})
*-

-- setting up input gates
-- "l" indices go (view, feature, coordinate)
varInSymbs = (flatten (for i from 2 to m list {w_i,x_i,y_i,z_i}))|
    (flatten (for i from 2 to m list {t_(i,1),t_(i,2),t_(i,3)}))

nvars = #varInSymbs -- should change variable name
inSymbs=varInSymbs|
    (flatten for j in indepLines list flatten for i from 1 to m list for k from 1 to 3 list l_(i,j,k))|
	(flatten for j in depLines list flatten for i from 1 to m list for k from 1 to 2 list a_(i,j,k))|
	(flatten (for i from 2 to m list(
		     {ct_(i,1), ct_(i,2), ct_(i,3)}
		     )))|{ct_0}|
	 (flatten (for i from 2 to m list toList(cq_(i,0)..cq_(i,4))))
scan(inSymbs,s->value(toString(s)|" = inputGate " | toString(s)))
inGates = toList(inSymbs)/value

-- variable groups
cameraVars = take(inGates,nvars)
centerVars = drop(cameraVars,4*(m-1))
qVars = take(cameraVars,4*(m-1))
varMatrix = matrix{cameraVars}

-- parameter groups
dataParams = drop(inGates,nvars)
chartParams = drop(dataParams,
    #dataParams-(3*(m-1)+1+5*(m-1)))
tchartParams = take(chartParams,3*(m-1)+1)
qchartParams = drop(chartParams,3*(m-1)+1)
paramMatrix = matrix{dataParams}
inputMatrix = varMatrix | paramMatrix

-- chart equations
tChart = matrix{tchartParams}*transpose (matrix{centerVars}|matrix{{1_CC}})
qCharts = rfold for i from 0 to m-2 list (
	matrix{take(qchartParams,{5*i,5*i+4})}*
    transpose (matrix{take(qVars,{4*i,4*i+3})}|matrix{{1_CC}})
    )
charts = tChart||qCharts

-- rotation and projection matrices
R = {gateMatrix(id_(FF^3))} | for i from 2 to m list Q2R(w_i,x_i,y_i,z_i)
T = {gateMatrix{{0},{0},{0}}} | for i from 2 to m list (w_i^2+x_i^2+y_i^2+z_i^2)*transpose matrix{{t_(i,1),t_(i,2),t_(i,3)}}
pMatrices = apply(R,T,(r,t)->r|t)

-- camera evaluator used for fabrication routine
C = pMatrices/(P->makeSLProgram(varMatrix, P))

--row vectors giving implicit equations of 2d visibleLines
-- better variable name?
visibleLines = for i from 1 to m list for j from 0 to nLines-1 list (
    if member(j,indepLines) then matrix{{l_(i,j,1),l_(i,j,2),l_(i,j,3)}} 
    else (
	whichIncs = select(last D, i -> member(j,i));
	assert(#whichIncs == 1);
	(k1,k2) = (whichIncs#0#0, whichIncs#0#1);
	a_(i,j,1)*matrix{{l_(i,k1,1),l_(i,k1,2),l_(i,k1,3)}} + a_(i,j,2)*matrix{{l_(i,k2,1),l_(i,k2,2),l_(i,k2,3)}}
	)
    )

-- backprojected planes of visiblelines for each view
planes = for i from 0 to m-1 list apply(visibleLines#i,l->l*pMatrices#i)

-- matrices whose rank deficicency encodes a common point on world lines
CoPmatrices = (last D)/(i->
    rfold(i/(l->
	rfold(planes/(p->
	    p#l
	    ))
	))
    )

-- matrices whose rank deficicency encodes a common world line
-- not needed for m==2
CoLmatrices = if (m<3) then {} else apply(D#0,l->
    rfold(planes/(p->
	    p#l
	    ))
    )

-- master system w/ all max minors from CoL, CoP matrices
elapsedTime F=charts || transpose gateMatrix{
    flatten(
	CoLmatrices/(M -> allMinors(M, 3, Laplace=>COFACTOR)) |
	CoPmatrices/(M -> maxMinors(M, Laplace=>COFACTOR))
	)
    };
<< " number of polynomials is " << numrows F << endl
-- filter path jumps during monodromy
filterEval = (p,x) -> (
    -- false iff residual small
    resid := norm evaluate(F,x||p);
--    << "residual: " << resid << endl;
    (resid > 4e-4)
    )

-- matrix evaluators: used for rank filter
PE=apply(#CoPmatrices,i->(
	G := CoPmatrices#i;
	(makeSLProgram(inputMatrix,G),
	mutableMatrix(FF,numrows G,numcols G)
	)
    )
)
LE=apply(#CoLmatrices,i->(
	G := CoLmatrices#i;
	(makeSLProgram(inputMatrix,G),
	mutableMatrix(FF,numrows G,numcols G)
	)
    )
)


filterRank = (p,x) -> (
    -- false iff residual small
    (cop,col) := ranks(matrix x,matrix p);
    not(all(cop,x->x==3) and all(col,x->x==2))
    )

filterRankCoP = (p,x) -> (
    -- false iff residual small
    a := PE/( m -> (
	    evaluate(first m, mutableMatrix(x||matrix p), last m);
	    numericalRank matrix last m
	    )
	    );
    not all(a,x->x==3)
    )

--setRandomSeed 31452345342
(y, c) = fabricateyc CC
filterRank(y,c)
varMatrix = gateMatrix{cameraVars}
paramMatrix = gateMatrix{dataParams}
masterGS = gateSystem(paramMatrix, varMatrix, F);
norm evaluate(masterGS,y,c) -- ~0?

-*
-- this block is very hacky!
if (instance(Jpivots, Symbol) and JACOBIAN) then (
    -- better to have this precomputed
    << "differentiating" << endl;
    elapsedTime J = diff(varMatrix,F);
    (M,N) = size J;
    elapsedTime JGS = gateSystem(paramMatrix, varMatrix, transpose matrix{flatten entries J});
    elapsedTime J0 = matrix(transpose evaluate(JGS,y,c),M,N);
    elapsedTime Jpivots = rowSelector(J0,Threshold=>1e-6);
    elapsedTime S = first SVD J0^Jpivots;
    << "pivot indices are " << toString Jpivots << endl;
    )
*-
elapsedTime GS= squareDown(y, c, numVariables masterGS, masterGS)

(y, c) = fabricateyc CC
filterRank(gammify y,c)
if not instance(NEDGES,ZZ) then NEDGES=4
if not instance(NNODES,ZZ) then NNODES=2
if not instance(ROOTCOUNT,ZZ) then ROOTCOUNT=null
if instance(SATURATE,Symbol) then SATURATE=true
monodromyGroup(GS, y, {c},
    "msOptions" => {
	Verbose=>true,
    	FilterCondition=>filterRank,
    	Randomizer=>gammify,
    	EdgesSaturated => SATURATE,
    	NumberOfEdges=>NEDGES,
    	NumberOfNodes=>NNODES,
    	TargetSolutionCount=>ROOTCOUNT
	},
    FileName => FILENAME
    );
-- clear symbols for next run
w=symbol w
x=symbol x
y=symbol y
z=symbol z
t=symbol t
l=symbol l
a=symbol a
ct=symbol ct
cq=symbol cq
Jpivots = symbol Jpivots
end
