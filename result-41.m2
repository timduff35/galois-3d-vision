restart

-- configure options for MonodromySolver
NNODES=3
NEDGES=3
RUNMONODROMY=true
SATURATE=true
COFACTOR=true
needsPackage "NumericalAlgebraicGeometry"
setDefault(tStepMin=>1e-8)
setDefault(maxCorrSteps=>2)

-- For a description of the problem encoding, see PLMP paper and supplementary material (D., Kohn, Lekyin, P., ICCV 2019)
LIST = {{ -- 2-camera problems    
    (5,0,{{0,1},{1,2},{2,3},{3,4},{4,0}}) => ("5000", 20), -- IMPRIMITIVE: C2 wr S10 ^ A20
    (4,1,{{0,1},{1,2},{2,3},{3,0},{0,4}}) => ("4100_3", 16), -- IMPRIMITIVE: C2 wr S8 ^ A16
    (2,4,{{0,1},{0,2},{0,3},{1,4},{1,5}}) => ("3200_3", 12) -- IMPRIMITIVE
},
{ -- 3-camera problems
    (4,2,{{4,5}}) => ("1040_0", 360), -- full symmetric
    (5,0,{{3,4}}) => ("1032_2", 552), -- full symmetric
    (6,0,{{2,3,4,5}}) => ("1024_4", 480), -- full symmetric
    (4,1,{{2,3},{2,4}}) => ("2021_1", 264), -- full symmetric
    (5,0,{{0,1,2},{0,3}}) => ("2013_2", 432), -- full symmetric
    (5,1,{{0,1,2,3},{0,5}}) => ("2013_3", 328),  -- full symmetric
    (6,0,{{0,1,3,4},{0,2,5}}) => ("2005_3", 480), -- full symmetric
    (6,0,{{0,1,2,3,4},{0,5}}) => ("2005_4", 240), -- full symmetric
    (6,1,{{0,6},{0,1,2,3,4,5}}) => ("2005_5", 64), -- IMPRIMITIVE
    (4,0,{{1,2},{2,3},{1,3}})=> ("3010_0", 216), -- full symmetric
    (5,0,{{0,1,3},{0,2,4},{1,2}}) => ("3002_1", 312), -- full symmetric
    (5,0,{{0,1,3,4},{1,2},{0,2}}) => ("3002_2", 224), -- full symmetric
    (3,2,{{1,2},{1,3},{1,4}})=> ("2111_1", 40), -- full symmetric
    (4,0,{{0,1},{0,2},{0,3}})=> ("2103_1", 144), -- full symmetric
    (4,1,{{0,1,2},{0,3},{0,4}})=>("2103_2", 144), -- full symmetric
    (4,2,{{0,1,2,3},{0,4},{0,5}})=>("2103_3", 144), -- fullsymmetric
    (1,5,{{0,1},{0,2},{0,3},{4,5}})=>("3100_0", 64) -- IMPRIMITIVE: C2 wr (C2 wr S16 ^ A32) ^ A64
    },
{  -- 4-camera problems (w/ fewer than 1K solutions)    
    (3,1,{{0,1},{0,2},{0,3}}) => ("2102_1", 544), -- full symmetric
    (3,2,{{0,1,2},{0,3},{0,4}}) => ("2102_2", 544), -- full symmetric
    (2,3,{{1,2},{1,3},{1,4}}) => ("2110_0", 32) -- IMPRIMITIVE: -*
    } 
}

-- loop computing the Galois/monodromy groups
for i in {0,1,2} do (
    m = i + 2; -- number of cameras
    for P in LIST#i do (
    	print P;
    	D = first P;
    	(file, rc) := (first last P, last last P);
    	DIR = "./groups-relative-pose/" |toString(m) |"-cameras/";
    	FILENAME = DIR|first last P|".txt";
    	if (not fileExists FILENAME) then (
	    << D << endl;
	    --use more nodes for the highly imprimitive 4-camera problem of degree 32 to make Galois group computation more reliable
            if (rc == 32) then NNODES = 4 else NNODES = 3; 
            ROOTCOUNT = rc;
	    setRandomSeed 2022;
    	    compTime = first elapsedTiming load("plmp-problem-builder.m2");
	    stdio << "computation time = " << compTime << "s" << endl;
	    stdio << "------------------------------------------" << endl;
	    );
	);        
)
