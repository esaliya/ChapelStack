use BlockDist,ReplicatedDist;
use Random;
config var numDim: int = 2,
	   numClusters: int = 4,
	   numPoints: int  = 2000,
	   numIter: int = 50,
	   threshold: real = 0.0001;

proc main() {
  writeln ("dimensions: ", numDim);
  writeln ("clusters ", numClusters);
  writeln ("points: ", numPoints);
  writeln ("iterations: ", numIter);
  writeln ("threshold: ", threshold);
  
  /* domains */
  const PointsSpace = {0..#numPoints,0..#numDim};
  const ClusterSpace = {0..#numClusters,0..#numDim};
  const ClusNumSpace = {0..#numClusters};
  /* distribute points and initialize to random values */
  var blockedLocaleView = {0..#numLocales,1..1};
  var blockedLocales: [blockedLocaleView] locale 
    = reshape(Locales, blockedLocaleView);
  // this in essence blocks along the 
  // zeroth dimension (i.e. range 0..#numPoints)
  const BlockedPointsSpace = PointsSpace 
    dmapped Block(boundingBox=PointsSpace,
    targetLocales=blockedLocales);
  var points: [BlockedPointsSpace] real;
  var rand: RandomStream = new RandomStream();
  rand.fillRandom(points);

  /* current centers initilized to 
     first N points where N = numClusters */
  var currCenters: [ClusterSpace] real;
  currCenters = points[ClusterSpace];
  
  /* new centers */
  var newCenters: [ClusterSpace] real;

  /* replicated local center */
  const ReplClusterSpace = ClusterSpace 
    dmapped ReplicatedDist();
  var localCurrCenters: [ReplClusterSpace] real;
  // work around to non existing atomic blocks
  var localCenterUpdates: [ReplClusterSpace] atomic real;
  const ReplClusNumSpace = ClusNumSpace dmapped ReplicatedDist();
  var localCenterPCounts: [ReplClusNumSpace] atomic int; 
 
  var converged: bool = false;
  var step: int = 0;
  // whish something like real.max exist
  const MaxDist = 9223372036854775807;
  while (step < numIter && !converged) {
    // copy currCenters to each replicand
    localCurrCenters = currCenters;
    // reset local updates and counts
    forall lcu in localCenterUpdates do lcu.write(0.0);
    forall lcpc in localCenterPCounts do lcpc.write(0);

    forall p in {0..#numPoints} {
      on points[p,0] {
        var closestC: int = -1;
	var closestDist: real = MaxDist;
	for c in {0..#numClusters} { // possibility to parallelize
	  var dist: atomic real;
	  dist.write(0.0);
	  forall d in {0..#numDim} {
	    var tmp = points[p,d] - localCurrCenters[c,d];
	    dist.add(tmp*tmp);
	  }

	  if (dist.read() < closestDist) {
	    closestDist = dist.read();
	    closestC = c; 
	  }
	}
        
	forall d in {0..#numDim} {
	  atomic { // would have been great, but not implemented yet in Chapel
	    localCenterUpdates[closestC,d].add(points[p,d]);
	  }
	}
	localCenterPCounts[closestC].add(1);
      }
    }
    
    var tmpLCU: [ClusterSpace] atomic real;
    forall tlcu in tmpLCU do tlcu.write(0.0);
    var tmpLCPC: [0..#numClusters] atomic int;
    forall tlcpc in tmpLCPC do tlcpc.write(0);

    forall loc in Locales {
      on loc do {
        cobegin {
	  forall (i,j) in ClusterSpace {
	    tmpLCU[i,j].add(localCenterUpdates[i,j].read());
	  }
	  forall i in {0..#numClusters} {
	    tmpLCPC[i].add(localCenterPCounts[i].read());
	  }
	}
      }
    }

    var b: atomic bool;
    forall (i,j) in ClusterSpace {
      newCenters[i,j] = tmpLCU[i,j].read()/tmpLCPC[i].read();
      if (abs(newCenters[i,j] - currCenters[i,j]) > threshold){
        b.write(false);
      }
    }
    converged = b.read();
    step += 1;
  }

  writeln("new centers");
  writeln(newCenters);

  writeln("currCenters");
  writeln(currCenters);
}
