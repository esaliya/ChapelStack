/* Author: Saliya Ekanayake 
           esaliya at gmail dot com */

use BlockDist,ReplicatedDist;
use Random;
config var numDim: int = 3 ,
	   numClusters: int = 4,
	   numPoints: int  = 20,
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
  const ClusterNumSpace = {0..#numClusters};
  /* distribute points and initialize to random values */
  var blockedLocaleView = {0..#numLocales,1..1};
  var blockedLocales: [blockedLocaleView] locale 
    = reshape(Locales, blockedLocaleView);
  // this in essence blocks along the 
  // zeroth dimension (i.e. range 0..#numPoints)
  const BlockedPointsSpace = PointsSpace 
    dmapped Block(boundingBox=PointsSpace,
    targetLocales=blockedLocales);
    var points: [BlockedPointsSpace] real = ((0.687341, 0.29802, 0.752728),
(0.419732, 0.37836, 0.0912831),
(0.248889, 0.0659724, 0.416977),
(0.711823, 0.0305533, 0.280535),
(0.52954, 0.818765, 0.224503),
(0.981434, 0.225217, 0.944875),
(0.383762, 0.150229, 0.318658),
(0.453673, 0.880605, 0.361181),
(0.90988, 0.276158, 0.974794),
(0.254215, 0.0159331, 0.792719),
(0.796911, 0.24941, 0.759422),
(0.513687, 0.243573, 0.625268),
(0.857396, 0.811445, 0.20519),
(0.539379, 0.386104, 0.446041),
(0.778445, 0.611402, 0.438312),
(0.333145, 0.0453983, 0.0177009),
(0.105817, 0.00589306, 0.269783),
(0.117595, 0.726253, 0.348978),
(0.359102, 0.646129, 0.795696),
(0.0122162, 0.309424, 0.420688));
  
  /* current centers initilized to 
     first N points where N = numClusters */
  var currCenters: [ClusterSpace] real;
  currCenters = points[ClusterSpace];
  writeln("initial centers");
  writeln(currCenters);

  /* replicated local center */
  const ReplClusterSpace = ClusterSpace 
    dmapped ReplicatedDist();
  var localCurrCenters: [ReplClusterSpace] real;
  // work around to non existing atomic blocks
  var localCenterUpdates: [ReplClusterSpace] atomic real;
  const ReplClusterNumSpace = ClusterNumSpace dmapped ReplicatedDist();
  var localCenterPCounts: [ReplClusterNumSpace] atomic int; 
 
  var converged: bool = false;
  var step: int = 0;
  // whish something like real.max exist
  const MaxDist = 9223372036854775807;
  while (step < numIter && !converged) {
    cobegin {
      // copy currCenters to each replicand
      localCurrCenters = currCenters;
      // reset local updates and counts
      forall lcu in localCenterUpdates do lcu.write(0.0);
      forall lcpc in localCenterPCounts do lcpc.write(0);
    }

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
    var tmpLCPC: [ClusterNumSpace] atomic int;
    forall tlcpc in tmpLCPC do tlcpc.write(0);

    forall loc in Locales {
      on loc do {
        cobegin {
	  forall (i,j) in ClusterSpace {
	    tmpLCU[i,j].add(localCenterUpdates[i,j].read());
	  }
	  forall i in ClusterNumSpace {
	    tmpLCPC[i].add(localCenterPCounts[i].read());
	  }
	}
      }
    }

    var b: atomic bool;
	b.write(true);
    forall (i,j) in ClusterSpace {
      var center: real = tmpLCU[i,j].read()/tmpLCPC[i].read();
      if (abs(center - currCenters[i,j]) > threshold){
        b.write(false);
      }
      currCenters[i,j] = center;
    }
    converged = b.read();
    step += 1;
  }

  writeln("final centers");
  writeln(currCenters);
}
