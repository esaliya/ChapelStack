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
  const PointsSpace = {0..numPoints-1,0..numDim-1};
  const ClusterSpace = {0..numClusters-1,0..numDim-1};
  const ClusNumSpace = {0..#numClusters};
  /* distribute points and initialize to random values */
  var blockedLocaleView = {0..#numLocales,1..1};
  var blockedLocales: [blockedLocaleView] locale = reshape(Locales, blockedLocaleView);
  // this in essence blocks along the zeroth dimension (i.e. range 0..#numPoints)
  const BlockedPointsSpace = PointsSpace dmapped Block(boundingBox=PointsSpace, 
                                                       targetLocales=blockedLocales);
  var points: [BlockedPointsSpace] real;
  var rand: RandomStream = new RandomStream();
  rand.fillRandom(points);

  /* current centers initilized to first N points where N = numClusters */
  var currCenters: [ClusterSpace] real;
  currCenters = points[ClusterSpace];
  
  /* new centers */
  var newCenters: [ClusterSpace] real;

  /* replicated local center */
  const ReplClusterSpace = ClusterSpace dmapped ReplicatedDist();
  var localCurrCenters: [ReplClusterSpace] real;
  var localCenterUpdates: [ReplClusterSpace] atomic real; // work around to non existing atomic blocks
  const ReplClusNumSpace = ClusNumSpace dmapped ReplicatedDist();
  var localCenterPCounts: [ReplClusNumSpace] atomic int; // work around to non existing atomic blocks
 
  var converged: bool = false;
  var step: int = 0;
  const MaxDist = 9223372036854775807; // whish something like real.max exist
  while (step < numIter && !converged) {
    // copy currCenters to each replicand
    localCurrCenters = currCenters;
    // reset local updates and counts
    forall lcu in localCenterUpdates do lcu.write(0.0);
    forall lcpc in localCenterPCounts do lcpc.write(0);
//    localCenterUpdates = 0.0; compiler says not supported when atomic values are used
//    localCenterPCounts = 0; compiler says not supported when atomic values are used

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

/*    writeln("tmpLCU");
    writeln(tmpLCU);
    writeln("tmpLCPC");
    writeln(tmpLCPC); */

    step += 1;

  }

  writeln("new centers");
  writeln(newCenters);


/* printout */
    /*writeln("lcc");
    writeln(localCurrCenters);
    writeln("lcu");
    writeln(localCenterUpdates);
    writeln("lcpc");
    writeln(localCenterPCounts);
  
  writeln("points");
  writeln(points);*/

  writeln("currCenters");
  writeln(currCenters);

  



 /*  
/* custom locale arrays */
  var pointsLocales : [PointsSpace] locale;
  forall ij in PointsSpace do {
    pointsLocales[ij] = Locales[ij(1)%numLocales];
  }
  var centerLocales : [CenterSpace] locale;
  forall ijk in CenterSpace do {
    centerLocales[ijk] = Locales[ijk(1)];
  }
  
 /* writeln("center locales");
  writeln(centerLocales); */

  /* create points array an initilize to random values */
  // block cyclic domain for points
  var blkSize : int = numPoints / numLocales;
  const BCPointsSpace = PointsSpace dmapped BlockCyclic(
                          startIdx=PointsSpace.low, blocksize=(blkSize,numDim),targetLocales=pointsLocales);
  // points array
  var points : [BCPointsSpace] real;
  // fill with random values
  var rand : RandomStream = new RandomStream();
  rand.fillRandom(points);

  /* arrays to hold current and new cluster centers */
  var currCenters : [ClusterSpace] real;
  // initially old centers are first N points, where N is number of clusters
  currCenters[ClusterSpace] = points[ClusterSpace];
  var newCenters : [ClusterSpace] real;
  

  /* points are distributed across locales. Therefore it's good to have 
   * local arrays containing the current centers and new center updates */

  // block cyclic domain for local center arrays
  const BCCenterSpace = CenterSpace dmapped BlockCyclic(
                          startIdx=CenterSpace.low, blocksize=(1,numClusters,numDim),targetLocales=centerLocales);
  // each locale will get a copy of the current centers
  var localCurrCenters : [BCCenterSpace] real;
  /*forall (i,j,k) in BCCenterSpace do{
    localCurrCenters[i,j,k] = here.id;
    writeln("(", i, ",", j, ",", k, ") is in", localCurrCenters[i,j,k].locale, " and its value is: ", localCurrCenters[i,j,k]);
  }*/

  for l in {0..numLocales-1} do {
//    on loc do {
      for (i,j) in ClusterSpace do {
        localCurrCenters[l, i, j] = currCenters[i,j];
        writeln("lcc in locale: ", l, " after i,j: ", i, " - ", j);
        writeln(localCurrCenters);

      }
  //  }
  }

  // each locale will do updates to separate center arrays
  var localCenterUpdates : [BCCenterSpace] real;
  
  /* iterate and refine current centers */
  var converged : bool = false;
  var step : int = 0;
  
/*
  while (step < numIter && !converged) do {
    writeln("assigning values");
    for i in 0..numLocales-1 do on Locales[i] do {
      for (j,k) in ClusterSpace do {
        localCurrCenters[i,j,k] = currCenters[j,k];
        writeln("on locale: ", i, " with id: ", here.id, " with j: ", j, " k: ", k, " lcc[i,j,k]: ", localCurrCenters[i,j,k], " cc[j,k]: ", currCenters[j,k]);
	// could have done localCenterUpdates = 0.0 outside all loops, but why waste this loop :)
	localCenterUpdates[i,j,k] = 0.0;
      }
    }
    step += 1;
  }
  */

 /* forall ijk in CenterSpace do {
    localCurrCenters[ijk]
  }*/


  

  /* Test code */
  /*forall p in points do p=here.id;*/

  /*forall c in localCenters do {
    c = here.id;
  }*/

  writeln("curr centers");
  writeln(currCenters);

  writeln("Local curr centers");
  writeln(localCurrCenters);
  /*
  for i in {0..numLocales-1} do {
    writeln("locale: ", i);
    for j in {0..numClusters-1} do {
      for k in {0..numDim-1} do {
        write(localCurrCenters[i,j,k], " ");
      }
      writeln();
    }
    writeln();
  }*/

  /*writeln("Local center updates")
  for i in {0..numLocales-1} do {
    writeln("locale: ", i);
    for j in {0..numClusters-1} do {
      for k in {0..numDim-1} do {
        write(localCenterupdates[i,j,k], " ");
      }
      writeln();
    }
    writeln();
  }*/
  /*for i in {0..numLocales-1} do {
     var ele = localClusterCenters[i,0,0];
     writeln(ele.locale.id == i);
  //   writeln(i);
  }*/


//  forall p in points do {
    // an implict locale shift due to the nature of standard distribution's implementaitons.
  //  p = here.id; 
  //}

/*  writeln("points");
  for i in {0..numPoints-1} do {
    for j in {0 .. numDim-1} do {
      write(points[i,j], " ");
    }
    writeln();
  }*/

//  writeln(points);

  
*/
}
