use BlockCycDist;
use Random;
config var numDim = 2,
	   numClusters = 4;

proc main() {
  writeln ("dimensions: ", numDim);
  writeln ("clusters ", numClusters);
  
  /* domains */
  const ClusterSpace = {0..numClusters-1,0..numDim-1};
  const CenterSpace = {0..numLocales-1,0..numClusters-1,0..numDim-1};
  
  /* arrays to hold current and new cluster centers */
  var currCenters : [ClusterSpace] real;
  for (i,j) in ClusterSpace do {
    currCenters[i,j] = i+0.1*j;
  }
  writeln("curr centers");
  writeln(currCenters);
  
  var centerLocales : [CenterSpace] locale;
  forall ijk in CenterSpace do {
    centerLocales[ijk] = Locales[ijk(1)];
  }

  // block cyclic domain for local center arrays
  const BCCenterSpace = CenterSpace dmapped BlockCyclic(
                          startIdx=CenterSpace.low, blocksize=(1,numClusters,numDim),targetLocales=centerLocales);
  // each locale will get a copy of the current centers
  var localCurrCenters : [BCCenterSpace] real;
  for e in localCurrCenters do {
    e = e.locale.id;
    writeln("lcc");
    writeln(localCurrCenters);
    writeln();
  }



/*  forall (i,j,k) in CenterSpace do{ 
    localCurrCenters[i,j,k] = ;
    writeln("(", i, ",", j, ",", k, ") is in", localCurrCenters[i,j,k].locale, " and its value is: ", localCurrCenters[i,j,k]);
  }*/

  

  /* Test code */
  /*forall p in points do p=here.id;*/

  /*forall c in localCenters do {
    c = here.id;
  }*/

  /*writeln("curr centers");
  writeln(currCenters);*/

/*  writeln("Local curr centers");
  writeln(localCurrCenters);*/
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

  

}
