module GRAPH500_Kronecker_Generator
{
  /*
   *  Kronecker Edge List Generator for Graph500 Kernel 1
   *  --------------
   *  Generate an edgelist according to the Graph500
   *  parameters.  In this sample, the edge list is
   *  returned in an array with two rows, where StartVertex
   *  is first row and EndVertex is the second.  The vertex
   *  labels start at zero.
   */

  proc generate_edge_list(SCALE : int(64), edge_factor : int(64), maxWeight: int) {
    use Random;
    use BlockDist;

    writeln("Generating edge list ... ");
    
    write("  Creating initial data structures ... ");

    // Number of vertices
    const N: int(64) = 2 ** SCALE;

    // Create initial tree edges
    var halfN: int(64) = N/2;
    var r1: range(int(64), BoundedRangeType.bounded, stridable=true);
    var r2: range(int(64), BoundedRangeType.bounded, stridable=true);
    r1 = 1..N;
    r2 = 1..2;
    const ijTreeD: domain(2) dmapped Block(boundingBox=[1..2,1..N]) = [r2,r1];
    const k: [1..halfN] int(64);
    const ijTree: [ijTreeD] int(64);

    // Number of edges
    const M: int(64) = edge_factor * N;

    // Set initiator probabilities
    const A: real = 0.55;
    const B: real = 0.1;
    const C: real = 0.1;

    // Reusable domains
    const DM: domain(1) dmapped Block(boundingBox=[1..M]) = [1..M];

    // Storage of intermediate values
    var ab: real = 0.0;
    var c_norm: real = 0.0;
    var a_norm: real = 0.0;

    var edgeDomain: domain(2) dmapped Block(boundingBox=[1..2,1..M]) = [1..2, 1..M];
    var edgeArray: [edgeDomain] int(64);

    var mPlusn: int(64) = M + N;
    var fEdgeDomain: domain(2)dmapped Block(boundingBox=[1..2,1..mPlusn]) = [1..2, 1..mPlusn];
    var fEdgeArray: [fEdgeDomain] int(64);

    writeln("done.");

    write("  Initializing k and ijTree ... ");
    // Init k and ijTree
    [i in [1..halfN]] k[i] = i;

   ijTree[1,1..halfN] = k;
   ijTree[1,(halfN+1)..N] = 2*k;
   ijTree[2,1..halfN] = k;
   ijTree[2,(halfN+1)..N] = 2*k+1;
  
  /*
    forall (i, j) in ijTreeD {
      if (i == 0) {
        if ( j <= N/2 ) {
          ijTree[i, j] = k[j];
        } else {
          ijTree[i, j] = 2 * k[j - N/2];
        }
      } else {
        if ( j <= N/2 ) {
          ijTree[i, j] = k[j];
        } else {
          ijTree[i, j] = 2 * k[j - N/2] + 1;
        }
      }
    }
    */
    
    writeln("done.");

    // Initialize all the elements of edgeArray to 1
    [(i, j) in edgeDomain] edgeArray[i, j] = 1;

    // Loop over each order of bit
    ab = A + B;
    c_norm = C / (1 - (A + B));
    a_norm = A / (A + B);

    write("  starting 'for ib in [1 .. SCALE]' ...");
    forall ib in [1..SCALE] {
      var randArray: [DM] real;
      var ii_bit: [DM] int;
      var jj_bit: [DM] int;
      var ii_jj: [edgeDomain] int;

      fillRandom(randArray);
      [j in DM] ii_bit[j] = randArray[j] > ab;

      fillRandom(randArray);
      [j in DM] jj_bit[j] = if (randArray[j] > c_norm * ii_bit[j] + a_norm*not(ii_bit[j])) then 1 else 0;

      [j in DM] ii_jj[1,j] = ii_bit[j];
      [j in DM] ii_jj[2,j] = jj_bit[j];

      edgeArray = edgeArray + 2**(ib - 1) * ii_jj;
    }
    writeln("done.");

    write("  Starting 'for (i,j) in ijTreeD' ... ");
    forall (i,j) in ijTreeD {
      fEdgeArray[i,j] = ijTree[i,j];
    }
    writeln("done.");

    write("  Starting 'for (i,j) in edgeDomain' ... ");
    forall (i,j) in edgeDomain {
      fEdgeArray[i, N + j] = edgeArray[i,j];
    }
    writeln("done.");

    // Adjust to zero based lables
    fEdgeArray = fEdgeArray - 1;

    var randNums = new RandomStream(parSafe=false);

    var seed0: uint(64) = (randNums.getNext()*N):uint ;
    var seed1: uint(64) = (randNums.getNext()*N):uint;

    write("  Skipping scramble ... ");
    // Scramble vertex labels
    // [(i,j) in fEdgeDomain] fEdgeArray[i,j] = scramble(fEdgeArray[i,j], SCALE, seed0, seed1);
    writeln("done.");

    write("  Creating weight array ... ");
    var weightD: domain(1) dmapped Block(boundingBox=[fEdgeArray.domain.dim(2)]) = [fEdgeArray.domain.dim(2)];
    var weightArray: [weightD] int;
    writeln(".done");

    write("  Filling weight array ... ");
    var tmpWeightArray : [weightD] real;
    fillRandom(tmpWeightArray);
    [j in weightD] weightArray[j] = (tmpWeightArray[j] * 255) : int(64);
    writeln("done.");

    write("  Creating weighted edge list ... ");
    var weightedEdgeListDomain: domain(2) dmapped Block(boundingBox=[1..3,fEdgeArray.domain.dim(2)]) = [1..3, fEdgeArray.domain.dim(2)];
    var weightedEdgeList: [weightedEdgeListDomain] int(64);

    [j in fEdgeArray.domain.dim(2)] weightedEdgeList[1,j] = fEdgeArray[1,j];
    [j in fEdgeArray.domain.dim(2)] weightedEdgeList[2,j] = fEdgeArray[2,j];
    [j in fEdgeArray.domain.dim(2)] weightedEdgeList[3,j] = weightArray[j];
    writeln("done.");

    writeln("  Returning weighted edge list ... done.");
    writeln("done.");
    return weightedEdgeList;
  }

  // SCALE should be at least 26 to get non-zero output from this procedure.
  proc scramble(val: int(64), SCALE: int, seed0: uint(64), seed1: uint(64)){
    var v: uint(64) = val:uint(64);
    v += seed0 + seed1;
    v *= (seed0 | 0x4519840211493211);
    v = (bitreverse(v) >> (64 - SCALE));
    v *= (seed1 | 0x3050852102C843A5);
    v = (bitreverse(v) >> (64 - SCALE));

    return v:int(64);
  }

  proc bitreverse(v: uint(64)): uint(64) {
    var x: uint(64) = v;

    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
    x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
    x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
    x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));

    return((x >> 16) | (x << 16));
  }

  // not operation for integers
  proc not(v: int): int {
    return if (v == 0) then 1 else 0;
  }

  /*proc main(){
    writeln(generate_edge_list(3, 2, 34));
  }*/

}


