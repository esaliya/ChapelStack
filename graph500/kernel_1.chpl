module GRAPH500_Kernel_1 {
        use GRAPH500_Kronecker_Generator;
	use Time;
	use BlockDist;
	// ijw is an array of 2 dimensions
	// representing the edge list of the graph along with weights
	proc kernel_1(ijw:[]){
		write("Starting kernel_1 ... ");
		// Improve parallelizing here
                var ijwCols = ijw.domain.dim(2).high;
                var innerD: subdomain(ijw.domain) = [1..2, 1..ijwCols];
 
                write("  Starting max reduce on ijw with innerD ... ");
		var maxVal = max reduce ijw[innerD];
		writeln("done.");

		/*for k in ijw.domain.dim(2) {
			var i = ijw[1,k];
			var j = ijw[2,k];
			if (i != j) then {
				if (i >= max) then max = i;
				if (j >= max) then max = j;
			}
		}*/
		
		maxVal += 1; // going away from zero labels

		var D: domain(2) dmapped Block(boundingBox=[1..maxVal,1..maxVal]) = [1 .. maxVal, 1..maxVal];
		var G: [D] int = -1;
		
		write("  Creating D matrix ... ");
		forall k in ijw.domain.dim(2){
			var i = ijw[1,k] + 1; // going away from zero labels
			var j = ijw[2,k] + 1; // going away from zero labels
			var w = ijw[3,k];
			if (i != j) then {
				G[i,j] = w;
				G[j,i] = w;				
			} 
		}
		writeln("done.");
		writeln("  Returning G ... done.");
		writeln("done.");
		return G;
	}
	
	config var scale = 15;
	// Uncomment for testing kernel_1 alone
	proc main(){
	    writeln("Starting main()");
	    var t: Timer;
	    
	    t.start();
	    var G = generate_edge_list(scale, 16, 255);
	    t.stop();
	    writeln("Graph generation took: ", t.elapsed());
	    t.clear();
	    t.start();
	    kernel_1(G);
	    t.stop();
	    writeln("Kernel_1 took: ", t.elapsed());
	    writeln("Ending main()");
        }
}
