use BlockDist;
config var dim = 2,
	   clusters = 4,
	   points = 2000,
	   iterations = 50;

proc main() {
  writeln("numLocales: ", numLocales);
  for loc in Locales do
    on loc do
      writeln("hello from locale: ", here.id);
  writeln();
  forall loc in Locales do
    writeln("now on: ", here.id);
  writeln();
  
  const Space = {1..8,1..8};
  const BlockSpace = Space dmapped Block(boundingBox=Space);
  var BA: [BlockSpace] int;
  forall (i,j) in Space do{
    writeln("where is BA[i,j]: ", BA[i,j].locale.id);
    writeln("i: ", i, " j: ", j);
  }

// Very nice.
//
// Chapel variables are stored using the memory of the locale
// executing the task that encounters the variable declaration.
// Thus, in the following code, x is declared on locale 0 and
// y is declared on locale 1 if numLocales > 1 (0 otherwise).
//

{
  var x: int = 2;
  on Locales[1 % numLocales] {
    var y: int = 3;
    writeln("From locale ", here.id, ", x is: ", x, " and y is: ", y);
    on Locales[0] {
      writeln("From locale 0, x is: ", x, " and y is: ", y);
    }
  }
  writeln();
}
 
  writeln ("dimension: ", dim);
  writeln ("clusters ", clusters);
  writeln ("points: ", points);
  writeln ("iterations: ", iterations);
  
  

}
