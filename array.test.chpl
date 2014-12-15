use Random;
proc main(){
  var D = {0..2,0..1};
  var A: [D] real;
  var rand: RandomStream = new RandomStream();
  rand.fillRandom(A);
  writeln(A);
  var B: [D] real = ((1.1, 1.2),(2.1,2.2),(3.1,3.2));
  writeln(B);
}
