// Copyright (c) 2014-2015 Alan Meeson
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

function  spectralDrawing(graph, dims) {
	//Layout a graph using the spectral algorithm in 
	//Koren, Y, (2005) Drawing graphs by eigenvectors: Theory and practice
	//
	// This function performs an eigen decomposition on the graph matrix provided.
	// Power iteration is used to calculate the eigenvectors sequentially so that
	// only those required to cover the specified number of dimensions (dims) are 
	// calculated.
	
	//constants
	var MAX_ITTERATIONS = 100; //We use power iteration, this is analogous to wall time to avoid infinite loops.
	var epsilon = Math.pow(10, -8); //stopping criterion. Smaller values will use more iterations (if allowed).
	
	var num_elements = graph.length; //number of nodes in graph
	var D = deg(graph); //degree of each node in graph (number of connections).
	
	dims = dims + 1; //add one to the dims to allow for the first eigen vector
	var u = new Array(dims);//declare the eigen vector matrix
	u[0] = normalize(ones(num_elements)); //create & normalize the first eigen vector
	for (var i = 1; i < dims; i++) u[i] = zeros(num_elements); //create empty space for the other eigen vectors
	
	//Power iteration to determine the remaining eigen vectors.
	for (var k=1; k < dims; k++) { //for each eigen vector after the first, 
		//initialize eigen vector with random values
		var uhk = normalize(rand(num_elements));
		
		var itt_count = 0; //we are allowing a max of 100 iterations, to avoid hanging and infinite loops. (specified above in constants)
		var stop = false; //stopping criterion flag.
		while (!stop) { // do...while I used flags to keep it consistent with my matlab implementation
			
			//D-orthogonalize against previous eigenvectors
			var uk = uhk.slice();
			for (var l = 0; l < k; l++) {						
				var ul = u[l]; //extract the l-th eigen vector
				
				//Calculate (uk'.D.ul)/(ul'.D.ul)
				var top = 0;
				var bottom = 0;
				for (var vmi = 0; vmi < uk.length; vmi++) {
					top += (uk[vmi] * D[vmi] * ul[vmi]);
					bottom += (ul[vmi] * D[vmi] * ul[vmi]);
				}
				var ratio = top / bottom;
				
				//uk = uk - ((uk' . D . ul) / (ul' . D ul)) . ul
				for (var vsi = 0; vsi < uk.length; vsi++) {
					uk[vsi] = uk[vsi] - (ratio * ul[vsi]);
				}
			}
				
			//multiply with .5(I+D^-1 A)
			for (var i = 0; i < uhk.length; i++) {
				uhk[i] = 0.5 * (uk[i] + (dot(graph[i], uk) / D[i]));
			}
			uhk = normalize(uhk);
			
			itt_count = itt_count + 1;
			stop = (itt_count > 100) | !(dot(uhk, uk) < (1-epsilon));
		}
		u[k] = uhk.slice();	
	}

	//discard the first eigenvector which should be [ones].
	var v = new Array(u.length);
	for (var i=0; i < u.length; i++) {
		v[i] = new Array(u[i].length);
		for (var j=0; j < u[i].length; j++) v[i][j] = u[i][j];
	}
	return v;
}

function deg(graph) {
	//Calculate the degree of each node from the graph matrix.
	var d = zeros(graph.length);
	
	//degree of node i is the sum of the weights of all edges connected to it.
	for (var i = 0; i < graph.length; i++) {
		var node_degree = 0;
		for (var j = 0; j < graph[i].length; j++) {
			node_degree += graph[i][j];
		}
		d[i] = node_degree;
	}
	
	return d;
}

function dot(a,b) {
	//inner product of two vectors
	var d = 0;
	for (var i = 0; i < a.length; i++) {
		d += a[i] * b[i];
	}
	return d;
}

function euclideanDistance(coordinates) {
	//calculate the euclidean distance between two points/vectors.
	// used for normalization.
	var d = 0;
	
	for (var i = 0; i < coordinates.length; i++) {
		d += Math.pow(coordinates[i], 2);
	}
	return Math.sqrt(d);
}

function normalize(arr) {
	//normalizes a vector = arr/||arr||
	var d = euclideanDistance(arr);
	var narr = new Array(arr.length);
	for (var i = 0; i < arr.length; i++) {
		narr[i] = arr[i] / d;
	}
	
	return narr;
}

function rand(n) {
	//create a vector of length n and fill with random numbers.
	var arr = new Array(n);
	for (var i = 0; i < n; i++) arr[i] = Math.random();
	return arr;
}

function zeros(n) {
	//create a vector filled with zeros
	var arr = new Array(n);
	for (var i = 0; i < n; i++) arr[i] = 0;
	return arr;
}

function ones(n) {
	//create a vector filled with ones
	var arr = new Array(n);
	for (var i = 0; i < n; i++) arr[i] = 1;
	return arr;
}