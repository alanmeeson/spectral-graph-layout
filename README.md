# spectral-graph-layout
A graph layout algorithm which uses eigen decomposition and stress majorization to generate edge-node graph layouts on the fly in a manner which could be embedded in HTML/SVG documents.

# Known Issues
Attempting to layout graphs in which two or more nodes share the same connectedness often results in the nodes in question being placed in the same location. An example of this issue occurs in the example used in the demo page.

# Future Development Possibilities
1. Investigate fudge factors as applied to the connectedness graph for separating otherwise identical nodes.
2. Investigate post processing algorithms for separating overlaped nodes.
3. Investigate sub-dividing the graph and laying out (or re-laying out) the sub-graphs

# Change Log
V0.2 Early commits to GitHub
+ Added Demo.html providing annotated example.
+ Added functions for generating random symmetric matrices, and matrix addition.

V0.1 Initial prototype as converted from the basic matlab code.
+ Implemented eigen decomposition and stress majorization
+ Basic example in HTML document