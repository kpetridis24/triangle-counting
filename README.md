# Vertex-wise-triangle-counting
Counting number of triangles that each node contributes to in large sparse matrices.
For this algorithm we use large sparse matrices with in CSC or CSR format only!

~ Function V3 ~

The logic followed in function V3 is as follows. First we access the
nodes-rows (row_ptr table) and their corresponding column elements (col_index).
For each possible col_index data pair of a given node, a search is started in the
line, the pointer of which is defined by the first element of the pair. Searhing the data
col_index contained in the respective row-node, if we locate its second component
pair, we have proven the existence of a triangle, and the appropriate pointers increase. 
For example, lets say node zero is adjacent to the elements {1, 2, 3, 4}. For
pair (1,2), a pointer will be referenced in line 1 and will
access its col_index data. In case the second
element of the pair is also located on the new line (line 1), there is a connection
of nodes (0-1), (0-2), (1-2) therefore they participate in a triangle.

~ Function V4 ~

The algorithm that implements V4 follows a similar logic to V3 with respect to
accessing nodes, indirectly using array multiplication for
finding the node contribution to existing triangles. This time we receive every
pair of neighboring nodes, and look for the number of their common col_index elements.
Considering the above example, node zero is adjacent to node 1, so
our pointer will access the col_index table data corresponding to the row
1, counting the common neighbors of the nodes (0,1). This result expresses
substantially the number of non-zero existing multiplications (Hadamard
product), therefore the final value of the index will be indicated in c3 [0,1]. 
The algorithm uses merge function to compare lists of neighboring nodes
in order to significantly reduce complexity. The process continues for
all neighbors of node zero, ending in the calculation of c3 [0, {1,2,3,4}]
by summing the individual c3 [0, x] and repeating the same for all of them
remaining nodes.


~ Code Weaknesses ~

none detected


