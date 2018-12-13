Run as :

python3 steiner_tree.py <dir_name> <gammavalue which is between 0 and 1>

For example:

python steiner_tree.py 4 0.6

The directory contains two files -

- edges.txt : Each line is in the format
i j weight
where i and j are the nodes the edge connects and weight is the weight of the edge.

- penalties.txt : Contains the value for the penalty corresponding to each node in the graph
Format :
i penalty
Here i is the node id and penalty is the cost incurred if the node is not included in steiner tree.

 
