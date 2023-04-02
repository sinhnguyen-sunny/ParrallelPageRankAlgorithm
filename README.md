# ParrallelPageRankAlgorithm

For share-memory parrallel pageRank algorithm
1.Install the required software: You will need to install the Metis library and OpenMP if they are not already installed on both laptops.

2.Download the Google web graph dataset: You can download the dataset from the Google Web Graph page. Make sure you download the appropriate version for your needs.

4.Partition the graph using Metis: Use the Metis library to partition the graph into smaller subgraphs. You can use the code provided above as a starting point.

5.Transfer data: Transfer the partitioned subgraphs and any other necessary data between the two laptops.

6.Run the code: Run the optimized code on each laptop using the partitioned subgraph data.

7.Combine the results: Combine the results from each laptop to get the final PageRank values for the entire graph.
