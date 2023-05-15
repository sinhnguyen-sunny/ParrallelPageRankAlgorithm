# ParrallelPageRankAlgorithm

## Alogrithms:

    * [pagerank_group_slide_link](https://docs.google.com/presentation/d/1NwhEI3SX1bi2esqgvTDS6Oaywb2_tXpbr0lyePAbMRY/edit#slide=id.g21039aed30c_0_7)

## For serial small graph

1. Install g++ by command: "sudo apt-get install g++"
2. Download the Eigen 3.3.9 source code from the official website: https://eigen.tuxfamily.org/index.php?title=Main_Page#Download
3. Build and install Eigen
4. cd to the folder "serial_small_graph"
5. Run "g++ -I eigen-3.3.9 pagerank.cpp -o pagerank"

## For serial big graph

1. cd to the folder "serial_big_graph"
2. Run "gcc -o pagerank serial_page_rank.c -lm" then run "./pagerank > ./output.txt"

## For paralell

1. cd to the folder "parallel"
2. Run "g++ -fopenmp page_rank_parallel.cpp -o pagerank" then run "./pagerank > ./output.txt"
