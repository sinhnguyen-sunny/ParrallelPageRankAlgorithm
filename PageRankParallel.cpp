#include <ctime>
#include <cstdio>
#include <direct.h>
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <sstream>
#include <string>

using namespace std;


int main() {

    /*************************** TIME, VARIABLES ***************************/

    ofstream resultFile;
    resultFile.open("./result.txt");

    // Keep track of the execution time
    clock_t begin, end;
    double time_spent;
    begin = clock();
    int numthreads = 16;
    int granularity = 8;

    /******************* OPEN FILE + NUM OF NODES/EDGES ********************/


    string fileName = "./web-BerkStan.txt";
    ifstream file(fileName);

    if (!file) { //Check open correctly
        resultFile << "Error cannot open file!" << endl;
        return EXIT_FAILURE;
    }

    int n = 0, e = 0;
    string line;

    //while (ch == '#') {
    //    fgets(str, 100 - 1, fp);
    //    //print title of the data set
    //    std::cout << str << endl;
    //    sscanf(str, "%*s %d %*s %d", &n, &e); //number of nodes
    //    ch = getc(fp);
    //}

    streampos oldpos;
    int counter = 0;
    while (getline(file, line) && line[0] == '#') {

        if (counter == 2) {
            istringstream ss(line);
            string str1, str2, str3;
            int a, b;
            ss >> str1 >> str2 >> a >> str3 >> b;
            n = a;
            e = b;
        }
        else if (counter == 0) {
            resultFile << line.substr(2);
        }
        oldpos = file.tellg();
        counter++;
    }

    file.seekg(oldpos);

    // DEBUG: Print the number of nodes and edges, skip everything else
    resultFile << "\n" << "Graph data:" << "\n\n" << " Nodes: " << n << ", Edges: " << e << "\n" << endl;



    float* val = new float[e];
    int* col_ind = new int[e];
    int* row_ptr = new int[n + 1];

    // The first row starts at position 0
    row_ptr[0] = 0;

    int fromnode, tonode;
    int cur_row = 0;
    int i = 0;
    int j = 0;
    // Elements for row
    int elrow = 0;
    // Cumulative numbers of elements
    int curel = 0;
    string str1;
    string contentLine;

    while (!file.eof()) {

        getline(file, contentLine);
        if (!contentLine.empty()) {
            fromnode = stoi(contentLine.substr(0, contentLine.find("	")));
            tonode = stoi(contentLine.substr(contentLine.find("	") + 1));

            // DEBUG: print fromnode and tonode
            //printf("From: %d To: %d\n",fromnode, tonode);

            if (fromnode > cur_row) { // change the row
                curel = curel + elrow;
                for (int k = cur_row + 1; k <= fromnode; k++) {
                    row_ptr[k] = curel;
                }
                elrow = 0;
                cur_row = fromnode;
            }
            val[i] = 1.0;
            col_ind[i] = tonode;
            elrow++;
            i++;
        }

    }
    row_ptr[cur_row + 1] = curel + elrow - 1;



    // Fix the stochastization
    int* out_link = new int[n];
    for (i = 0; i < n; i++) {
        out_link[i] = 0;
    }

    /* DEBUG: row pointer test
    printf("\nRow_ptr:\n");
    for (i=0; i<n; i++){
      printf("%d ", row_ptr[i]);
    }
    printf("\n");
    */

    int rowel = 0;
    for (i = 0; i < n; i++) {
        if (row_ptr[i + 1] != 0) {
            rowel = row_ptr[i + 1] - row_ptr[i];
            out_link[i] = rowel;
        }
    }


    int curcol = 0;
    for (i = 0; i < n; i++) {
        rowel = row_ptr[i + 1] - row_ptr[i];
        for (j = 0; j < rowel; j++) {
            val[curcol] = val[curcol] / out_link[i];
            curcol++;
        }
    }

    /* DEBUG: val print test
    for(i=0; i<e; i++){
        printf("%f ", val[i]);
    }*/

    /******************* INITIALIZATION OF P, DAMPING FACTOR ************************/

    // Set the damping factor 'd'
    float d = 0.85f;

    // Initialize p[] vector
    float* p = new float[n];
    for (i = 0; i < n; i++) {
        p[i] = 1.0f / n;
    }

    // Set the looping condition and the number of iterations 'k'
    int looping = 1;
    int k = 0;

    // Set 'parallel' depending on the number of threads
    int parallel = 0;
    if (numthreads >= 2) {
        parallel = 1;
    }

    // Initialize new p vector
    float* p_new = new float[n];

    /*************************** PageRank LOOP  **************************/

    while (looping) {

        // Initialize p_new as a vector of n 0.0 cells
        for (i = 0; i < n; i++) {
            p_new[i] = 0.0;
        }

        int rowel = 0;
        int curcol = 0;

        // Page rank modified algorithm + parallelization
#pragma omp parallel for schedule(static) if(parallel) num_threads(numthreads)
        for (i = 0; i < n; i = i + granularity) {
            rowel = row_ptr[i + 1] - row_ptr[i];
            for (j = 0; j < rowel; j++) {
                p_new[col_ind[curcol]] = p_new[col_ind[curcol]] + val[curcol] * p[i];
                curcol++;
            }
        }

        /*DEBUG: print pnew
        for (i=0; i<n; i++){
          printf("%f ", p_new[i]);
        }*/

        // Adjustment to manage dangling elements 
        for (i = 0; i < n; i++) {
            p_new[i] = d * p_new[i] + ((1.0f - d) / n);
        }

        /*DEBUG: print pnew after the damping factor multiplication
        for (i=0; i<n; i++){
          printf("%f ", p_new[i]);
        }*/

        // TERMINATION: check if we have to stop
        float error = 0.0f;
        for (i = 0; i < n; i++) {
            error = error + abs(p_new[i] - p[i]);
        }
        //if two consecutive instances of pagerank vector are almost identical, stop
        if (error < 0.000001) {
            looping = 0;
        }

        // Update p[]
        for (i = 0; i < n; i++) {
            p[i] = p_new[i];
        }

        // Increase the number of iterations
        k = k + 1;
    }

    /*************************** CONCLUSIONS *******************************/

    // Stop the timer and compute the time spent
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    // Sleep a bit so stdout is not messed up
    //Sleep(500);

    // Print results
    resultFile << "\n" << "Number of iteration to converge: " << k << "\n" << endl;

    resultFile << "Final Pagerank values:" << "\n\n";

    for (i = 0; i < n; i++) {
        resultFile << "Node " << i << ": " << fixed << setprecision(12) << p[i] << endl;
    }
    resultFile << "\n\n" << "Time spent: " << time_spent << " seconds." << endl;

    resultFile.close();
    return EXIT_SUCCESS;
}