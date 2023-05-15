#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

/* Code for:
	 sequential implementation of the PageRank algorithm with
	 Compressed Sparse Row (CSR) representation of matrix A
*/

int main() {

		/*************************** TIME, VARIABLES ***************************/

		// Keep track of the execution time
		clock_t begin, end;
		double time_spent;
		begin = clock();

		/******************* OPEN FILE + NUM OF NODES/EDGES ********************/

		// Open the data set
		// char filename[] = "/home/sinhnta/page_rank/dataset/web-NotreDame.txt";
		// char filename[] = "/home/sinhnta/page_rank/dataset/web-Stanford.txt";
		char filename[] = "/home/sinhnta/page_rank/dataset/web-BerkStan.txt";
		FILE* fp;
		if ((fp = fopen(filename, "r")) == NULL) {
				fprintf(stderr, "[Error] Cannot open the file");
				exit(1);
		}

		// Read the data set and get the number of nodes (n) and edges (e)
		int n, e;
		char ch;
		char str[100];
		ch = getc(fp);
		while (ch == '#') {
				fgets(str, 100 - 1, fp);
				//Debug: print title of the data set
				//printf("%s",str);
				sscanf(str, "%*s %d %*s %d", &n, &e); //number of nodes
				ch = getc(fp);
		}
		ungetc(ch, fp);

		// DEBUG: Print the number of nodes and edges, skip everything else
		printf("\nGraph data:\n\n  Nodes: %d, Edges: %d \n\n", n, e);

		/************************* CSR STRUCTURES *****************************/

		/* Compressed sparse row format:
			 - Val vector: contains 1.0 if an edge exists in a certain row
			 - Col_ind vector: contains the column index of the corresponding value in 'val'
			 - Row_ptr vector: points to the start of each row in 'col_ind'
		*/

		float* val = calloc(e, sizeof(float));
		int* col_ind = calloc(e, sizeof(int));
		int* row_ptr = calloc(n + 1, sizeof(int));

		// The first row always starts at position 0
		row_ptr[0] = 0;

		int fromnode, tonode;
		int cur_row = 0;
		int i = 0;
		int j = 0;
		// Elements for row
		int elrow = 0;
		// Cumulative numbers of elements
		int curel = 0;

		while (!feof(fp)) {

				fscanf(fp, "%d%d", &fromnode, &tonode);

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
		row_ptr[cur_row + 1] = curel + elrow;

		fclose(fp);

		// DEBUG: Print row_ptr, col_ind, and val vectors
		printf("\nCSR Structures:\n\n");
		printf("  Row_ptr: ");
		for (int i = 0; i < n+1; i++) {
				printf("%d ", row_ptr[i]);
		}
		printf("\n");
		printf("  Col_ind: ");
			for (int i = 0; i < e; i++) {
					printf("%d ", col_ind[i]);
			}
			printf("\n");

			printf("  Val: ");
			for (int i = 0; i < e; i++) {
					printf("%.1f ", val[i]);
			}
			printf("\n\n");

			/************************* PAGERANK ALGORITHM **************************/

			// PageRank parameters
			double alpha = 0.85;
			double tol = 1.0e-6;
			double pr_init = 1.0 / n;

			// PageRank vectors
			double* pr = calloc(n, sizeof(double));
			double* pr_next = calloc(n, sizeof(double));

			// Initialize the vectors with pr_init
			for (int i = 0; i < n; i++) {
					pr[i] = pr_init;
					pr_next[i] = pr_init;
			}

			// Calculate the dangling nodes vector
			double* dn = calloc(n, sizeof(double));
			for (int i = 0; i < n; i++) {
					dn[i] = 0;
			}
			for (int i = 0; i < e; i++) {
					dn[col_ind[i]] += val[i];
			}

			// Set the dangling nodes probability
			double dn_prob = 0;
			for (int i = 0; i < n; i++) {
					if (dn[i] == 0) {
							dn_prob += pr[i] / n;
					}
			}

			// Main loop of the algorithm
			int iter = 0;
			double err = tol + 1;
			while (err > tol) {
					// Perform matrix-vector multiplication
					for (int i = 0; i < n; i++) {
							double sum = 0;
							for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
									int col = col_ind[j];
									sum += val[j] * pr[col] / dn[col];
							}
							pr_next[i] = (1 - alpha) / n + alpha * (dn_prob + sum);
					}

					// Calculate the L1 norm of the error
					err = 0;
					for (int i = 0; i < n; i++) {
							err += fabs(pr[i] - pr_next[i]);
					}

					// Swap the vectors
					double* temp = pr;
					pr = pr_next;
					pr_next = temp;

					iter++;
			}

			end = clock();
			time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

			// DEBUG: Print the PageRank scores
			printf("PageRank Scores:\n\n");
			for (int i = 0; i < n; i++) {
					printf("  Node %d: %.6f\n", i, pr[i]);
			}
			printf("\n");

			// DEBUG: Print the number of iterations and the execution time
			printf("Number of iterations: %d\n", iter);
			printf("Execution time: %f seconds\n", time_spent);

			// Free the memory
			free(val);
			free(col_ind);
			free(row_ptr);
			free(pr);
			free(pr_next);
			free(dn);

			return 0;
}
