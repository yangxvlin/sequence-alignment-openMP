// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/ and
// fixed an error when initializing the dp array :-)
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>

using namespace std;

int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans);

// Return current time, for performance measurement
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t) 1000000 + tv.tv_usec;
}


// Driver code
int main() {
    int misMatchPenalty;
    int gapPenalty;
    std::string gene1;
    std::string gene2;
    std::cin >> misMatchPenalty;
    std::cin >> gapPenalty;
    std::cin >> gene1;
    std::cin >> gene2;
    std::cout << "misMatchPenalty=" << misMatchPenalty << std::endl;
    std::cout << "gapPenalty=" << gapPenalty << std::endl;

    int m = gene1.length(); // length of gene1
    int n = gene2.length(); // length of gene2
    int l = m + n;
    int xans[l + 1], yans[l + 1];

    uint64_t start = GetTimeStamp();

    // calling the function to calculate the result
    int penalty = getMinimumPenalty(gene1, gene2,
                                    misMatchPenalty, gapPenalty,
                                    xans, yans);

    // print the time taken to do the computation
    printf("Time: %ld us\n", (uint64_t)(GetTimeStamp() - start));

    // postprocessing of the answer, for printing results

    // Since we have assumed the answer to be n+m long,
    // we need to remove the extra gaps in the starting
    // id represents the index from which the arrays
    // xans, yans are useful
    int id = 1;
    int i;
    for (i = l; i >= 1; i--) {
        if ((char) yans[i] == '_' && (char) xans[i] == '_') {
            id = i + 1;
            break;
        }
    }

    // Printing the final answer
    std::cout << "Minimum Penalty in aligning the genes = ";
    std::cout << penalty << std::endl;
	std::cout << "The aligned genes are :" << std::endl;
	for (i = id; i <= l; i++)
	{
		std::cout<<(char)xans[i];
	}
	std::cout << "\n";
	for (i = id; i <= l; i++)
	{
		std::cout << (char)yans[i];
	}
	std::cout << "\n";

    return 0;
}

/******************************************************************************/
/* Do not change any lines above here.            */
/* All of your changes should be below this line. */
/******************************************************************************/
#include <math.h>       /* ceil */
#include <omp.h>
#include <iomanip>      // std::setw

// uncomment to enable debug mode
// #define DEBUG 0

inline int min3(int a, int b, int c) {
    if (a <= b && a <= c) {
        return a;
    } else if (b <= a && b <= c) {
        return b;
    } else {
        return c;
    }
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
inline int **new2d(int width, int height) {
    int **dp = new int *[width];
    size_t size = width;
    size *= height;
    int *dp0 = new int[size];
    if (!dp || !dp0) {
        std::cerr << "getMinimumPenalty: new failed" << std::endl;
        exit(1);
    }
    dp[0] = dp0;
    for (int i = 1; i < width; i++)
        dp[i] = dp[i - 1] + height;

    return dp;
}

// function to find out the minimum penalty
// return the maximum penalty and put the aligned sequences in xans and yans
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap,
                      int *xans, int *yans) {
    int i, j; // intialising variables

    int m = x.length(); // length of gene1
    int n = y.length(); // length of gene2
    int row = m + 1, col = n + 1;

    // table for storing optimal substructure answers
    int **dp = new2d(row, col);
//	size_t size = m + 1;
//	size *= n + 1;
//	memset (dp[0], 0, size);

    // intialising the table
    for (i = 0; i <= m; i++) {
        dp[i][0] = i * pgap;
    }
    for (i = 0; i <= n; i++) {
        dp[0][i] = i * pgap;
    }

    // calcuting the minimum penalty
//	for (i = 1; i <= m; i++)
//	{
//		for (j = 1; j <= n; j++)
//		{
//			if (x[i - 1] == y[j - 1])
//			{
//				dp[i][j] = dp[i - 1][j - 1];
//			}
//			else
//			{
//				dp[i][j] = min3(dp[i - 1][j - 1] + pxy ,
//						dp[i - 1][j] + pgap ,
//						dp[i][j - 1] + pgap);
//			}
//		}
//	}

    #ifdef DEBUG
		cout.fill(' ');
        for (i = 0; i < row; i++) {
            for (j = 0; j < col; j++) {
                // Prints ' ' if j != n-1 else prints '\n'           
                cout << setw(3) << dp[i][j] << " "; 
			}
			cout << "\n";
		}
        cout << ">>>> \n";
    #endif

    // Tile parallel
    int n_threads = 18;
    omp_set_num_threads(n_threads);

    int tile_width = (int) ceil((1.0*m) / n_threads), tile_length = (int) ceil((1.0*n) / n_threads);
    int num_tile_in_width = (int) ceil((1.0*m) / tile_width);
    int num_tile_in_length = (int) ceil((1.0*n) / tile_length);;

    // cout << "tile_width :" << tile_width << " "<< "tile_length :" << tile_length<< endl;
    // cout << "num_tile_in_width :" << num_tile_in_width << " num_tile_in_length :" << num_tile_in_length<<endl;

    // There will be tile_width+COL-1 lines in the output
    for (int line = 1; line <= (num_tile_in_width + num_tile_in_length - 1); line++) {
        /* Get column index of the first element in this line of output.
           The index is 0 for first tile_width lines and line - tile_width for remaining
           lines  */
        int start_col = max(0, line - num_tile_in_width);

        /* Get count of elements in this line. The count of elements is
           equal to minimum of line number, tile_length-start_col and num_tile_in_width */
        int count = min(line, min((num_tile_in_length - start_col), num_tile_in_width));

        /* Print elements of this line */
        #pragma omp parallel for
        for (int z = 0; z < count; z++) {
            // cout << (min(num_tile_in_width, line)-z-1)  << " " << (start_col+z)  << "->" << (min(num_tile_in_width, line)-z-1)*tile_width +1<< " " << (start_col+z)*tile_length +1<< endl;

            int tile_i_start = (min(num_tile_in_width, line)-z-1)*tile_width +1,
                tile_j_start = (start_col+z)*tile_length +1;

            for (int i = tile_i_start; i < min(tile_i_start + tile_width, row); i++) {
                for (int j = tile_j_start; j < min(tile_j_start + tile_length, col); j++) {
                    // cout << "(i, j) ("<< i << ", " << j << ")" << endl;

                    if (x[i - 1] == y[j - 1]) {
                        dp[i][j] = dp[i - 1][j - 1];
                    } else {
                        dp[i][j] = min3(dp[i - 1][j - 1] + pxy ,
                                dp[i - 1][j] + pgap ,
                                dp[i][j - 1] + pgap);
                    }
                }
            }
            // cout << "tile end" << endl;


            // cout << (min(num_tile_in_width, line)-j-1)  << " " << (start_col+j)  << "->" << (min(num_tile_in_width, line)-j-1)*tile_width +1<< " " << (start_col+j)*tile_length +1<< endl;
        }
//            printf("%5d ", matrix[minu(tile_width, line)-j-1][start_col+j]);

        /* Ptint elements of next diagonal on next line */
//        printf("\n");
    }

    #ifdef DEBUG
		cout.fill(' ');
        for (i = 0; i < row; i++) {
            for (j = 0; j < col; j++) {
                // Prints ' ' if j != n-1 else prints '\n'           
                cout << setw(3) << dp[i][j] << " "; 
			}
			cout << "\n";
		}
        cout << ">>>> \n";
    #endif

//    for (i = 1; i <= (tile_width + tile_length - 1); i++) {
//        // std::cout << "* i " << i << "\n";
//        int z1 = i < tile_length ? 0 : i - tile_length + 1;
//        int z2 = i < tile_width ? 0 : i - tile_width + 1;
//        // std::cout << "* num_threads " << i - z2 - z1 + 1 << "\n";
////#pragma omp parallel for private(j)
//        for (j = z1 + 1; j <= i - z2 + 1; j++) {
//            int dpy = i - j + 2;
//            cout << i * tile_width << " " << dpy * tile_length << endl;
//        }
//    }


//    for (int ii = 1; ii < row; ii += tile_width) {
//        // pragma for
//        for (int jj = 1; jj < col; jj += tile_length) {
//            // Go over each block sequentially
//            for (int i = ii; i < min(ii + tile_width, row); i++) {
//                for (int j = jj; j < min(jj + tile_length, col); j++) {
//
//                }
//            }
//        }
//    }
    // exit(1);
    // Reconstructing the solution
    int l = n + m; // maximum possible length

    i = m;
    j = n;

    int xpos = l;
    int ypos = l;

    while (!(i == 0 || j == 0)) {
        // cout << "(i, j) ("<< i << ", " << j << ")" << endl;

        if (x[i - 1] == y[j - 1]) {
            xans[xpos--] = (int) x[i - 1];
            yans[ypos--] = (int) y[j - 1];
            i--;
            j--;
        } else if (dp[i - 1][j - 1] + pxy == dp[i][j]) {
            xans[xpos--] = (int) x[i - 1];
            yans[ypos--] = (int) y[j - 1];
            i--;
            j--;
        } else if (dp[i - 1][j] + pgap == dp[i][j]) {
            xans[xpos--] = (int) x[i - 1];
            yans[ypos--] = (int) '_';
            i--;
        } else if (dp[i][j - 1] + pgap == dp[i][j]) {
            xans[xpos--] = (int) '_';
            yans[ypos--] = (int) y[j - 1];
            j--;
        }
    }
    while (xpos > 0) {
        if (i > 0) xans[xpos--] = (int) x[--i];
        else xans[xpos--] = (int) '_';
    }
    while (ypos > 0) {
        if (j > 0) yans[ypos--] = (int) y[--j];
        else yans[ypos--] = (int) '_';
    }

    int ret = dp[m][n];

    delete[] dp[0];
    delete[] dp;

    return ret;
}

// g++ -fopenmp -o xuliny-seqalignomp xuliny-seqalignomp.cpp -O3