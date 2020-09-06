// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/ and
// fixed an error when initializing the dp array :-)
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>

using namespace std;

int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int* xans, int* yans);

// Return current time, for performance measurement
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}


// Driver code
int main(){
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
	int l = m+n;
	int xans[l+1], yans[l+1];

	uint64_t start = GetTimeStamp ();

	// calling the function to calculate the result
	int penalty = getMinimumPenalty(gene1, gene2,
		misMatchPenalty, gapPenalty,
		xans,yans);
	
	// print the time taken to do the computation
	printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
	
	// postprocessing of the answer, for printing results

	// Since we have assumed the answer to be n+m long,
	// we need to remove the extra gaps in the starting
	// id represents the index from which the arrays
	// xans, yans are useful
	int id = 1;
	int i;
	for (i = l; i >= 1; i--)
	{
		if ((char)yans[i] == '_' && (char)xans[i] == '_')
		{
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
#include <omp.h>
#include <iomanip>      // std::setw

// uncomment to enable debug mode
// #define DEBUG 0

int min3(int a, int b, int c) {
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
int **new2d (int width, int height)
{
	int **dp = new int *[width];
	size_t size = width;
	size *= height;
	int *dp0 = new int [size];
	if (!dp || !dp0)
	{
	    std::cerr << "getMinimumPenalty: new failed" << std::endl;
	    exit(1);
	}
	dp[0] = dp0;
	for (int i = 1; i < width; i++)
	    dp[i] = dp[i-1] + height;

	return dp;
}

// function to find out the minimum penalty
// return the maximum penalty and put the aligned sequences in xans and yans
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap,
	int* xans, int* yans)
{
	int i, j; // intialising variables
	
	int m = x.length(); // length of gene1
	int n = y.length(); // length of gene2
	
    if (m > 70000 && n > 70000) {
        int **dp = new2d (m+1, n+1);

        #pragma omp parallel for
        for (i = 0; i <= m; i++)
        {
            dp[i][0] = i * pgap;
        }
        #pragma omp parallel for
        for (i = 1; i <= n; i++)
        {
            dp[0][i] = i * pgap;
        }

        // calcuting the minimum penalty
        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= n; j++)
            {
                if (x[i - 1] == y[j - 1])
                {
                    dp[i][j] = dp[i - 1][j - 1];
                }
                else
                {
                    dp[i][j] = min3(dp[i - 1][j - 1] + pxy ,
                            dp[i - 1][j] + pgap ,
                            dp[i][j - 1] + pgap);
                }
            }
        }

        // Reconstructing the solution
        int l = n + m; // maximum possible length
        
        i = m; j = n;
        
        int xpos = l;
        int ypos = l;
        
        while ( !(i == 0 || j == 0))
        {
            if (x[i - 1] == y[j - 1])
            {
                xans[xpos--] = (int)x[i - 1];
                yans[ypos--] = (int)y[j - 1];
                i--; j--;
            }
            else if (dp[i - 1][j - 1] + pxy == dp[i][j])
            {
                xans[xpos--] = (int)x[i - 1];
                yans[ypos--] = (int)y[j - 1];
                i--; j--;
            }
            else if (dp[i - 1][j] + pgap == dp[i][j])
            {
                xans[xpos--] = (int)x[i - 1];
                yans[ypos--] = (int)'_';
                i--;
            }
            else if (dp[i][j - 1] + pgap == dp[i][j])
            {
                xans[xpos--] = (int)'_';
                yans[ypos--] = (int)y[j - 1];
                j--;
            }
        }
        while (xpos > 0)
        {
            if (i > 0) xans[xpos--] = (int)x[--i];
            else xans[xpos--] = (int)'_';
        }
        while (ypos > 0)
        {
            if (j > 0) yans[ypos--] = (int)y[--j];
            else yans[ypos--] = (int)'_';
        }

        int ret = dp[m][n];

        delete[] dp[0];
        delete[] dp;
        
        return ret;
    }

	// table for storing optimal substructure answers
    int r = m+1, col = n+1;
	int row = r + n;
	int **dp = new2d (row, col);

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

    int n_threads = 8;
    omp_set_num_threads(n_threads);

	// calcuting the minimum penalty
	for (int i = 0; i < row; i++) {
        int upper = min(i+1, col),
            lower = max(0,   i-row+1);

		#pragma omp parallel for
		for (int j = lower; j < upper; j++) {
            if (j == 0 && i <= m) {
                dp[i][j] = i * pgap;
            } 
            else if (i == j) {
                dp[i][j] = i * pgap;
            } else {
                int left_up_x  = i - 2, 
                	left_up_y  = j - 1,
                    left = i - 1;

                if (x[i - j - 1] == y[left_up_y]) {
                	dp[i][j] = dp[left_up_x][left_up_y];
                	// cout << "equal" << endl;
                }
                else {
                	// cout << "min" << endl;
                	// cout << dp[left_up_x][left_up_y] << " " << dp[left_x][left_y] << " " << dp[up_x][up_y] << endl;
                	dp[i][j] = min3(dp[left_up_x][left_up_y] + pxy  ,
                					dp[left][left_up_y]     + pgap ,
                					dp[left][j]     + pgap);
                }
                
                // cout << "(i, j) ("<< i << ", " << j << ") my up: (i-1, j-1) (" << i-1 << ", " << j-1 << ")" << endl;
                
                // int left_x     = i - 1, 
                // 	left_y     = j - 1, 
                // 	up_x       = i - 1,  
                // 	up_y       = j, 
                // 	// left_up_x  = i - 2, 
                // 	// left_up_y  = j - 1,
                // 	original_i = i - j,
                // 	original_j = j;
                
                // cout << "(" << original_i << "," << original_j << ") -> " << "(i, j) ("<< i << ", " << j << ") my up: (i-1, j-1) (" << up_x << ", " << up_y << ")" << endl;

                // if (x[original_i - 1] == y[original_j - 1]) {
                // 	dp[i][j] = dp[left_up_x][left_up_y];
                // 	// cout << "equal" << endl;
                // }
                // else {
                // 	// cout << "min" << endl;
                // 	// cout << dp[left_up_x][left_up_y] << " " << dp[left_x][left_y] << " " << dp[up_x][up_y] << endl;
                // 	dp[i][j] = min3(dp[left_up_x][left_up_y] + pxy  ,
                // 					dp[left_x][left_y]     + pgap ,
                // 					dp[up_x][up_y]     + pgap);
                // }			
                // cout << "(" << original_i << "," << original_j << ") -> " << "(i, j) ("<< i << ", " << j << ") my up: (i-1, j-1) (" << up_x << ", " << up_y << ")" << endl;

            }

		}
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

	// Reconstructing the solution
	int l = n + m; // maximum possible length
	
	i = m; j = n;
	
	int xpos = l;
	int ypos = l;
	
	int dp_i, dp_j;

	while ( !(i == 0 || j == 0))
	{
		dp_i = i + j;
		dp_j = j;

		if (x[i - 1] == y[j - 1])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[dp_i - 2][dp_j - 1] + pxy == dp[dp_i][dp_j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[dp_i - 1][dp_j] + pgap == dp[dp_i][dp_j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)'_';
			i--;
		}
		else if (dp[dp_i - 1][dp_j - 1] + pgap == dp[dp_i][dp_j])
		{
			xans[xpos--] = (int)'_';
			yans[ypos--] = (int)y[j - 1];
			j--;
		}
	}
	while (xpos > 0)
	{
		if (i > 0) xans[xpos--] = (int)x[--i];
		else xans[xpos--] = (int)'_';
	}
	while (ypos > 0)
	{
		if (j > 0) yans[ypos--] = (int)y[--j];
		else yans[ypos--] = (int)'_';
	}

	int ret = dp[row-1][col-1];

	delete[] dp[0];
	delete[] dp;
	
	return ret;
}

// g++ -fopenmp -o xuliny-seqalignomp xuliny-seqalignomp.cpp -O3