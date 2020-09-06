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
	int* xans, int* yans) {
	int i, j; // intialising variables
	
	int m = x.length(); // length of gene1
	int n = y.length(); // length of gene2

	// swap two string if m < n
	if (m < n) {
		string temp { x };
		x = y;
		y = temp;

		m = x.length(); // length of gene1
		n = y.length(); // length of gene2

		int * temp_ans = xans;
		xans = yans;
		yans = temp_ans;
	}
	
	// table for storing optimal substructure answers
    int row = m+1, col = n+1;  // size of original & swapped dp matrix
	int r = row + n; // #rows of diamond matrix
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
	// unntil last n rows
	int until_last_n_row = r - n;
	for (int i = 0; i < until_last_n_row; i++) {
		// int j_skew = i-r+1;
        int upper  = min(i+1, col);
		// int i_swaped = i - row;
		// cout << i << " -> " << i_swaped << " j_skew: "<< j_skew << endl;
		// #pragma omp parallel for
		for (int j = 0; j < upper; j++) {
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

	// 1st row to be swapped
	i = 0;
	// int upper  = min(until_last_n_row+i, col);
	for (j = 1; j < col; j++) {
		int left_x     = until_last_n_row - 1, 
			left_y     = j - 1, 
			up_x       = until_last_n_row - 1,  
			up_y       = j, 
			left_up_x  = until_last_n_row - 2, 
			left_up_y  = j - 1;

			if (x[i - j - 1] == y[left_up_y]) {
				dp[i][j] = dp[left_up_x][left_up_y];
			}
			else {
				dp[i][j] = min3(dp[left_up_x][left_up_y] + pxy  ,
								dp[left_up_x][left_up_y]     + pgap ,
								dp[up_x][up_y]     + pgap);
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

	// 2nd row to be swapped
	i = 1;
	// upper  = min(until_last_n_row+i, col);
	for (j = 2; j < col; j++) {
		int left_x     = 0, 
			left_y     = j - 1, 
			up_x       = 0,  
			up_y       = j, 
			left_up_x  = until_last_n_row - 1, 
			left_up_y  = j - 1;

		if (x[i - j - 1] == y[left_up_y]) {
			dp[i][j] = dp[left_up_x][left_up_y];
		}
		else {
			dp[i][j] = min3(dp[left_up_x][left_up_y] + pxy  ,
							dp[left_up_x][left_up_y]     + pgap ,
							dp[up_x][up_y]     + pgap);
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

	// 3rd ~ n row to be swapped
	for (i = 2; i < col; i++) {
		// int upper  = min(until_last_n_row+i, col);
		for (int j = i+1; j < col; j++) {
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

	cout << "answer: " << dp[(r - 1) % row][col-1] << endl;

	// exit(1);

	while ( !(i == 0 || j == 0)) {


		dp_i = i + j;
		dp_j = j;

		// cout << "(" << dp_i << "," << dp_j << ") -> " << "(i, j) ("<< i << ", " << j << ")" << endl;

		if (x[i - 1] == y[j - 1]) {
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[(dp_i - 2) % row][dp_j - 1] + pxy == dp[dp_i % row][dp_j]) {
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[(dp_i - 1) % row][dp_j] + pgap == dp[dp_i % row][dp_j]) {
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)'_';
			i--;
		}
		else if (dp[(dp_i - 1) % row][dp_j - 1] + pgap == dp[dp_i % row][dp_j]) {
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

	int ret = dp[(r - 1) % row][col-1];

	delete[] dp[0];
	delete[] dp;
	
	return ret;
}
