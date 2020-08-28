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

int min3(int a, int b, int c) {
    if (a < b && a < c) {
        return a;    
    } else if (b < a && b < c) {
        return b;
    } else {
        return c;    
    }
}

/******************************************************************************/
/* Do not change any lines above here.            */
/* All of your changes should be below this line. */
/******************************************************************************/

#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define MIN(a, b) ((a) < (b)) ? (a) : (b)

// uncomment to enable debug mode
// #define DEBUG 0

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
int **new2d (int width, int height)
{
    int **dp = new int *[width];
    int *dp0 = new int [width * height];
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
    
    // table for storing optimal substructure answers
    //int dp[m+1][n+1] = {0};
    int **dp = new2d (m+1, n+1);
    memset (dp[0], 0, (m+1) * (n+1));

    // intialising the table
    #pragma omp parallel for
    for (i = 0; i <= m; i++) {
        dp[i][0] = i * pgap;
    }
    #pragma omp parallel for
    for (i = 0; i <= n; i++) {
        dp[0][i] = i * pgap;
    }

    #ifdef DEBUG
        for (i = 0; i <= m; i++) 
            for (j = 0; j <= n; j++) 
                // Prints ' ' if j != n-1 else prints '\n'           
                std::cout << dp[i][j] << " \n"[j == n]; 
        std::cout << ">>>> \n";
    #endif

    // std::cout << "* m " << m << "\n";
    // std::cout << "* n " << n << "\n";
    // std::cout << "* max " << (MAX(m, n)+1) << "\n";
    
    // calcuting the minimum penalty
    // int row = m;
    // int col = n;
    // int max_mn = MAX(row, col);

    for (i = 0; i <= (m + n -1); i++) {
        // std::cout << "* i " << i << "\n";

        int z1 = i < n ? 0 : i - n + 1;
        int z2 = i < m ? 0 : i - m + 1;

        // std::cout << "* num_threads " << i - z2 - z1 + 1 << "\n";
        #pragma omp parallel for
        for (int j = i - z2; j >= z1; --j) {
            int dpx = j+1;
            int dpy = i-j+1;
            // std::cout << "* i, j " << dpx << "," << dpy << "\n";
            if (x[dpx - 1] == y[dpy - 1]) {
                dp[dpx][dpy] = dp[dpx - 1][dpy - 1];
            }
            else {
                dp[dpx][dpy] = min3(dp[dpx - 1][dpy - 1] + pxy  ,
                                    dp[dpx - 1][dpy]     + pgap ,
                                    dp[dpx][dpy - 1]     + pgap);
            }
        }
    }

    #ifdef DEBUG
        for (i = 0; i <= m; i++) 
        for (j = 0; j <= n; j++) 
            // Prints ' ' if j != n-1 else prints '\n'           
            std::cout << dp[i][j] << " \n"[j == n]; 
    #endif

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
