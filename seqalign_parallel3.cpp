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
    if (a <= b && a <= c) {
        return a;    
    } else if (b <= a && b <= c) {
        return b;
    } else {
        return c;    
    }
}

/******************************************************************************/
/* Do not change any lines above here.            */
/* All of your changes should be below this line. */
/******************************************************************************/
#include "omp.h"
#include "assert.h"
#include <bits/stdc++.h>

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

// uncomment to enable debug mode
// #define DEBUG 0

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
    for (i = 0; i <= m; i++)
    {
        dp[i][0] = i * pgap;
    }
    for (i = 0; i <= n; i++)
    {
        dp[0][i] = i * pgap;
    }

    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            dp[i][j] = -1;
        }
    }

    #ifdef DEBUG
        for (i = 0; i <= m; i++) 
            for (j = 0; j <= n; j++) 
                // Prints ' ' if j != n-1 else prints '\n'           
                std::cout << dp[i][j] << " \n"[j == n]; 
        std::cout << ">>>> \n";
    #endif

    vector< pair<int,int> > vect;

    for (i = 0; i <= (m + n -1); i++) {
        // std::cout << "* i " << i << "\n";
        if (i <= n) {
            for(j = i; j > 0; j--) {
                int dpx = i+1-j;
                if (dpx <= m) {
                    // std::cout << "* i, j " << dpx << "," << j << " if \n";
                    vect.push_back( make_pair(dpx, j) ); 
                }
            }
        } else {
            for(j = n; j >= i-m+1; j--) {
                if (j > 0) {
                    int dpx = i+1-j;
                    // std::cout << "* i, j " << dpx << "," << j << " else\n";
                    vect.push_back( make_pair(dpx, j) ); 
                }
            }
        }
    }
    
    int n_threads = 4;
    #pragma omp parallel num_threads(n_threads) private(i) shared(dp)
    {
    int id = omp_get_thread_num();

    for (i = id; i < vect.size(); i += n_threads) {
        pair<int,int> dp_xy = vect.at(i);
        int left_x    = dp_xy.first-1, 
            left_y    = dp_xy.second, 
            up_x      = dp_xy.first, 
            up_y      = dp_xy.second-1, 
            left_up_x = dp_xy.first-1, 
            left_up_y = dp_xy.second-1;

        #pragma omp critical
        cout << "i = " << i << "; thread: " << omp_get_thread_num() << endl;

        #pragma omp critical
        cout << "*waiting (i,j) (" << dp_xy.first << "," << dp_xy.second << ")" << endl;
        // wait for synchronization
        // #pragma omp flush(dp)
        int flag = 0;
        while (dp[left_x][left_y] == -1 || dp[up_x][up_y] == -1 || dp[left_up_x][left_up_y] == -1) {
            if (flag == 0) {
                #pragma omp critical
                cout << "-1 (i,j) (" << dp_xy.first << "," << dp_xy.second << ")" << ", (left)" << dp[left_x][left_y] << ", (up)" << dp[up_x][up_y] << ", (up left)" << dp[left_up_x][left_up_y] << endl;
                flag = 1;
            }
            // cout << "waiting (i,j) (" << dp_xy.first << "," << dp_xy.second << ")" << endl;
            // #pragma omp flush(dp)
            // if (dp[left_x][left_y] == -1) {
            //     // #pragma omp critical
            //     cout << "-1 (i,j) (" << dp_xy.first << "," << dp_xy.second << ")" << endl;
            //     cout << "*-1 (i,j) (" << left_x << "," << left_y << ")" << endl;
            // } else if (dp[up_x][up_y] == -1) {
            //     // #pragma omp critical
            //     cout << "-1 (i,j) (" << dp_xy.first << "," << dp_xy.second << ")" << endl;
            //     cout << "*-1 (i,j) (" << up_x << "," << up_y << ")" << endl;
            // } else if (dp[left_up_x][left_up_y] == -1) {
            //     // #pragma omp critical
            //     cout << "-1 (i,j) (" << dp_xy.first << "," << dp_xy.second << ")" << endl;
            //     cout << "*-1 (i,j) (" << left_up_x << "," << left_up_y << ")" << endl;
            // }
        }

        #pragma omp critical
        cout << "  *running (i,j) (" << dp_xy.first << "," << dp_xy.second << ")" << endl;
        
        // #pragma omp task depend(in:dp[left_up_x][left_up_y],dp[left_x][left_y],dp[up_x][up_y]) depend(out:dp[dp_xy.first][dp_xy.second])
        if (x[left_up_x] == y[left_up_y]) {
            // #pragma omp flush(dp)
            #pragma omp critical
            dp[dp_xy.first][dp_xy.second] = dp[left_up_x][left_up_y];
        } else {
            // #pragma omp flush(dp)
            #pragma omp critical
            dp[dp_xy.first][dp_xy.second] = min3(dp[left_up_x][left_up_y] + pxy  ,
                                                 dp[left_x][left_y]     + pgap ,
                                                 dp[up_x][up_y]     + pgap);
        }

        #pragma omp critical
        cout << "    *finished (i,j) (" << dp_xy.first << "," << dp_xy.second << ")" << endl;
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
