using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace TSLCSS
{
    /// <summary>
    /// This is implementation of time series longest common subsequence algorithm
    /// 
    /// Published reports of research 
    /// using this code (or a modified version) should cite the 
    /// article that describes the algorithm:
    /// 
    /// M. Vlachos, M. Hadjieleftheriou, D. Gunopulos, E. Keogh:  
    /// "Indexing Multi-Dimensional Time-Series with Support for Multiple Distance Measures",
    /// In Proc.of 9th SIGKDD, Washington, DC, 2003
    /// </summary>
    public class MyTslcss
    {
        /// <summary>
        /// Time Series matching using the Longest Common Subsequence within region of delta and epsilon.
        /// </summary>
        /// <param name="dataA">One time series data</param>
        /// <param name="dataB">The second time series data</param>
        /// <param name="delta">Time matching region (left & right)</param>
        /// <param name="epsilon">Spatial matching region (up & down)</param>
        /// <param name="transpose">[optional parameter] how much to shift the time series vertically, so as the matching is better visualized</param>
        /// <returns>Similarity</returns>
        public double Match(Vector<double> dataA, Vector<double> dataB, int delta, double epsilon, int transpose = 0  )
        {
            int m = dataA.Count;
            int n = dataB.Count;

            // put the shorter first
            if (n < m)
            {
                Vector<double> temp = dataA;
                dataA = dataB;
                dataB = temp;
                m = dataA.Count;
                n = dataB.Count;
            }

            Matrix<double> lcstable = Matrix<double>.Build.Dense(m + 1, n + 1, 0);

            for (int i = 0; i < m; i++)
            {
                for (int j = i - delta -1; j < i + delta; j++)
                {
                    if (j < 0 || j > n)
                        continue;
                    if ((dataB[j] + epsilon >= dataA[i]) && (dataB[j] - epsilon <= dataA[i]))
                    {
                        lcstable[i + 1, j + 1] = lcstable[i, j] + 1;
                    }
                    else if (lcstable[i, j + 1] > lcstable[i + 1, j])
                    {
                        lcstable[i + 1, j + 1] = lcstable[i, j + 1];
                    }
                    else
                    {
                        lcstable[i + 1, j + 1] = lcstable[i + 1, j];
                    }
                }
            }

            // Get rid of initial conditions
            lcstable = lcstable.SubMatrix(1, lcstable.RowCount - 1, 1, lcstable.ColumnCount - 1);

            // LCS similarity
            double lcs = lcstable.Row(m-1).Max();

            return m > n ? (lcs / m) : (lcs / n);
        }
    }
}
