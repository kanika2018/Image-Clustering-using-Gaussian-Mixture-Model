package gaussianmixturemodel;

public class Matrix {
	// return B = A^T
    public static double[][]transpose(double[][] a) {
    int m = a.length;
    int n = a[0].length;
    double[][] b =new double[n][m];
    for(int i =0; i < m; i++)
    {
        for(int j =0; j < n; j++)
        {
            b[j][i]= a[i][j];
        }   
    }
    return b;
    
    }
    public static double[][]transpose(double[]a) {
    int m = a.length;
    int n = m;
    double[][] b =new double[n][m];
    for(int i =0; i < m; i++)
    {

            b[i][0]= a[i];
          
    }
    return b;
    
    }


    // return c = a + b
    public static double[][]add(double[][] a,double[][] b) {
    int m = a.length;
    int n = a[0].length;
    double[][] c =new double[m][n];
    for(int i =0; i < m; i++)
    {
        for(int j =0; j < n; j++)
        {
           c[i][j]= a[i][j]+ b[i][j]; 
        }
    }    
    return c;
    }

    // return c = a - b
    public static double[][]subtract(double[][] a,double[][] b) {
    int m = a.length;
    int n = a[0].length;
    double[][] c =new double[m][n];
    for(int i =0; i < m; i++)
    {
        for(int j =0; j < n; j++)
        {
            c[i][j]= a[i][j]- b[i][j];
        }
                        
    }    
    return c;
    }
    
    public static double[]subtract(double[]a,double[]b) {
    int m = a.length;
    
    double[]c =new double[m];
    for(int i =0; i < m; i++)
    {

        c[i]= a[i]- b[i];
                  
    }    
    return c;
    }

    // return c = a * b
    public static double[][]multiply(double[][] a,double[][] b) {
    int m1 = a.length;
    int n1 = a[0].length;
    //int m2 = b.length;
    int n2 = b[0].length;
    //if(n1 != m2)thrownewRuntimeException("Illegal matrix dimensions.");
    double[][] c =new double[m1][n2];
    for(int i =0; i < m1; i++)
    {
        for(int j =0; j < n2; j++)
        {
            c[i][j]=0;
            for(int k =0; k < n1; k++)
                        c[i][j]+= a[i][k]* b[k][j];
        }
    
    }
    
    return c;
        }

    // matrix-vector multiplication (y = A * x)
    public static double[]multiply(double[][] a,double[] x) {
    int m = a.length;
    int n = a[0].length;
    //if(x.length != n)thrownewRuntimeException("Illegal matrix dimensions.");
    double[] y =new double[m];
    for(int i =0; i < m; i++)
    {
       for(int j =0; j < n; j++)
                    y[i]+= a[i][j]* x[j]; 
    }
    
    return y;
    }

    public static double[][] multiply(double[][]b,double a)
    {
        int m=b.length;
        int n=b[0].length;
        double[][] y=new double[m][n];
        for(int i=0;i<m;i++)
        {
           for (int j=0;j<n;j++)
                    y[i][j]=b[i][j]*a; 
        }
        
        return y; 
    }
    public static double determinant(double[][]a,int n)
    {
        double det=0.0;
        double sign=1;
        int p=0;
        int q=0;
        
        if(n==1)
        {
            det=a[0][0];
        }
        else
        {
           double b[][]=new double[n-1][n-1];
                for(int x=0;x<n;x++)
                {
                    p=0;q=0;
                    for(int i=1;i<n;i++)
                    {
                        for(int j=0;j<n;j++)
                        {
                            if(j!=x)
                            {
                                b[p][q++]=a[i][j];
                                if(q%(n-1)==0)
                                {
                                    p++;
                                    q=0;
                                }
                            }
                        }
                    }
                    det=det+a[0][x]*determinant(b,n-1)*sign;
                    sign=-sign;
                } 
        }
        
    return det;  
    }
    
            


    // vector-matrix multiplication (y = x^T A)
    public static double[]multiply(double[] x,double[][] a) {
    int m = a.length;
    int n = a[0].length;
    //if(x.length != m)thrownewRuntimeException("Illegal matrix dimensions.");
    double[] y =new double[n];
    for(int j =0; j < n; j++)
    {
        for(int i =0; i < m; i++)
                    y[j]+= a[i][j]* x[i];
    }
    return y;
    }
    
    public static double[][] inverse(double a[][])
    {
        int n =a.length;
	double x[][]=new double[n][n];
	double b[][]=new double[n][n];
	int index[]=new int[n];
	for(int i=0; i<n;++i)
	            b[i][i]=1;
	 
	// Transform the matrix into an upper triangle
	gaussian(a, index);
	 
	// Update the matrix b[i][j] with the ratios stored
	for(int i=0; i<n-1;++i)
	for(int j=i+1; j<n;++j)
	for(int k=0; k<n;++k)
	                    b[index[j]][k]
		-= a[index[j]][i]*b[index[i]][k];
	 
	// Perform backward substitutions
	for(int i=0; i<n;++i)
	{
	            x[n-1][i]= b[index[n-1]][i]/a[index[n-1]][n-1];
	for(int j=n-2; j>=0;--j)
	{
	                x[j][i]= b[index[j]][i];
	for(int k=j+1; k<n;++k)
	{
	                    x[j][i]-= a[index[j]][k]*x[k][i];
	}
	                x[j][i]/= a[index[j]][j];
	}
	}
	return x;
	}
	 
	// Method to carry out the partial-pivoting Gaussian
	// elimination.  Here index[] stores pivoting order.
	 
	public static void gaussian(double a[][], int index[])
	{
	int n =index.length;
	double c[]=new double[n];
	 
	// Initialize the index
	for(int i=0; i<n;++i)
	            index[i]= i;
	 
	// Find the rescaling factors, one from each row
	for(int i=0; i<n;++i)
	{
	double c1 =0;
	for(int j=0; j<n;++j)
	{
	double c0 =Math.abs(a[i][j]);
	if(c0 > c1) c1 = c0;
	}
	            c[i]= c1;
	}
	 
	// Search the pivoting element from each column
	int k =0;
	for(int j=0; j<n-1;++j)
	{	double pi1 =0;
	for(int i=j; i<n;++i)
	{
	double pi0 =Math.abs(a[index[i]][j]);
	                pi0 /= c[index[i]];
	if(pi0 > pi1)
	{
	                    pi1 = pi0;
	                    k = i;
	}
	}
	 
	// Interchange rows according to the pivoting order
	int itmp= index[j];
	            index[j]= index[k];
	            index[k]=itmp;
	for(int i=j+1; i<n;++i)	
	{
	double pj= a[index[i]][j]/a[index[j]][j];
	 
	// Record pivoting ratios below the diagonal
	                a[index[i]][j]=pj;
	 
	// Modify other elements accordingly
	for(int l=j+1; l<n;++l)
	                    a[index[i]][l]-=pj*a[index[j]][l];
}
}
}
}


