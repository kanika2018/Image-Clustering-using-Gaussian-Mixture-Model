package gaussianmixturemodel;

public class GaussianMixtureModel {
	/**
     * @param p_no_of_data_points
     * @param p_no_of_bands
     * @param p_input_matrix
     * @param p_threshold
     * @param p_no_of_clusters
     * @return 
     */
    public static double[][] fgaussianMixtureModel(int p_no_of_data_points,int p_no_of_bands,int p_no_of_clusters,double[][]p_input_matrix,double p_threshold) {
        
        // TODO code application logic here
        int no_of_data_points=p_no_of_data_points;
        int no_of_bands=p_no_of_bands;
        int no_of_clusters=p_no_of_clusters;
        double[][] input_matrix=p_input_matrix;
        double threshold=p_threshold;
        double[][] responsibility_matrix=new double[no_of_data_points][no_of_clusters];
        double[] pi_k=new double[no_of_clusters];
        double[][][] covariance_matrix=new double[no_of_clusters][no_of_bands][no_of_bands];
        double[][] mean_matrix=new double[no_of_clusters][no_of_bands];
        
        
//        System.out.println("printing input matrix");
//        for(int counter1=0;counter1<no_of_data_points;counter1++)
//        {
//            for(int counter2=0;counter2<no_of_bands;counter2++)
//            {
//                System.out.print(input_matrix[counter1][counter2] + "  ");
//            }
//            System.out.println();
//        }
//         System.out.println();
        

        
        //initializing parameters
        
        //initializing pi_k[] i.e mixing_coefficient_of_clusters
        double one_by_k;
        one_by_k=(double)1/no_of_clusters;
        
        for(int counter=0;counter<no_of_clusters;counter++)
        {
            pi_k[counter]=one_by_k;
        }
        
//        System.out.println("Initialization Step");
//        System.out.println("printing pi_k matrix");
//        for(int counter=0;counter<no_of_clusters;counter++)
//        {
//            System.out.print(mixing_coefficient_of_clusters[counter] + " ");
//        }
//        System.out.println();
//        System.out.println();
        
        
        //initializing covariance_matrix[k][b][b]
        double[][] dataset_covariance_matrix=new double[no_of_bands][no_of_bands];
        double[][] deviation_score_matrix=new double[no_of_data_points][no_of_clusters];
        double[][] unity_matrix= new double[no_of_data_points][no_of_data_points];
        for(int counter1=0;counter1<no_of_data_points;counter1++)
        {
            for(int counter2=0;counter2<no_of_data_points;counter2++)
            {
                unity_matrix[counter1][counter2]=1;
            }
        }
        double one_by_n;
        one_by_n=(double)1/no_of_data_points;
        deviation_score_matrix=Matrix.subtract(input_matrix,Matrix.multiply(Matrix.multiply(unity_matrix, input_matrix),one_by_n));
        dataset_covariance_matrix=Matrix.multiply(Matrix.multiply(Matrix.transpose(deviation_score_matrix),deviation_score_matrix ),one_by_n);
        for(int counter1=0;counter1<no_of_clusters;counter1++)
        {
            for(int counter2=0;counter2<no_of_bands;counter2++)
            {
                System.arraycopy(dataset_covariance_matrix[counter2], 0, covariance_matrix[counter1][counter2], 0, no_of_bands);
            }
        }
        
        double mean[]=new double[no_of_bands];
        for(int i=0;i<no_of_data_points;i++)
        {
        	for(int j=0;j<no_of_bands;j++)
        	{
        		mean[j]=mean[j]+input_matrix[i][j];
        	}
        }
        for(int j=0;j<no_of_bands;j++)
        {
        	mean[j]=mean[j]*one_by_n;
        }
        double[][]z=new double[no_of_data_points][no_of_bands];
        for(int i=0;i<no_of_data_points;i++)
        {
        	for(int j=0;j<no_of_bands;j++)
        	{
        		z[i][j]=input_matrix[i][j]-mean[j];
        	}
        }
        double[][]c=new double[no_of_bands][no_of_bands];
        c=Matrix.multiply(Matrix.transpose(z), z);
        double one_by_n_minus_one=(double)1/(no_of_data_points);
        c=Matrix.multiply(c, one_by_n_minus_one);
        for(int i=0;i<no_of_bands;i++)
        {
        	for(int j=0;j<no_of_bands;j++)
        	{
        		System.out.print(c[i][j]+" ");
        	}
        	System.out.println();
        	
        }
//        for(int counter1=0;counter1<no_of_clusters;counter1++)
//        {
//            for(int counter2=0;counter2<no_of_bands;counter2++)
//            {
//                System.arraycopy(c[counter2], 0, covariance_matrix[counter1][counter2], 0, no_of_bands);
//            }
//        }
        
        System.out.println("printing covariance matrix for each cluster");
        for(int counter1=0;counter1<no_of_clusters;counter1++)
        {
            for(int counter2=0;counter2<no_of_bands;counter2++)
            {
                for(int counter3=0;counter3<no_of_bands;counter3++)
                {
                    System.out.print(covariance_matrix[counter1][counter2][counter3] + "  ");
                }
                System.out.println();
            }
            System.out.println();
        }
        System.out.println();
        
        
        
        //initializing mu_k
//        double Dgap=(double)no_of_data_points/(double)no_of_clusters;
//        int gap=(int) Math.round(Dgap);
//        
//        int nth_data_point=0;
//        for(int counter1=0;counter1<no_of_clusters;counter1++)
//        {
//            System.arraycopy(input_matrix[nth_data_point], 0, mean_matrix[counter1], 0, no_of_bands);
//            nth_data_point+=gap;
//        }
        
//        System.out.println("printing mu_k matrix");
//        for(int counter1=0;counter1<no_of_clusters;counter1++)
//        {
//            for(int counter2=0;counter2<no_of_bands;counter2++)
//            {
//                System.out.print(mean_matrix[counter1][counter2] + "  ");
//            }
//            System.out.println();
//        }
//         System.out.println();
         
        //inital logLikelihood calculation
        Kmeans KM=new Kmeans();
        KM.clustering(input_matrix,no_of_clusters,-1,null);
        mean_matrix=KM.getCentroids();
        
        double loglikelihood=0.0;
        double probOfXn=0.0;
        double[] temp_xn=new double[no_of_bands];
        double[] temp_mu_k=new double[no_of_bands];
        double[][] temp_covariance_matrix=new double[no_of_bands][no_of_bands];
        double gaussian_value;
        double b_by_2=(double)no_of_bands/2;
        double[] temp_matrix_computation;
        double temp_computed_value;
        double new_loglikelihood=0.0;
        double determinant;
        

        for(int counter1=0;counter1<no_of_data_points;counter1++ )
            {
                for(int counter2=0;counter2<no_of_clusters;counter2++)
                {
                    for(int counter3=0;counter3<no_of_bands;counter3++)
                    {
                        
                        temp_xn[counter3]=input_matrix[counter1][counter3];
                        
                        temp_mu_k[counter3]=mean_matrix[counter2][counter3];
                        
                        System.arraycopy(covariance_matrix[counter2][counter3], 0, temp_covariance_matrix[counter3], 0, no_of_bands);
                        
                        
                    }
                    
                    temp_matrix_computation=Matrix.multiply(Matrix.subtract(temp_xn, temp_mu_k), Matrix.multiply(Matrix.inverse(temp_covariance_matrix), Matrix.transpose(Matrix.subtract(temp_xn, temp_mu_k))));
                    temp_computed_value=temp_matrix_computation[0];
                    gaussian_value=Math.pow(2*3.14, -b_by_2)*Math.pow(Matrix.determinant(temp_covariance_matrix,no_of_bands), -0.5)* Math.exp(-0.5*temp_computed_value);
                   
                    probOfXn=probOfXn+ pi_k[counter2]*gaussian_value;
                    
                }
                
                loglikelihood=loglikelihood+Math.log(probOfXn)/Math.log(2.0);
                
                
            }
            
            System.out.println("the log likelihood is" + loglikelihood );
            
        
        
        double change_in_loglikelihood=1.0;
        
        while(change_in_loglikelihood>threshold)
        {

            //Estimation Step
            //calculating responsibility matrix
            for(int counter1=0;counter1<no_of_data_points;counter1++)
            {
                probOfXn=0.0;
                for(int counter2=0;counter2<no_of_clusters;counter2++)
                {
                    for(int counter3=0;counter3<no_of_bands;counter3++)
                    {
                        
                        temp_xn[counter3]=input_matrix[counter1][counter3];
                        
                        temp_mu_k[counter3]=mean_matrix[counter2][counter3];
                        System.arraycopy(covariance_matrix[counter2][counter3], 0, temp_covariance_matrix[counter3], 0, no_of_bands);
                        
                        
                    }
                    Matrix.subtract(temp_xn,temp_mu_k);
                    

                    temp_matrix_computation=Matrix.multiply(Matrix.subtract(temp_xn, temp_mu_k), Matrix.multiply(Matrix.inverse(temp_covariance_matrix), Matrix.transpose(Matrix.subtract(temp_xn, temp_mu_k))));
                    temp_computed_value=temp_matrix_computation[0];
                    gaussian_value=Math.pow(2*3.14, -b_by_2)*Math.pow(Matrix.determinant(temp_covariance_matrix,no_of_bands), -0.5)* Math.exp(-0.5*temp_computed_value);
                    
                    responsibility_matrix[counter1][counter2]=pi_k[counter2]*gaussian_value;
                    probOfXn=probOfXn+responsibility_matrix[counter1][counter2];
                }
                for(int counter4=0;counter4<no_of_clusters;counter4++)
                {
                    responsibility_matrix[counter1][counter4]=responsibility_matrix[counter1][counter4]/probOfXn;
                }
                
            }
            
//            System.out.println("printing responsibility matrix");
//            for(int counter1=0;counter1<no_of_data_points;counter1++)
//            {
//                for(int counter2=0;counter2<no_of_clusters;counter2++)
//                {
//                    System.out.print(responsibility_matrix[counter1][counter2] + "  ");
//                }
//                System.out.println();
//            }
//             System.out.println();
            
            //Maximization Step
            //calculating changed mu_k,pi_k,sigma_k
            
            //calculating N_k i.e sum of responsibility of each cluster
            
            double[] cluster_responsibility_vector=new double[no_of_clusters];
            for(int counter1=0;counter1<no_of_clusters;counter1++)
            {
                for(int counter2=0;counter2<no_of_data_points;counter2++)
                {
                    cluster_responsibility_vector[counter1]+= responsibility_matrix[counter2][counter1];
                }
            }
            
//            System.out.println("printing N_k");
//            for(int counter4=0;counter4<no_of_clusters;counter4++)
//            {
//                System.out.print(cluster_responsibility_vector[counter4]+ " ");
//            }
//            System.out.println();
//            System.out.println();
            
            //calculating new_pi_k
            
            for(int counter1=0;counter1<no_of_clusters;counter1++)
            {
                pi_k[counter1]=cluster_responsibility_vector[counter1]/no_of_data_points;
            }
            
//            System.out.println("printing new pi_k");
//            for(int counter4=0;counter4<no_of_clusters;counter4++)
//            {
//                System.out.print(pi_k[counter4]+ " ");
//            }
//            System.out.println();

            //calculating new_mu_k
            for(int counter1=0;counter1<no_of_clusters;counter1++)
            {
                for(int counter2=0;counter2<no_of_bands;counter2++)
                {
                    mean_matrix[counter1][counter2]=0;
                    for(int counter3=0;counter3<no_of_data_points;counter3++)
                    {
                        for(int counter4=0;counter4<no_of_bands;counter4++)
                        {
                            temp_xn[counter4]=input_matrix[counter3][counter4];
                        }
                        
                    
                        mean_matrix[counter1][counter2]+=responsibility_matrix[counter3][counter1]*temp_xn[counter2];
                        
                    }
                }
                
                for(int counter5=0;counter5<no_of_bands;counter5++)
                {
                    mean_matrix[counter1][counter5]=mean_matrix[counter1][counter5]/cluster_responsibility_vector[counter1];
                }
            }
            
//            System.out.println();
//            System.out.println("printing new mu_k matrix");
//            for(int counter1=0;counter1<no_of_clusters;counter1++)
//            {
//                for(int counter2=0;counter2<no_of_bands;counter2++)
//                {
//                    System.out.print(mean_matrix[counter1][counter2]+ " ");
//                }
//                System.out.println();
//            }
//            System.out.println();
            
            //calculating new sigma_k
            
            double one_by_N_k;
            double temp_subtract_matrix[][]=new double[no_of_bands][no_of_bands];
            double temp_subtract_vector[];
            for(int counter1=0;counter1<no_of_clusters;counter1++)
            {
                double[][] kth_covariance_matrix=new double[no_of_bands][no_of_bands];
                for(int counter2=0;counter2<no_of_data_points;counter2++)
                {
                    for(int counter3=0;counter3<no_of_bands;counter3++)
                    {
                        temp_xn[counter3]=input_matrix[counter2][counter3];
                        temp_mu_k[counter3]=mean_matrix[counter1][counter3];
                    }
                    temp_subtract_vector=Matrix.subtract(temp_xn, temp_mu_k);
                    for(int counter4=0;counter4<no_of_bands;counter4++)
                    {
                        temp_subtract_matrix[0][counter4]=temp_subtract_vector[counter4];
                    }
                    double[][] temp_kth_covariance_matrix=Matrix.multiply(Matrix.multiply(Matrix.transpose(Matrix.subtract(temp_xn,temp_mu_k)),temp_subtract_matrix),responsibility_matrix[counter2][counter1]);
                    
                    kth_covariance_matrix=Matrix.add(kth_covariance_matrix,temp_kth_covariance_matrix);
                    
                }
                one_by_N_k=1/cluster_responsibility_vector[counter1];
                kth_covariance_matrix=Matrix.multiply(kth_covariance_matrix,one_by_N_k);
                for(int counter5=0;counter5<no_of_bands;counter5++)
                {
                    System.arraycopy(kth_covariance_matrix[counter5], 0, covariance_matrix[counter1][counter5], 0, no_of_bands);
//                    for(int counter6=0;counter6<no_of_bands;counter6++)
//                        covariance_matrix[counter1][counter5][counter6]=kth_covariance_matrix[counter5][counter6];
                }
            }
            
//            System.out.println("printing new covariance matrix for each cluster");
//            for(int counter1=0;counter1<no_of_clusters;counter1++)
//            {
//                for(int counter2=0;counter2<no_of_bands;counter2++)
//                {
//                    for(int counter3=0;counter3<no_of_bands;counter3++)
//                    {
//                        System.out.print(covariance_matrix[counter1][counter2][counter3] + "  ");
//                    }
//                    System.out.println();
//                }
//                System.out.println();
//            }
//            System.out.println();
            
            //calculating new loglikelihood
            gaussian_value=0.0;
            probOfXn=0.0;
            new_loglikelihood=0;
            for(int counter1=0;counter1<no_of_data_points;counter1++ )
            {
                for(int counter2=0;counter2<no_of_clusters;counter2++)
                {
                    for(int counter3=0;counter3<no_of_bands;counter3++)
                    {
                        
                        temp_xn[counter3]=input_matrix[counter1][counter3];
                        
                        temp_mu_k[counter3]=mean_matrix[counter2][counter3];
                        
                        System.arraycopy(covariance_matrix[counter2][counter3], 0, temp_covariance_matrix[counter3], 0, no_of_bands);
                        
                        
                    }
//                    double[][] temp_inverse=Matrix.inverse(temp_covariance_matrix);
//                    for(int i=0;i<no_of_bands;i++)
//                    {
//                        for(int j=0;j<no_of_bands;j++)
//                        {
//                            System.out.print(temp_inverse[i][j]+" ");
//                        }
//                        System.out.println();
//                    }
                    //System.out.println();
                    temp_matrix_computation=Matrix.multiply(Matrix.subtract(temp_xn, temp_mu_k), Matrix.multiply(Matrix.inverse(temp_covariance_matrix), Matrix.transpose(Matrix.subtract(temp_xn, temp_mu_k))));
//                    System.out.println("temp_matrix");
//                    for(int i=0;i<no_of_bands;i++)
//                    {
//                        System.out.print(temp_matrix_computation[i]+" ");
//
//                    }
//                    System.out.println();
                    temp_computed_value=temp_matrix_computation[0];
                    //System.out.println("temp: "+temp_matrix_computation[0]);
                    //System.out.println("det: "+Matrix.determinant(temp_covariance_matrix,no_of_bands));
                    //added math.absolute
                    determinant=Matrix.determinant(temp_covariance_matrix,no_of_bands);
//                    if(determinant==0)
//                    {
//                        determinant=0.01;
//                    }
                    gaussian_value=Math.pow(2*3.14, - b_by_2)*Math.pow(determinant,-0.5)* Math.exp(-0.5*temp_computed_value);
                    
                    //System.out.println("gaussian: "+ gaussian_value);
                    probOfXn=probOfXn+ pi_k[counter2]*gaussian_value;
                    //System.out.println("prob: "+ probOfXn);
                    
                    
                }
                
                //System.out.println("prob: "+probOfXn);
                //System.out.println("log of prob"+Math.log(probOfXn));
                new_loglikelihood=new_loglikelihood+Math.log(probOfXn)/Math.log(2.0);
                
            }
            System.out.println("loglikelihood is :"+ loglikelihood);
            System.out.println("new loglikelihood is :"+ new_loglikelihood);
            change_in_loglikelihood=Math.abs(new_loglikelihood-loglikelihood);
            loglikelihood=new_loglikelihood;
            
            
            System.out.println("change in loglikelihood is :"+ change_in_loglikelihood);
  
        }

        return responsibility_matrix;
 
    }
    
}


