package gaussianmixturemodel;

import java.util.Scanner;

public class Main {
	public static void main(String[] args) {
	      
        int no_of_data_points;
        int no_of_clusters;
        int no_of_bands;
        double threshold;
        Scanner scan=new Scanner(System.in);
        System.out.println("Enter the number of data points: ");
        no_of_data_points=scan.nextInt();
        System.out.println("Enter the number of bands: ");
        no_of_bands=scan.nextInt();
        double[][] input_matrix=new double[no_of_data_points][no_of_bands];
        System.out.println("Enter input matrix: ");
        for(int counter1=0;counter1<no_of_data_points;counter1++)
        {
            for(int counter2=0;counter2<no_of_bands;counter2++)
            {
                input_matrix[counter1][counter2]=scan.nextDouble();
            }
        }
        
        System.out.println("Enter the number of clusters: ");
        no_of_clusters=scan.nextInt();
        System.out.println("Enter the threshold value: ");
        threshold=scan.nextDouble();
                
        double[][] responsibility_matrix=GaussianMixtureModel.fgaussianMixtureModel(no_of_data_points,no_of_bands,no_of_clusters,input_matrix,threshold);
        
        //printing the final responsibility matrix
        
        System.out.println("Responsibility matrix");
        for(int counter1=0;counter1<no_of_data_points;counter1++)
        {
            for(int counter2=0;counter2<no_of_clusters;counter2++)
            {
                System.out.print(responsibility_matrix[counter1][counter2] + "  ");
            }
            System.out.println();
        }
        System.out.println();
        //System.out.println("Showing sum of each row is 1");
        double sum=0;
         for(int counter=0;counter<no_of_data_points;counter++)
             
        {
           
            for(int counter2=0;counter2<no_of_clusters;counter2++)
            {
                if(responsibility_matrix[counter][counter2]>0.9)
                {
                	System.out.print(counter2);
                }
                
            }
            System.out.println();
            
        }
         System.out.println();

}
}
