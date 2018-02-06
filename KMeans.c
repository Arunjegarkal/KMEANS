#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/*Function for random number generation*/
float* randomGenerator(int total_datapoints) {
    float *data = (float *)malloc(sizeof(float) * total_datapoints);
    for (int i = 0; i < total_datapoints; i++) 
    {
        //data[i] = ((rand()%1000/ (float)(RAND_MAX)) * 100);
        data[i] =((rand()%1000)+(rand() / (float)RAND_MAX));
    }
    return data;
}

/*Function for calculating distance */
float findDistance(float *datapoint1, float *datapoint2, int dim) {
    float distance = 0.0;
    for (int i = 0; i<dim; i++) {
        float temp = datapoint1[i] - datapoint2[i];
        temp = temp * temp;
        distance += temp;
    }
    return distance;
}
/*Function to display the centroids*/
void centroids(float * centroids, int cluster_num, int dim) 
{
    printf("Centroid is ");
    for (int i = 0; i < cluster_num * dim; i++) 
    {
        printf("[%f]", centroids[i]);
    }
    printf("\n");
}

/*Function to search datapoints in a cluster*/
void searchDatapoint(float *query, float *datapoints, int dim){
    float shortestDistance=0.0,dest=0.0;
    //printf("\n\ntest");
    float *nearest_datapoint;
    int near = 0;
    for(int i=0; i< 100; i+=dim)
    {

        dest= findDistance(query,&datapoints[i], dim);
        //printf("\n distance=%f\n",dest);
        if(i==0)
        {
            shortestDistance=dest;

                
                near =i;
                


        }
        if(shortestDistance > dest)
        {
           // printf("shortestDistance > dest /t");
            shortestDistance=dest;
            near =i;

        }


    }
    printf("Nearest DataPoint ");
    for(int k=near;k<dim+near;k++ )
    printf("[%f]",datapoints[k]);
    printf("\n\n\nShortest Distance is %f\n",shortestDistance);
}



int main(int argc, char **argv) {
    int no_of_datapoints=10000,dim=15;
    MPI_Init(NULL, NULL); //initializing the MPI process
    int cluster_rank, no_of_clusters;
    MPI_Comm_rank(MPI_COMM_WORLD, &cluster_rank); //rank of the MPI process
    MPI_Comm_size(MPI_COMM_WORLD, &no_of_clusters); //size of the MPI process
    int nData = no_of_datapoints / no_of_clusters; //data points per cluster
    int total_of_cluster = no_of_clusters -1;// total number 0f clusters.
    /*Dynamic Memory allocation*/
    float* Data = malloc(dim * nData * sizeof(float));
    float* cluster_size = malloc(dim * total_of_cluster * sizeof(float));
    int* cluster = malloc(total_of_cluster * sizeof(int));
    float* initial_centroids = malloc(dim * total_of_cluster * sizeof(float));
    int* cluster_process = malloc(nData * sizeof(int));


    float* initial_dataset = NULL; //dataset of each cluster
    float* total_data = NULL; //total data points in cluster
    int* total_clusters = NULL; //total no of clusters
    int* new_centroids = NULL;

    if (cluster_rank == 0) 
    {
        initial_dataset = randomGenerator(dim * nData * no_of_clusters);
        for (int i = 0; i < dim * total_of_cluster; i++) {
            initial_centroids[i] = initial_dataset[i];
        }
        centroids(initial_centroids, total_of_cluster, dim);
        total_data = malloc(dim * total_of_cluster * sizeof(float));
        total_clusters = malloc(total_of_cluster * sizeof(int));
        new_centroids = malloc(no_of_clusters * nData * sizeof(int));
    }
    /*Distibute data points to all process.*/
    MPI_Scatter(initial_dataset, dim *nData, MPI_FLOAT, Data,dim*nData, MPI_FLOAT, 0, MPI_COMM_WORLD);


    float min_dist = 1.0;

    while (min_dist > 0.00001) { 

        /*Broadcast the current cluster centroids to all processes.*/
        MPI_Bcast(initial_centroids, dim * total_of_cluster, MPI_FLOAT, 0, MPI_COMM_WORLD);


        for (int i = 0; i < total_of_cluster; i++) {
            cluster[i] = 0;
        }

        for (int i = 0; i < dim * total_of_cluster; i++) {
            cluster_size[i] = 0.0;
        }

        float* current_data_point = Data;


        /*Calculate the nearest centroid and update the sum of the centroid*/
        for (int i = 0; i < nData; i++) 
        {
            int clusterIndex = 0;
            float dist = findDistance(current_data_point, initial_centroids, dim);
            float* nextCentroid = initial_centroids + dim;
            for (int c = 1; c < total_of_cluster; c++, nextCentroid += dim) 
            {
                float new_dist = findDistance(current_data_point, nextCentroid, dim);
                if (dist > new_dist) 
                {
                    clusterIndex = c;
                    dist = new_dist;
                }
            }

            cluster[clusterIndex]++;
            //update the sum vector for centroid
            float *sum_of_cluster = &cluster_size[clusterIndex*dim];
            for (int i = 0; i<dim; i++) {
                sum_of_cluster[i] += current_data_point[i];
            }

            current_data_point += dim;
        }

        MPI_Reduce(cluster_size, total_data, dim * total_of_cluster, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(cluster, total_clusters, total_of_cluster, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if (cluster_rank == 0) 
        {
                /*Find the new centroids*/
            for (int i = 0; i<total_of_cluster; i++) 
            {
                for (int j = 0; j<dim; j++) {
                    int tmp = dim*i + j;
                    total_data[tmp] = total_data[tmp] / total_clusters[i];
                }
            }
        
            min_dist = findDistance(total_data, initial_centroids, dim * total_of_cluster);
           for (int i = 0; i<dim * total_of_cluster; i++) {
                initial_centroids[i] = total_data[i];
            }
            centroids(initial_centroids, total_of_cluster, dim);
        }
        MPI_Bcast(&min_dist, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

// centroid are fixed .
    float* data_point = Data;
    for (int i = 0; i < nData; i++, data_point += dim) {
        int clusterIndex = 0;
        float new_dist = findDistance(data_point, initial_centroids, dim);
        float* nextCentroid = initial_centroids + dim;
        for (int c = 1; c < total_of_cluster; c++, nextCentroid += dim) {
            float dist = findDistance(data_point, nextCentroid, dim);
            if (dist < new_dist) {
                clusterIndex = c;
                new_dist = dist;
            }
        }
        cluster_process[i] = clusterIndex;
    }

    MPI_Gather(cluster_process, nData, MPI_INT,
               new_centroids, nData, MPI_INT, 0, MPI_COMM_WORLD);

    if (cluster_rank == 0) {
        float* data = initial_dataset;
        for (int i = 0; i < no_of_clusters * nData; i++) {
          data += dim;
        }
    }


    int i,j;
    float datapoint[dim];

    for(i=0; i< dim; i++)
    {
        datapoint[i]= rand() % 100 + 1;

    }

    printf("\nSearch DataPoint  ");
    for(i=0; i< dim; i++)
    {
       // datapoint[i]= rand() % 10 + 1;
        printf("[%f]",datapoint[i]);
    }
    printf("\n");

    searchDatapoint(datapoint,data_point,dim);

    MPI_Finalize();

}
