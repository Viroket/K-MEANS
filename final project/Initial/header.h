#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <time.h>
#include <float.h>

#define MAX_POINTS 10000
#define MAX_CLUSTERS 100
#define NUMBER_OF_BLOCKS 10
#define THREADS_NUMBER = 2

#define TOTAL_SIZE 4
#define CUDA_THREADS 4
#define RANGE 20
#define CUDA_SIZE 5
#define MASTER 0

//every point have x and y and now what claster he need to be in 
typedef struct
{
	int id;
	double x;
	double y;
	double centerX;
	double centerY;
	int clusterId;
}Point;

//my cluster point 
typedef struct
{
	int id;
	double x;
	double y;
}Cluster;


#define FILE_TO_OPEN "D:\\coda final lo makbili\\final project\\Initial\\cluster1.txt"
#define RESULT "D:\\coda final lo makbili\\final project\\Initial\\Result2.txt"
MPI_Datatype clusterDataType();
MPI_Datatype pointDataType();
void createClust(Point* points, Cluster* clusters , int numClusters);
void setCenter(Point* points, Cluster* clusters, int numberOfPoints, int numOfClust);
double DistanceOfPoints(double x1, double x2, double y1, double y2);
bool settingTheClusterGravity(Point* points, Cluster* clusters ,int currentNumOfClusters , int numOfPoints);
void printClusters(Cluster* clusters ,int currentNumOfClusters);
double ClustDiameter(Point* points  , int clusterId , int numOfPoints);
double getingTheClustersQuality(Cluster* clusters , double clusterDiameter, int clusterId ,int currentNumOfClusters);

