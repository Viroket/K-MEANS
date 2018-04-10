#include "header.h"



int main(int argc,char *argv[])
{
	//----------Variables----------//
	FILE* file;
	double quality;
	int maxIter;//maximum iterations for inner loop
	int loop = 0;//counter of iterations
	int numberOfPoints;//the total number of points i got
	Point* points;//the array which we store the points from the file
	int chunkSize;//the size of the sub array master sends to the slaves
	Point* slavePoints;//the sub array which the slaves uses in the slaves block
	Point* masterPoints;//the sub array which the master uses in the master block
	Cluster* clusters;//the array in which the master holdes the clusters
	int clustersNum;//the bumber of clusters in each iteration
	Cluster* slaveClusters;//the sub array of clusters used in the slaves block
	int maxNumOfClust;//maximum number of clusters
	int gravityFlag = 1;//the flag that signals the master and slaves to stop the inner loop
	int qualityFlag = 1;//the flag that signals the master and slaves the program has ended
	int i = 0 , j;
	double diameter;//distance from every other cluster
	double gettingQuality = DBL_MAX;//quality check for every iteration
	double t1, t2;//parameters to check runtime

	int  namelen, numprocs, myid;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);	

	MPI_Get_processor_name(processor_name,&namelen);
	
	MPI_Status status;

	MPI_Datatype MPI_POINT = pointDataType();//creating point data Type
	
	MPI_Datatype MPI_CLUSTER = clusterDataType();//creating cluster data Type
	
	t1 = MPI_Wtime();//start time check

	//----------Read from files----------//
	if(myid == MASTER)
	{
		//MASTER procces read data from file
		file = fopen(FILE_TO_OPEN, "r");//opening the file to read from it
		if (file == NULL)
		{
			printf("File did not opened");
		}
		// reading the first line and getting the data from the file 
		fscanf(file, "%d,", &numberOfPoints); //getting all my points
		fscanf(file, "%d,", &maxNumOfClust); //the maximun number of clusters that are allowed
		fscanf(file, "%d,", &maxIter); //max iterations befor we stop
		fscanf(file, "%lf",  &quality); //the quality neasur 
		clusters = (Cluster*)malloc(sizeof(Cluster) * maxNumOfClust);// array of clusters
		points = (Point*)malloc(sizeof(Point) * numberOfPoints);// array of points
		// init each point
		for (i = 0; i < numberOfPoints; i++)
		{
			fscanf(file, "%d,", &points[i].id);// read id of point
			fscanf(file, "%lf,", &points[i].x);// setting the x of point
			fscanf(file, "%lf", &points[i].y);// setting the y of point

			// setting my point center x and y of the claster -1 for now (and the cluster that he needs to be in)
			points[i].centerX = -1;
			points[i].centerY = -1;
			points[i].clusterId = -1;
		}
		fclose(file);//closing the file after we finished reading from it
	}

	clustersNum = 1;//the starting number of clusters that we want that ites 2 (minimum)
	chunkSize = numberOfPoints/(numprocs-1);//calculating the chunck size by deviding the total points to number of process
	
	//----------Start main program----------//
	if(myid == MASTER)
	{
		while (clustersNum < maxNumOfClust && gettingQuality > quality)// check if we have reached the wanted quality or the maximum number of clusters 
		{
			gettingQuality = 0;
			clustersNum++;
			gravityFlag = 1;
			createClust(points ,clusters ,clustersNum);// create new center point for 'k' cluster and returning the new current num of clusters
			//sending neccesary data to slaves to calculate setCenter(Point* points, Cluster* clusters, int numberOfPoints, int clustersNum)
			for(i=1; i<numprocs; i++)
			{
				MPI_Send(&qualityFlag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&chunkSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&points[(i-1)*chunkSize], chunkSize, MPI_POINT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&clustersNum, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&clusters[0], clustersNum, MPI_CLUSTER, i, 0, MPI_COMM_WORLD);
				MPI_Send(&gravityFlag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
			do
			{
				int nextChunk = 1; //index counter for points array
				while(nextChunk % numprocs != 0)//iteration for slaves responses
				{
					masterPoints = (Point*)malloc(sizeof(Point)*chunkSize);
					MPI_Recv(masterPoints, chunkSize, MPI_POINT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);//gets the chuncks from the slaves to the masterPoints array
					for(j=0; j<chunkSize; j++)
						points[(nextChunk-1)*chunkSize+j] = masterPoints[j];// insert the updated points back to the big points array 
					nextChunk++;
					free(masterPoints);
				}
				if(!settingTheClusterGravity(points, clusters , clustersNum , numberOfPoints))//if the new centrois of the clusters haven't changed signal to the slaves to stop the inner loop
				{
					gravityFlag = 0;
					break;
				}
				else// continue the inner loop of the slaves
				{
					for (i=1;i<numprocs;i++)
					{
						MPI_Send(&clusters[0], clustersNum, MPI_CLUSTER, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
						MPI_Send(&gravityFlag, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
					}
				}
				loop++;
			}while (loop < maxIter);//if points stopped moving or we have reached max iterations

			for (i=1;i<numprocs;i++)//send stop flag for inner loop
			{
				MPI_Send(&clusters[0], clustersNum, MPI_CLUSTER, i, 0, MPI_COMM_WORLD);
				MPI_Send(&gravityFlag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}

			for (j = 0; j < clustersNum; j++) 
			{
				diameter = ClustDiameter(points , j , numberOfPoints);//calculate the distances between every other cluster	
				gettingQuality += getingTheClustersQuality(clusters, diameter, j , clustersNum);//calculates the quality
			}
		}
		qualityFlag = 0;
		for (i=1;i<numprocs;i++)
		{
			MPI_Send(&qualityFlag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);//stop slaves outer loop
		}
		t2 = MPI_Wtime();
		printf("\ntotal time: %lf", t2-t1);
	}
	else//Slaves
	{
		do
		{
			MPI_Recv(&qualityFlag, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
			if(!qualityFlag)//check if quality has reached and stop the slaves if flase
				break;
			MPI_Recv(&chunkSize, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);//recivies the chunck size
			slavePoints = (Point*)malloc(sizeof(Point)*chunkSize);
			MPI_Recv(&slavePoints[0], chunkSize, MPI_POINT, MASTER, 0, MPI_COMM_WORLD, &status);//recivies the chunck of points
			MPI_Recv(&clustersNum, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);//recivies the clusters size
			slaveClusters = (Cluster*)malloc(sizeof(Cluster)*clustersNum);
			do
			{
				MPI_Recv(&slaveClusters[0], clustersNum, MPI_CLUSTER, MASTER, 0, MPI_COMM_WORLD, &status);//recivies the clusters array
				MPI_Recv(&gravityFlag, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
				if(!gravityFlag)//check if points stopped moving and if true stop the inner loop
					break;
				setCenter(slavePoints, slaveClusters , chunkSize , clustersNum);
				MPI_Send(&slavePoints[0], chunkSize, MPI_POINT, MASTER, 0, MPI_COMM_WORLD);//sends the updated chunk of points to the master
			}while(1);
			free(slavePoints);//free memory
			free(slaveClusters);//free memory
		}while(1);
	}

	//----------Read from files----------//
	if(myid == MASTER)
	{
		file = fopen(RESULT, "w"); // open file to rhite inside of it 
		if (file == NULL)
		{
			printf("File did not opened");
		}
		fseek(file, 0, SEEK_SET); // set file point to the start 
		// write data by format
		fputs("Number of best measure clusters ", file);
		fprintf(file, "K = %d \n Quality = %lf\n", clustersNum, gettingQuality);
		fputs("Number of best measure clusters \n", file);

		for (i = 0; i < clustersNum; i++)
		{
			fprintf(file, "x = %lf, y = %lf\n", clusters[i].x, clusters[i].y);
		}
	
		fclose(file); //closeing the file
		free(points);
		free(clusters);
	}
	printf("\n End!");
	MPI_Finalize();
	return(0);
}




///////////////////////////////////////////////initiat the clusters with first points//////////////////////////////////////////////////////////
void createClust(Point* points, Cluster* clusters , int numClusters)
{
	int i ;
	double new_centerVal_x, new_centerVal_y;

	#pragma omp parallel
	for (i = 0; i < numClusters; i++)
	{
		new_centerVal_x = points[i].x;
		new_centerVal_y = points[i].y;
		// set the cluster with the new x,y values
		clusters[i].x = new_centerVal_x;
		clusters[i].y = new_centerVal_y;
		clusters[i].id = i;

		// increase num of current clusters
		
	}
}

////////////////////////////////////////////////determends the closest cluster to every point/////////////////////////////////////////////////////////
void setCenter(Point* points, Cluster* clusters, int numberOfPoints, int clustersNum)
{
	int i, j, minClusterId;
	double minDistance, distance;

	#pragma omp parallel for private(j, distance, minDistance, minClusterId)
	for (i = 0; i < numberOfPoints; i++) // for each point
	{
		minDistance = DistanceOfPoints(points[i].x, points[i].y, clusters[0].x, clusters[0].y);// set first cluster as minDistance
		minClusterId = 0;
		for (j = 1; j < clustersNum; j++)// for each cluster
		{
			distance = DistanceOfPoints(points[i].x, points[i].y, clusters[j].x, clusters[j].y);// get distance of the point from the 'j' center

			if (distance < minDistance)// if distance smaller than minDistance - change to its values
			{
				minDistance = distance;
				minClusterId = j;
			}

		}
		// set closest cluster as point's center
		points[i].centerX = clusters[minClusterId].x;
		points[i].centerY = clusters[minClusterId].y;
		points[i].clusterId = minClusterId;


	}
}

/////////////////////////////////////////////////help function to calculate the distance////////////////////////////////////////////////////////
double DistanceOfPoints(double x1, double y1, double x2, double y2)
{
	return sqrt(((x1 - x2)*(x1 - x2))+((y1 - y2)*(y1 - y2)));
}

/////////////////////////////////////////////////calculates the new centroids////////////////////////////////////////////////////////
bool settingTheClusterGravity(Point* points, Cluster* clusters ,int currentNumOfClusters , int numOfPoints)
{
	int i, j;
	double sumX, sumY, numX, numY;
	double newCenterX, newCenterY;
	int numOfDoneMovingClusters = 0;
	bool pointsMoved = true;


	for (i = 0; i < currentNumOfClusters; i++)// for each cluster
	{
		// init sum computing values
		numX = 0;
		sumX = 0;
		numY = 0;
		sumY = 0;

		for (j = 0; j < numOfPoints; j++)// for avery point
		{

			if (points[j].clusterId == i)// if that point belong to that cluster
			{
				sumY += points[j].y;
				numY++;

				sumX += points[j].x;
				numX++;
			}
		}


		// compute avarege for new center values
		newCenterY = sumY / numY;
		newCenterX = sumX / numX;

		// check if center points didn't move
		if (clusters[i].x == newCenterX && clusters[i].y == newCenterY)
		{
			numOfDoneMovingClusters++;
		}

		// set new center points for cluster
		clusters[i].y = newCenterY;
		clusters[i].x = newCenterX;

	}

	// if all cluster's center points stop moving - inform the main loop
	if (numOfDoneMovingClusters == currentNumOfClusters)
		pointsMoved = false;


	return pointsMoved;

}

/////////////////////////////////////////////////calculates the diameter///////////////////////////////////////////////////////////////////////////////////
double ClustDiameter(Point* points , int clusterId , int numOfPoints)
{
	double maxDis=0, tempDis;
	int j;

	for (int i = 0; i < numOfPoints-2; i++) // for avey point
	{
		if (points[i].clusterId == clusterId) // if the point is belong to the cluster
		{
			for (j = i+1; j < numOfPoints-1; j++) // for avery one of the other points
			{
				if (points[j].clusterId == clusterId)// if my new point belong to cluster that we want
				{
					tempDis = DistanceOfPoints(points[i].x, points[i].y, points[j].x, points[j].y);// get the distance between the point and the seacond point
					
					if (tempDis > maxDis) // if the new distanse is bigger then ower distance we will change the max distance to save it
					{
						maxDis = tempDis;
					}
				}
			}
		}
	}
	return maxDis;
}

//////////////////////////////////////////////////calculates the qualityS///////////////////////////////////////////////////////
double getingTheClustersQuality(Cluster* clusters , double clusterDiameter, int clusterId ,int numOfClusters)
{
	int i;
	double clusterDistence, clusterPartQuality = 0;

	#pragma omp parallel for
	for (i = 0; i < numOfClusters; i++) // for avery cluster 
	{
		if (clusters[i].id != clusterId) // looking for current computed cluster
		{
			clusterDistence = DistanceOfPoints(clusters[clusterId].x, clusters[clusterId].y, clusters[i].x, clusters[i].y);// get distance between neighbor one and the current one
			
			clusterPartQuality += (clusterDiameter/clusterDistence); // adding it to current clusters part
		}
	}

	return clusterPartQuality;
}

MPI_Datatype clusterDataType()
{
	Cluster c;
	int array_of_blocklengths[3] = {1, 1, 1};
	MPI_Aint array_of_displacements[3];
	MPI_Datatype array_of_types[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Datatype MPI_CLUSTER;

	array_of_displacements[0] = (char*)&(c.id) - (char*)&(c);
	array_of_displacements[1] = (char*)&(c.x) - (char*)&(c);
	array_of_displacements[2] = (char*)&(c.y) - (char*)&(c);
	
	MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_CLUSTER);
	MPI_Type_commit(&MPI_CLUSTER);
	return MPI_CLUSTER;
}

MPI_Datatype pointDataType()
{
	Point p;
	int array_of_blocklengths[6] = {1, 1, 1, 1, 1, 1};
	MPI_Aint array_of_displacements[6];
	MPI_Datatype array_of_types[6] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
	MPI_Datatype MPI_POINT;

	array_of_displacements[0] = (char*)&(p.id) - (char*)&(p);
	array_of_displacements[1] = (char*)&(p.x) - (char*)&(p);
	array_of_displacements[2] = (char*)&(p.y) - (char*)&(p);
	array_of_displacements[3] = (char*)&(p.centerX) - (char*)&(p);
	array_of_displacements[4] = (char*)&(p.centerY) - (char*)&(p);
	array_of_displacements[5] = (char*)&(p.clusterId) - (char*)&(p);
	
	MPI_Type_create_struct(6, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_POINT);
	MPI_Type_commit(&MPI_POINT);
	return MPI_POINT;
}