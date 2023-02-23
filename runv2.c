//This is the file that should be run with Open MPI

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

//Include our local files
#include "mpi_structs.c"
//Global Variables
#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1
#define NDIMS 2
//Specify the interval of generating readings(in seconds)
#define INTERVAL 2
//Specify the magnitude to trigger a request for data from neighbours
#define MAGNITUDETHRESHOLD 10.0
//Specify this in KM
#define DISTANCETHRESHOLD 3000.0
//Specify the magnitude difference
#define MAG_DIFF_THRESHOLD 5.0

//Message tags
//Request for message
#define MSG_REQ 1
//Tag for receiving data
#define MSG_DATA 2
//Sensor node reporting to base station
#define MSG_REPORT 3
//Termination
#define MSG_TERM 4
// Base station checking if the node is alive
#define MSG_CHECK 5
//Node is alive
#define MSG_ALIVE 6

//How long we should sleep for
#define SLEEP_MICRO_SEC	100000

//How big our balloon readings array is
#define BALLOONARRAYSIZE 5


//Communicator containing all the seismic processes
MPI_Comm seismic_Comm;
//Communicator for the 2D Cartesian grid
MPI_Comm comm2D;
//Instantiate the data type
MPI_Datatype SeismicData;
//Instantiate the data type
MPI_Datatype SensorReport;

//Variable to indicate whether we are terminating or not
int termination = 0;
//Variable indicating whether we need neighbours data
int needNeighbourData = 0;
//Keeps track of last seismic reading (The one that will be sent to other sensors if they request)
struct SeismicReading lastReading;
//Global array consisting of SeismicReadings from the balloon sensor
struct SeismicReading balloonReadings[BALLOONARRAYSIZE];
int readingIndex = 0;
//Log file to debug
FILE *logFile;
//Rank inside the Cartesian grid
int my_cart_rank;
//Rank inside the sensor comm
int sensor_rank;
//Upper and lower neighbour ranks
int upperNeighbour, lowerNeighbour = -2;
//Left and right neighbour ranks
int leftNeighbour, rightNeighbour = -2;
//The number of iterations we should do
int iterations = 0;
//Size of the MPI_WORLD
int size;
//Flag indicating whether there is a pending report that needs processing or not
int needsProcessing = 0;
//Alert Reported time
struct tm* reported_Time;
//Storing the dimensions
int nrows, ncols;
//Keeping track of which iteration we are at
int iterationCount = 0;

// Checking node if alive code
int check = 100;

// Node alive code
int alive = 200;

// Node failure error code
int error = -1;


//Received report
SeismicReport receivedReport;

//PThread function
void* sensorFunc(void *pArg);

//PThread function for the balloon sensor
void* balloonFunction(void *pArg);

//PThread function for the send and receive for the base station
void* baseSendRecv(void *pArg);

int seismic_sensor(MPI_Comm world_comm);

int base_station();

void writeAdjacentNodeReadings(SeismicReading nodeReading, SeismicReading reportingNodeReading, int rank);

int ifFileExists(const char *filename);

int main(int argc, char **argv)
{

    //Create our own MPI data type for easy sending and receiving
    int my_rank, provided, reorder, ierr;
    //Create the dimensions array, size 2 for 2D
    int dims[NDIMS];
    //Initialise the array of wrap around values, 0 for no and 1 for yes FOR EACH DIMENSION
    int wrap_around[NDIMS];

    //Initialise MPI with Multi Threading
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Split process 0 (our base station) from other processes
    MPI_Comm_split( MPI_COMM_WORLD,my_rank == 0, 0, &seismic_Comm);


    //Ensure that the user has specified 4 arguments
    if(argc == 4){
        //Set the number of rows to the the first argument
        nrows = atoi(argv[1]);
        //Set the number of columns to be the second argument
        ncols = atoi(argv[2]);
        //Set the number of iterations to the third argument
        iterations = atoi(argv[3]);

        //Set the dimensions array
        dims[0] = nrows;
        dims[1] = ncols;
        //Ensure that there are enough processes
        if((nrows * ncols) + 1 != size){ //Number of processes should be (n x m) + 1, extra one for the base station
            if (my_rank == 0){
                printf("ERROR: nrows*ncols + 1 = %d *%d + 1 = %d != %d\n", nrows, ncols, (nrows*ncols) + 1,size);
            }
            MPI_Finalize();
            return 0;

        }

    } else { //The user did not specify 4 arguments
        //Only rank 0 prints
        if (my_rank == 0){
            printf("ERROR: Please specify the right number of arguments, you have specified %d argument(s)\n", argc);
        }
        MPI_Finalize();
        return 0;
    }

    //Create the SeismicReading MPI_DATATYPE
    createSeismicReadingDatatype(&SeismicData);
    //Create the SeismicReport MPI_DATATYPE
    createSeismicReportStruct(&SensorReport, SeismicData);

    //Split the processes
    if(my_rank == 0){ //Root process
        //Call root process function here
        base_station();
    } else{ //Sensor process
        //Create the Cartesian grid for the sensor processes
        //No wrap around
        wrap_around[0] = wrap_around[1] = 0;
        //Set reordering to true so that processes might get different ranks in order to optimise communications
        reorder = 1;
        //Default no error code
        ierr =0;
        //Create the card
        ierr = MPI_Cart_create(seismic_Comm, NDIMS, dims, wrap_around, reorder, &comm2D);
        //If we have an error code
        if(ierr != 0) {
            printf("ERROR[%d] creating CART\n",ierr);
        }
        //Call the Seismic function and pass in the communicators
        seismic_sensor(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

//The base station function
int base_station(){
    //Create two more threads
    pthread_t tid[2];
    //Balloon thread
    pthread_create(&tid[0], 0, balloonFunction, &termination); // Create the thread, pass in the thread ID
    //Send/Receive thread
    pthread_create(&tid[1], 0, baseSendRecv, &termination); // Create the thread, pass in the thread ID

    //Main thread (Writes reports and stuff), write code below
    //Open the file
    logFile = fopen("baseStationLog.txt", "w");
    //Seismic Reading Variable
    SeismicReading currentReading;
    SeismicReading adjacentReading;
    //Seismic Reading variable for balloon sensor
    SeismicReading lastBalloonReading;
    //Flag indicating conclusive or not report
    int conclusive = 0;
    //While we have not received the termination signal
    while(termination == 0){
        //Check if needs processing flag is true, if so we process the report, Add Mutex lock here later
        if(needsProcessing == 1){
            //Set it to false
            needsProcessing = 0;
            //Reset our conclusive flag
            conclusive = 0;
            //Time variables
            time_t s, val;
            struct tm* current_time;
            //Get current time
            s = time(NULL);

            //Uses the current time to create a tm struct
            current_time = localtime(&s);

            //Begin writing the report
            fprintf(logFile, "===========================================================\n");
            //Extract the reporting node's readings into a Seismic Reading
            arrayToReading(&currentReading, receivedReport.reportingNodeReading);

            fprintf(logFile, "Iteration : %d\n", receivedReport.iteration);
            fprintf(logFile, "Logged time : %d-%02d-%02d %02d:%02d:%02d\n", current_time->tm_year + 1900, current_time->tm_mon, current_time->tm_mday, current_time->tm_hour, current_time->tm_min, current_time->tm_sec);
            fprintf(logFile, "Alert reported time : %d-%02d-%02d %02d:%02d:%02d\n", reported_Time->tm_year + 1900, reported_Time->tm_mon, reported_Time->tm_mday, reported_Time->tm_hour, reported_Time->tm_min, reported_Time->tm_sec);

            //Get our current balloon reading
            if(readingIndex == 0){
                //If we are at the first index, get the last element of the array
                lastBalloonReading = balloonReadings[BALLOONARRAYSIZE - 1];
            } else {
                //Else, get the element before our current reading index
                lastBalloonReading = balloonReadings[(readingIndex -1) ];
            }
            //Check whether the report is conclusive or not
            double distance = getDistanceBetweenReadings(currentReading.latitude, currentReading.longitude, currentReading.depth, lastBalloonReading.latitude, lastBalloonReading.longitude, lastBalloonReading.depth);

            if(distance <= DISTANCETHRESHOLD){ //If the reading is within our threshold
                //Check if the magnitude is also within our threshold
                //Calculate the absolute difference between our reporting node's reading and balloon reading
                double magDiff = currentReading.magnitude - lastBalloonReading.magnitude;
                magDiff = fabs(magDiff);
                //Check if the magnitude is greater than our threshold
                if(magDiff <= MAG_DIFF_THRESHOLD){
                    //It is a match
                    conclusive = 1;
                }
            }
            //Check our flag
            if(conclusive == 1){
                fprintf(logFile, "Alert type: Conclusive\n");
            }else {
                fprintf(logFile, "Alert type: Inconclusive\n");
            }
            //Get the coordinates
            int x = getX(receivedReport.reportingNode, ncols);
            int y = getY(receivedReport.reportingNode, nrows);
            fprintf(logFile, "Reporting Node              Seismic Coord                     Magnitude\n");
            fprintf(logFile, "%d (%d,%d)                     (%02.2f, %02.2f)                   %.2f\n", receivedReport.reportingNode, x, y, currentReading.latitude, currentReading.longitude, currentReading.magnitude);

            fprintf(logFile, "Adjacent Nodes              Seismic Coord      Diff(KM)       Magnitude    Diff(Mag)\n");
            //Check if there is an upper neighbour
            if(receivedReport.neighbours[0] >= 0){ //If there is an upper neighbour
                arrayToReading(&adjacentReading, receivedReport.upperNeighbourReading);
                //Write it
                writeAdjacentNodeReadings(adjacentReading, currentReading, receivedReport.neighbours[0]);
             }
            if(receivedReport.neighbours[1] >= 0){ //If there is a lower neighbour
                arrayToReading(&adjacentReading, receivedReport.lowerNeighbourReading);
                //Write it
                writeAdjacentNodeReadings(adjacentReading, currentReading, receivedReport.neighbours[1]);
            }
            if(receivedReport.neighbours[2] >= 0){ //If there is a left neighbour
                arrayToReading(&adjacentReading, receivedReport.leftNeighbourReading);
                //Write it
                writeAdjacentNodeReadings(adjacentReading, currentReading, receivedReport.neighbours[2]);
            }
            if(receivedReport.neighbours[3] >= 0){ //If there is an right neighbour
                arrayToReading(&adjacentReading, receivedReport.rightNeighbourReading);
                //Write it
                writeAdjacentNodeReadings(adjacentReading, currentReading, receivedReport.neighbours[3]);
            }

            fprintf(logFile, "Balloon seismic reporting time : %d-%02d-%02d %02d:%02d:%02d\n", lastBalloonReading.year, lastBalloonReading.month, lastBalloonReading.day, lastBalloonReading.hour, lastBalloonReading.minute, lastBalloonReading.second);
            //printf("Main process loon reading(%.2f,%.2f)\n",lastBalloonReading.latitude, lastBalloonReading.longitude);
            fprintf(logFile, "Balloon seismic reporting Coord : (%.2f,%.2f)\n", lastBalloonReading.latitude, lastBalloonReading.longitude);
            //Find the distance and magnitude difference
            //Get the distance
            distance = getDistanceBetweenReadings(currentReading.latitude, currentReading.longitude, currentReading.depth, lastBalloonReading.latitude, lastBalloonReading.longitude, lastBalloonReading.depth);
            //Get the magnitude difference
            double magDiff = currentReading.magnitude - lastBalloonReading.magnitude;
            magDiff = fabs(magDiff);
            fprintf(logFile, "Balloon seismic reporting Coord Diff with reporting node %.2fkm\n", distance);

            fprintf(logFile, "Balloon seismic reporting magnitude %.2f\n", lastBalloonReading.magnitude);
            fprintf(logFile, "Balloon seismic reporting magnitude diff with reporting node %.2f\n", magDiff);
            fprintf(logFile, "Communication time: calculate later s\n");
            fprintf(logFile, "Total messages sent between reporting node and base station: 1\n");
            fprintf(logFile, "Number of adjacent matches to reporting node: %d\n", receivedReport.numMatches);
            fprintf(logFile, "Coordinate difference threshold: %fkm\n", DISTANCETHRESHOLD);
            fprintf(logFile, "Magnitude difference threshold: %f\n", MAG_DIFF_THRESHOLD);
            fprintf(logFile, "Earthquake magnitude threshold: %f\n", MAGNITUDETHRESHOLD);
            fprintf(logFile, "===========================================================\n");
        }
        //Continue with while loop

    }
    return 0;

}


void writeAdjacentNodeReadings(SeismicReading nodeReading, SeismicReading reportingNodeReading, int rank){
    //Get the distance
    double distance = getDistanceBetweenReadings(reportingNodeReading.latitude, reportingNodeReading.longitude, reportingNodeReading.depth, nodeReading.latitude, nodeReading.longitude, nodeReading.depth);
    //Get the magnitude difference
    double magDiff = reportingNodeReading.magnitude - nodeReading.magnitude;
    magDiff = fabs(magDiff);
    //Get the coordinates
    int x = getX(rank, ncols);
    int y = getY(rank, nrows);
    fprintf(logFile, "%d (%d,%d)                     (%02.2f, %02.2f)    %.2fkm     %.2f         %.2f\n", rank, x, y, nodeReading.latitude, nodeReading.longitude, distance, nodeReading.magnitude, magDiff);
}

//Pthread function for the balloon
void* balloonFunction(void *pArg){

    printf("Main process balloon thread running\n");
    //Seed the random number generator
    srand((unsigned) time(NULL) * sensor_rank);
    while (termination == 0)
    {
        SeismicReading balloonReading = generateReading(1);
        balloonReadings[readingIndex] = balloonReading;
        readingIndex++;
        readingIndex = readingIndex % BALLOONARRAYSIZE;

        sleep(INTERVAL);
    }
    return NULL;
}

//Pthread function for the send and receive thread of the main function
void* baseSendRecv(void *pArg){
    //Log down the start time
    time_t start_t, end_t,s, val;
    time(&start_t);
    //Total run time is number of iterations times the interval that we are running at for reading generation
    double totalRunTime = iterations * INTERVAL;
    //Flag
    int flag = 0;
    //Probe status
    MPI_Status probeStatus;
    //Receive buffer
    int reportBuffer[6];

    MPI_Status status;

    MPI_Request send_requests;

    //Received reading
    SeismicReading reportingReading;
    printf("Main process send and receive thread running\n");

    //Sentinel file name
    char sentinel[] = "sentinel";

    struct timespec clock_time;

	double elapsed_time;

	int i;

	FILE* ptr;

	char word[100];

	clock_gettime(CLOCK_MONOTONIC, &clock_time);
    struct timespec *node_last_alive_timestamp_array = (struct timespec*)malloc(size * sizeof(struct timespec));
        for(i = 1; i < size; i++){
			node_last_alive_timestamp_array[i].tv_sec =  clock_time.tv_sec;
			node_last_alive_timestamp_array[i].tv_nsec =  clock_time.tv_nsec;
		}




    //While loop that keeps running till termination is 1
    while(termination == 0){
        if(ifFileExists(sentinel)){
            ptr = fopen("sentinel", "r");
            fgets(word,100,ptr);
                if (strcmp(word,"stop\n") == 0){
                    printf("Keyword stop detected from the sentinel file, terminating program now\n");
                    int term = 1;
                    termination = 1;
                    fclose(ptr);
                    fprintf(logFile, "Keyword stop detected from the sentinel file\n");
                    fclose(logFile);
                    //Send data to all
                for(int i = 1; i < size; i++){
                    MPI_Send(&term, 1, MPI_INT,i, MSG_TERM, MPI_COMM_WORLD);
                }
            break;
        }

            fclose(ptr);
        }

        //Check if enough time has passed
        time(&end_t);
        if(difftime(end_t, start_t) >= totalRunTime){
            printf("Enough time has passed, terminating program now\n");
            int term = 1;
            termination = 1;
            //Send data to all
            for(int i = 1; i < size; i++){
                MPI_Send(&term, 1, MPI_INT,i, MSG_TERM, MPI_COMM_WORLD);
            }
            //Break out of the loop
            break;
        }


        for(int i = 1; i < size; i++){
                MPI_Isend(&check, 1, MPI_INT,i, MSG_CHECK, MPI_COMM_WORLD,&send_requests);
            }

        //Check for messages from the sensor nodes (MPI_COMM_WORLD)
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &probeStatus);
        while(flag == 1){
            //While we have a message pending,check what it is
            switch(probeStatus.MPI_TAG){
                case MSG_REPORT: //Terminating message pending
                    //Receive it
                    printf("Main server receiving\n");
                    //Set the alert time
                    //Get current time
                    s = time(NULL);
                    //Uses the current time to create a tm struct
                    reported_Time = localtime(&s);

                    MPI_Recv(&receivedReport, 1, SensorReport, probeStatus.MPI_SOURCE, MSG_REPORT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    printf("Main server received report from %d, U: %d, D: %d, L: %d, R: %d\n", receivedReport.reportingNode, receivedReport.neighbours[0], receivedReport.neighbours[1], receivedReport.neighbours[2],receivedReport.neighbours[3]);
                    //Set the needs processing flag to true so our other thread can process this report
                    needsProcessing = 1;
                    break;
                    case MSG_ALIVE:
                    MPI_Recv(&alive, 1, MPI_INT, probeStatus.MPI_SOURCE, MSG_ALIVE, MPI_COMM_WORLD, &status);
                    clock_gettime(CLOCK_MONOTONIC, &clock_time);
                    node_last_alive_timestamp_array[probeStatus.MPI_SOURCE].tv_sec =  clock_time.tv_sec;
                    node_last_alive_timestamp_array[probeStatus.MPI_SOURCE].tv_nsec =  clock_time.tv_nsec;
                    for(i = 1; i < size; i++){
                        elapsed_time = (clock_time.tv_sec - node_last_alive_timestamp_array[i].tv_sec) * 1e9;
                        elapsed_time = (elapsed_time + (clock_time.tv_nsec - node_last_alive_timestamp_array[i].tv_nsec)) * 1e-9;
                        printf("Process - %d. Last active: %lf seconds ago\n", i, elapsed_time);
                    if (elapsed_time > 1.5) {
                        printf("Node failure detected, exiting...\n");
                        fprintf(logFile, "Node failure\n");
                        fprintf(logFile, "Error code: %d\n", error);
                        fclose(logFile);
                        int term = 1;
                        termination = 1;
                        flag = 0;
                        free(node_last_alive_timestamp_array);
                        for(int i = 1; i < size; i++){
                            MPI_Send(&term, 1, MPI_INT,i, MSG_TERM, MPI_COMM_WORLD);
                        }
                        printf("Gracefully shutting down...\n");
                        return NULL;

                    }


                }
                    break;

            }

            MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &probeStatus);
        }

    }


    return NULL;

}


int seismic_sensor(MPI_Comm world_comm){

    //Coordinates array to store the current process's coordinate
    int coords[2];
    //Get the sensor rank of this process
    MPI_Comm_rank(seismic_Comm, &sensor_rank);
    //Get the current process's coordinates
    MPI_Cart_coords(comm2D, sensor_rank, NDIMS, coords); // coordinated is returned into the coord array
    //Using the coordinates to find the current rank
    MPI_Cart_rank(comm2D, coords, &my_cart_rank);

    //Storing the file anme
    char fileName[32] = {0};
    //Get the file name
    sprintf(fileName, "sensor%d_Debug.txt", my_cart_rank);
    //Open the file
    logFile = fopen(fileName, "w");

    //Spawn the extra thread
    pthread_t tid;
    pthread_create(&tid, 0, sensorFunc, world_comm); // Create the thread, pass in the world communicator

    //This thread is used to generate the random readings
    //While we have not received a termination message
    //Seed the random number generator
    srand((unsigned) time(NULL) * sensor_rank);

    while(termination == 0){
        //Increment our iteration count
        iterationCount += 1;
        //Generate the readings, USE MUTEX LOCK LATER
        lastReading = generateReading(0);
        //Set need neighbouring data to 1 for a new reading
        if(lastReading.magnitude >= MAGNITUDETHRESHOLD){
            needNeighbourData = 1;
            fprintf(logFile,"ALERT: %d/%d/%d, %d:%d:%d, (Lat: %f, Long: %f), Magnitude: %f, Depth: %f \n", lastReading.year, lastReading.month, lastReading.day, lastReading.hour, lastReading.minute, lastReading.second, lastReading.latitude, lastReading.longitude, lastReading.magnitude, lastReading.depth);

        } else{
            //Log this
            fprintf(logFile,"Reading: %d/%d/%d, %d:%d:%d, (Lat: %f, Long: %f), Magnitude: %f, Depth: %f \n", lastReading.year, lastReading.month, lastReading.day, lastReading.hour, lastReading.minute, lastReading.second, lastReading.latitude, lastReading.longitude, lastReading.magnitude, lastReading.depth);
        }

        //Sleep for "interval" seconds
        sleep(INTERVAL);

    }


    pthread_join(tid, NULL);
    fclose(logFile);
    //Free the communicator
    MPI_Comm_free(&seismic_Comm);
    return 0;
}

//Pthread function, used to handle send and receives
void* sensorFunc(void *pArg){
    //Default Buffer
    int defaultBuffer;
    //Keeps track of the status
    MPI_Status status;
    //Create the arrays of statuses
    MPI_Request send_requests[4];

    MPI_Request receive_requests[4];
    MPI_Status receive_status[4];

    MPI_Status probeStatus;
    //Receive request
    MPI_Request receive_request;
    //Send request
    MPI_Request send_request;

    int flag = 1;

    //Upper and lower neighbour readings
    struct SeismicReading upperReading;
    struct SeismicReading lowerReading;
    //Left and right neighbour readings
    struct SeismicReading leftReading;
    struct SeismicReading rightReading;

    //Received reading
    struct SeismicReading receivedReading;


    //Given a cartesian communicator, the dimension and displacement, get the neighbours
    //Assuming the cartShift actually returns the neighbour's ranks inside the communicator
    MPI_Cart_shift( comm2D, SHIFT_ROW, DISP, &upperNeighbour, &lowerNeighbour);
    MPI_Cart_shift( comm2D, SHIFT_COL, DISP, &leftNeighbour, &rightNeighbour);
    //printf("Process %d U: %d, D: %d, R: %d, L: %d\n", my_cart_rank, upperNeighbour, lowerNeighbour, rightNeighbour, leftNeighbour);
    //Keeps track of how many we have received
    int receiveCount = 0;
    //Keep track of the position of our buffer

    //Count the number of neighbours so we know when
    int numberOfNeighbours = 0;
    if (lowerNeighbour >= 0){
        numberOfNeighbours += 1;
    }
    if (rightNeighbour >= 0){
        numberOfNeighbours += 1;
    }
    if (upperNeighbour >= 0){
        numberOfNeighbours += 1;
    }
    if (leftNeighbour >= 0){
        numberOfNeighbours += 1;
    }
    //Count the number of neighbours that had recorded seismic events greater than our threshol
    int neighbourTriggerCount = 0;


    //While we have not received a termination message
    while(termination == 0){
        //printf("Process %d in while\n", my_cart_rank);
        //Let it sleep for a cycle so the other thread can generate a reading
        usleep(SLEEP_MICRO_SEC);
        //If our last reading's magnitude is greater than our threshold and we have not requested data from neighbiouring nodes
        if(lastReading.magnitude >= MAGNITUDETHRESHOLD && needNeighbourData == 1){
            printf("ALERT: Process %d %d/%d/%d, %d:%d:%d, (Lat: %f,Lon: %f), Magnitude: %f, Depth: %f \n", my_cart_rank, lastReading.year, lastReading.month, lastReading.day, lastReading.hour, lastReading.minute, lastReading.second, lastReading.latitude, lastReading.longitude, lastReading.magnitude, lastReading.depth);
            //Set it to false
            needNeighbourData = 0;
            //Set receive count to 0
            receiveCount = 0;
            //Reset the trigger count
            neighbourTriggerCount = 0;
            //Send requests
            //Send a request to its four neighbours
            MPI_Isend(&my_cart_rank, 1, MPI_INT, lowerNeighbour, MSG_REQ, comm2D, &send_requests[0]);
            MPI_Isend(&my_cart_rank, 1, MPI_INT, upperNeighbour, MSG_REQ, comm2D, &send_requests[1]);
            MPI_Isend(&my_cart_rank, 1, MPI_INT, rightNeighbour, MSG_REQ, comm2D, &send_requests[2]);
            MPI_Isend(&my_cart_rank, 1, MPI_INT, leftNeighbour, MSG_REQ, comm2D, &send_requests[3]);
        }

        //Check for messages from the base station (MPI_COMM_WORLD)
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &probeStatus);
        while(flag == 1){
            //While we have a message pending,check what it is
            switch(probeStatus.MPI_TAG){
                case MSG_TERM: //Terminating message pending
                    //Receive it
                    MPI_Recv(&termination, 1, MPI_INT, 0, MSG_TERM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //Sleep for some time
                    usleep(SLEEP_MICRO_SEC);
                break;
                case MSG_CHECK:
                    MPI_Recv(&check, 1, MPI_INT, 0, MSG_CHECK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&alive, 1, MPI_INT, 0, MSG_ALIVE, MPI_COMM_WORLD);
                break;
            }
            //Check for messages from the base station (MPI_COMM_WORLD)
            MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &probeStatus);
        }

        //Check if we have any pending receives, if not just ignore this loop
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, comm2D, &flag, &probeStatus);
        while(flag == 1){
            //printf("Process %d probing\n", my_cart_rank);
            //We have pending messages, call a blocking receive
            //Check the tag to know what we should be receiving
            switch(probeStatus.MPI_TAG){
                //Other node sent a request for data
                case MSG_REQ:
                    //Call MPI_RECV with default buffer
                    MPI_Recv(&defaultBuffer, 10, MPI_INT, probeStatus.MPI_SOURCE, MSG_REQ, comm2D, &status);
                    //printf("Process %d received request from %d \n", my_cart_rank, probeStatus.MPI_SOURCE);
                    //Send the data to them, MAKE A MUTEX LOCK LATER HERE!!!!!!!!!!!!!!!!!!!!!
                    MPI_Send(&lastReading, 1, SeismicData, probeStatus.MPI_SOURCE, MSG_DATA, comm2D);
                    //Case when other nodes are requesting for messages
                    //printf("Process %d sent data to %d \n", my_cart_rank, probeStatus.MPI_SOURCE);
                    fprintf(logFile, "Request for data received from sensor %d, data has been sent\n", probeStatus.MPI_SOURCE);
                    break;
                case MSG_DATA: //We are receiving data
                    //Check the source and store it in the right variable
                    if(probeStatus.MPI_SOURCE == lowerNeighbour){
                        MPI_Recv(&lowerReading, 1, SeismicData, probeStatus.MPI_SOURCE , MSG_DATA, comm2D, &status);
                        //Assign the received reading to point to the reading we just received
                        receivedReading = lowerReading;
                    }
                    else if(probeStatus.MPI_SOURCE == upperNeighbour){
                        MPI_Recv(&upperReading, 1, SeismicData, probeStatus.MPI_SOURCE , MSG_DATA, comm2D, &status);
                        //Assign the received reading to point to the reading we just received
                        receivedReading = upperReading;
                    }
                    else if(probeStatus.MPI_SOURCE == rightNeighbour){
                        MPI_Recv(&rightReading, 1, SeismicData, probeStatus.MPI_SOURCE , MSG_DATA, comm2D, &status);
                        //Assign the received reading to point to the reading we just received
                        receivedReading = rightReading;
                    }
                    else if(probeStatus.MPI_SOURCE == leftNeighbour){
                        MPI_Recv(&leftReading, 1, SeismicData, probeStatus.MPI_SOURCE , MSG_DATA, comm2D, &status);
                        //Assign the received reading to point to the reading we just received
                        receivedReading = leftReading;
                    }
                    //Calculate the distance between our last reading and the received reading
                    //Check if it is above the threshhold, if it is, store it
                    double distance = getDistanceBetweenReadings(lastReading.latitude, lastReading.longitude, lastReading.depth, receivedReading.latitude, receivedReading.longitude, receivedReading.depth);
                    fprintf(logFile, "Received from sensor %d, %d:%d:%d, %d:%d:%d, (Lat: %f, Long: %f), Magnitude: %f, Depth: %f distance = %fkm\n", status.MPI_SOURCE, receivedReading.year, receivedReading.month, receivedReading.day, receivedReading.hour, receivedReading.minute, receivedReading.second, receivedReading.latitude, receivedReading.longitude, receivedReading.magnitude, receivedReading.depth , distance);
                    //Increment the receive count
                    receiveCount += 1;
                    //Check if the distance is within threshold
                    if(distance <= DISTANCETHRESHOLD){
                        //Calculate the absolute difference between our reading and received reading
                        double magDiff = lastReading.magnitude - receivedReading.magnitude;
                        magDiff = fabs(magDiff);
                        //Check if the magnitude is greater than our threshold
                        if(magDiff <= MAG_DIFF_THRESHOLD){
                            //Our neighbour triggers it as well, increment our counter
                            neighbourTriggerCount += 1;
                        }
                    }
                    //Check if we have received all the neighbours
                    if(receiveCount == numberOfNeighbours){
                        //Reset the value
                        receiveCount = 0;
                        //We have, now check if we passed the threshhold for number of neighbours
                        if(neighbourTriggerCount >= 2){ //If we have 2 or more neighbours also having similar readings
                            //Create a Seismic Report
                            SeismicReport report;
                            report.reportingNode = my_cart_rank;
                            report.iteration = iterationCount;
                            report.numMatches = neighbourTriggerCount;
                            report.neighbours[0] = upperNeighbour;
                            report.neighbours[1] = lowerNeighbour;
                            report.neighbours[2] = leftNeighbour;
                            report.neighbours[3] = rightNeighbour;
                            //Flatten our readings into an array
                            flattenToArray(lastReading, report.reportingNodeReading);
                            flattenToArray(upperReading, report.upperNeighbourReading);
                            flattenToArray(lowerReading, report.lowerNeighbourReading);
                            flattenToArray(leftReading, report.leftNeighbourReading);
                            flattenToArray(rightReading, report.rightNeighbourReading);

                            // printf("Process %d before sending upper reading %f %f\n", my_cart_rank,upperReading.magnitude, upperReading.depth);

                            // printf("Process %d before sending upper array %f %f\n", my_cart_rank,report.upperNeighbourReading[8], report.upperNeighbourReading[9]);

                            //Sends the report with MSG_REPORT as tag to the root process
                            MPI_Send(&report, 1, SensorReport, 0, MSG_REPORT, MPI_COMM_WORLD);
                            fprintf(logFile, "Report sent to main\n");
                        }
                    }
                    break;
                    //Other switch case for checking if this node is alive.. implement later
            }
            //Probe again
            MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, comm2D, &flag, &probeStatus);
        }

    }

    printf("Process %d terminated\n", my_cart_rank);
    return NULL;
}

int ifFileExists(const char* filename){
    struct stat buffer;
    return (stat (filename, &buffer) == 0);
}
