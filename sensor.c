//This file contains all helper functions and data types that a seismic sensor needs
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define pi 3.14159265358979323846

//Structure for the randomly generated values
typedef struct SeismicReading{
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second;
    double latitude;
    double longitude;
    double magnitude;
    double depth;
} SeismicReading;

//Structure for reports that will be sent to the base station
typedef struct SeismicReport{
    int reportingNode;
    int iteration;
    int neighbours[4];
    int numMatches;
    double reportingNodeReading[10];
    double upperNeighbourReading[10];
    double lowerNeighbourReading[10];
    double leftNeighbourReading[10];
    double rightNeighbourReading[10];
} SeismicReport;


double randomNumber(double min, double max){
    //Our range
    double range = max - min;
    //Work out the ratio that we need to divide the rand() value by
    double ratio = RAND_MAX / range;
    //Get our random value
    double randomValue = min + (rand() / ratio);
    return randomValue;
}

//This function takes in a seismic reading and a pointer to an array and flattens the reading into an array form
void flattenToArray(SeismicReading reading, double *array){
    array[0] = reading.year;
    array[1] = reading.month;
    array[2] = reading.day;
    array[3] = reading.hour;
    array[4] = reading.minute;
    array[5] = reading.second;
    array[6] = reading.latitude;
    array[7] = reading.longitude;
    array[8] = reading.magnitude;
    array[9] = reading.depth;
}

//This function takes in an array and converts it back into a seismic reading
void arrayToReading(SeismicReading *reading, double *array){
    reading->year = (int)array[0];
    reading->month = (int)array[1];
    reading->day = (int)array[2];
    reading->hour = (int)array[3];
    reading->minute = (int)array[4];
    reading->second = (int)array[5];
    reading->latitude = array[6];
    reading->longitude = array[7];
    reading->magnitude = array[8];
    reading->depth = array[9];

}


//This function generates a seismic reading and returns it
struct SeismicReading generateReading(int sensorType){
    //Instantiate the structure
    struct SeismicReading reading;


    //Time variables
    time_t s, val;
    struct tm* current_time;
    //Get current time
    s = time(NULL);

    //Uses the current time to create a tm struct
    current_time = localtime(&s);

    //Set the reading's time values
    reading.year = current_time->tm_year + 1900;
    reading.month = current_time->tm_mon + 1;
    reading.day = current_time->tm_mday;
    reading.hour = current_time->tm_hour;
    reading.minute = current_time->tm_min;
    reading.second = current_time->tm_sec;


    //Set the reading's seismic data
    //Generate a value between 0 and negative 33
    reading.latitude = -randomNumber(0.0, 33.0);
    //Generate a value between 151 and 184
    reading.longitude = randomNumber(151.0, 180.0);
    //Generate a magnitude between 0 - 12
    if (sensorType == 0)
        reading.magnitude = randomNumber(0.0, 12.0);
    else
        reading.magnitude = randomNumber(2.5, 12.0);
    //Generate a depth between 0 - 10
    reading.depth = randomNumber(0.0, 10.0);
    //Return the reading
    return reading;
}


//The following functions are based off or copied from of https://www.geodatasource.com/developers/c

/*::  This function converts decimal degrees to radians             :*/
double deg2rad(double deg) {
    return (deg * pi / 180);
}

/*::  This function converts radians to decimal degrees             :*/
double rad2deg(double rad) {
    return (rad * 180 / pi);
}


double getDistanceBetweenReadings(double lat1, double lon1, double depth1, double lat2, double lon2, double depth2){
    double theta, dist, depthDiff;

    theta = lon1 - lon2;
    dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
    dist = acos(dist);
    dist = rad2deg(dist);
    dist = dist * 60 * 1.1515;
    //Get the distance in KM
    dist = dist * 1.609344;
    //Calculate the distance with depth as the third factor using Pythagoras theorem
    //Get the depth difference (absolute value)

    depthDiff = depth1 - depth1;
    depthDiff = fabs(depthDiff);

    //c^2 = a^2 + b^2, we are using that equation to get the distance by using our 2D distance and the depth difference
    dist = sqrt((dist * dist) + (depthDiff * depthDiff));

    return (dist);
}

//This function takes in a rank and returns the Y Coordinate inside the cartesian grid
int getY(int rank, int dimension){
  int y = rank / dimension;
  return y;
} 

//This function takes in a rank and returns the X Coordinate inside the cartesian grid
int getX(int rank, int dimension){
  int x = rank % dimension;
  return x;
}
