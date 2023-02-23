//This file contains functions to call for creating MPI Datatypes from Structs
#include <mpi.h>
#include "sensor.c"

//This function takes in a pointer to the datatype variable and creates a MPI datatype stored at that address
void createSeismicReadingDatatype(MPI_Datatype *datatype){
    //Keeps track of the received seismic reading
    struct SeismicReading receivedReading;
    //Create the data type
    //Specify the data types inside the struct
    MPI_Datatype type[10] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    //Specifies how many of each type there are in each block (a block represents a type)
    int blocklen[10] = {1,1,1,1,1,1,1,1,1,1};
    //Array used to represent addresses of each block
    MPI_Aint disp[10];

    //Get the base address
    MPI_Aint base_address;
    MPI_Get_address(&receivedReading, &base_address);


    //Get the addresses and store it inside our array
    MPI_Get_address(&receivedReading.year, &disp[0]);
    MPI_Get_address(&receivedReading.month, &disp[1]);
    MPI_Get_address(&receivedReading.day, &disp[2]);
    MPI_Get_address(&receivedReading.hour, &disp[3]);
    MPI_Get_address(&receivedReading.minute, &disp[4]);
    MPI_Get_address(&receivedReading.second, &disp[5]);
    MPI_Get_address(&receivedReading.latitude, &disp[6]);
    MPI_Get_address(&receivedReading.longitude, &disp[7]);
    MPI_Get_address(&receivedReading.magnitude, &disp[8]);
    MPI_Get_address(&receivedReading.depth, &disp[9]);
    //Make relative, so the address of the second variable is the address of the second - address of the first.
    disp[0] = MPI_Aint_diff(disp[0], base_address);
    disp[1] = MPI_Aint_diff(disp[1], base_address);
    disp[2] = MPI_Aint_diff(disp[2], base_address);
    disp[3] = MPI_Aint_diff(disp[3], base_address);
    disp[4] = MPI_Aint_diff(disp[4], base_address);
    disp[5] = MPI_Aint_diff(disp[5], base_address);
    disp[6] = MPI_Aint_diff(disp[6], base_address);
    disp[7] = MPI_Aint_diff(disp[7], base_address);
    disp[8] = MPI_Aint_diff(disp[8], base_address);
    disp[9] = MPI_Aint_diff(disp[9], base_address);
    // Create MPI struct, 10 represents how many blocks there are (how many variables there are inside our struct), blocklen represents how many elements are in each block
    // e.g if its just an int its 1, if its an array of ints it would be the length of that. The type signals the value types. We then store that inside our valuetype
    // variable
    MPI_Type_create_struct(10, blocklen, disp, type, datatype);
    //Commits the data type
    MPI_Type_commit(datatype);
}


//This function takes in a pointer to the datatype variable and creates a MPI datatype stored at that address
void createSeismicReportStruct(MPI_Datatype *datatype, MPI_Datatype SeismicData){
    //Keeps track of the received seismic reading
    struct SeismicReport report;
    //Create the data type
    //Specify the data types inside the struct
    MPI_Datatype type[9] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT,MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    //Specifies how many of each type there are in each block (a block represents a type)
    int blocklen[9] = {1,1,4,1,10,10,10,10,10};
    //Array used to represent addresses of each block
    MPI_Aint disp[9];
    //Get the base address
    MPI_Aint base_address;
    MPI_Get_address(&report, &base_address);
    //Get the addresses and store it inside our array
    MPI_Get_address(&report.reportingNode, &disp[0]);
    MPI_Get_address(&report.iteration, &disp[1]);
    MPI_Get_address(&report.neighbours, &disp[2]);
    MPI_Get_address(&report.numMatches, &disp[3]);
    MPI_Get_address(&report.reportingNodeReading, &disp[4]);
    MPI_Get_address(&report.upperNeighbourReading, &disp[5]);
    MPI_Get_address(&report.lowerNeighbourReading, &disp[6]);
    MPI_Get_address(&report.leftNeighbourReading, &disp[7]);
    MPI_Get_address(&report.rightNeighbourReading, &disp[8]);
    //Make relative, so the address of the second variable is the address of the second - address of the first.
    disp[0] = MPI_Aint_diff(disp[0], base_address);
    disp[1] = MPI_Aint_diff(disp[1], base_address);
    disp[2] = MPI_Aint_diff(disp[2], base_address);
    disp[3] = MPI_Aint_diff(disp[3], base_address);
    disp[4] = MPI_Aint_diff(disp[4], base_address);
    disp[5] = MPI_Aint_diff(disp[5], base_address);
    disp[6] = MPI_Aint_diff(disp[6], base_address);
    disp[7] = MPI_Aint_diff(disp[7], base_address);
    disp[8] = MPI_Aint_diff(disp[8], base_address);
    // Create MPI struct, 10 represents how many blocks there are (how many variables there are inside our struct), blocklen represents how many elements are in each block
    // e.g if its just an int its 1, if its an array of ints it would be the length of that. The type signals the value types. We then store that inside our valuetype
    // variable
    MPI_Type_create_struct(9, blocklen, disp, type, datatype);
    //Commits the data type
    MPI_Type_commit(datatype);
}