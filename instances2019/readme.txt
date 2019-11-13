

This document explains the instance file for the truck and freighter routing problem. 

Three fines are necessary to run an instance:
- the vertex list
- the parameter file
- the distanceMatrix

all files are made from the instanceNantes.txt file from which you can pick locations to create larger instances if you need

==========================
VERTEX LIST:

example
CUSTNO	DEMAND	VERTEX_TYPE	LATITUDE	LONGITUDE	TW_OPENING
1	0	D	47.282609	-1.520303	0
15	0	P	47.234018	-1.535018	0
...
3	40	S	47.214248	-1.557261	10800
...
2	20	LS	47.256911	-1.513408	23400
...

FIELDS:
1 - CUSTNO index of the vertex (also the index in the time matrix)
2 - DEMAND	quantity ordered by customers, 
3 - VERTEX_TYPE D for depot, P for parkings, S for customers delivered only by the small truck, L for customers of the big truck and LS if the customer can be served by both types. 
4 - LATITUDE
5 - LONGITUDE
6 - TW_OPENING Time at which the service can start at the vertex. 

The depot and the parkings have a time window width equal to the time horizon, customers have the width that is specified in the header

==========================
PARAMETERS: 
A single parameter file is given: parameters.txt

HEADER :

example
# SMALL VEHICLE CAPACITY	Speed ratio small vs large	Time horizon	TW width	service duration
400	2	32400	7200	300

FIELDS: 
1 - The first number in the header is the capacity of the small vehicle
2 - Speed ratio small vs large: a ratio of 2 means that the time needed by the small vehicle to travel between two points is twice the time indicated in the time matrix
3 - Time horizon: time at which all vehicles should be back at the depot - they are supposed to leave the depot at time 0.
4 - TW width: with for all customer time windows in the problem
5 -	duration of the service at each customer (same for all customers and also for coupling/decoupling/transfers at parking spaces). 

You can create separate files to run instances with various hypothesis

DISTANCE MATRIX

this file contains the large vehicle travelling times in seconds between all pair of vertices of the instanceNantes.txt file. 





