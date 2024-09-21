The algorithm takes a CSV dataset as input with some changes in order for the csv parser to correctly parse.
1)on the first line of the csv you need to have the value types of the attributes that will follow on next lines.
2)all the categorical attributes should be encoded to integers
3)the class of the datapoint should be the last attribute (last column of dataset) 
e.g. we have a dataset for dogs, at the first column there is the gender of the dog(M,F) ,at the second one the weight(in grams ),at the third the age,and the target class race at the forth column
the data set would look like this at start : 

M,13455,2,Labrador
.....
F,18654,4,GermanShepherd

it should be modified like the following in order to be correctly parsed :

int,float,int,int
1,13455,2,1
...
2,18654,4,7

at the above example gender is encoded as m->1,f->2 and the target class as Labrador->1 GermanShepherd->7

THE ALGORITHM MISSES THE BAGGING AND BOOSTING FOR THE DATASET SO IT DOES NOT ACHIEVE THE EXPECTED ACCURACY.
