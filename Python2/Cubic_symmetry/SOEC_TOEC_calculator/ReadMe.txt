This file calculates SOEC and TOEC from the ab initio quantum energy-strain calculations.
An example input is provided in this folder. Separate numbers in input file by a space or a new line, no commas or semicolons.
#Serves as a start of a comment and the rest of the line is ignored
The output folder contains two files: Constants.txt where the results were calculated from all available data.
LeaveOneOut.txt - where one set of data was omitted each time - Consistency check.

Running the script:

arguments:
Input file - a simple txt file with numbers in the order specified in the example file
Output folder -  a name of a folder into which the results will be stored. The folder must be/will be created into the working directory. 
	If the folder already exists it will store the data into that folder. However if there are any files with the same name, they will be overwritten!!
example:

python Constants.py Example-input Results

This will take the Example-input.txt which is in the inputdata folder and calculates the constants. 
The result will be stored in folder Results that will be created in Inputdata folder.
However if the python script is in the same folder as the input data and thus the working directory the path to the script can be omitted:
python Constants.py Example-input Results



