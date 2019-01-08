TOECalc is a collection of Python scripts for calculation of Second and Third order elastic constants from quantum Energy-Strain calculations.

If you are using this script for publishing, please acknowledge the use of this software by citing the article TBC.

To run this code, you need to have a working installation of Python (version 2.7.6). 
The code then can be run through the interpreter or in a command line(the second option is recommended).

To run it in the interpreter change the working directory to the folder that contains the input files with the numerical data you would like to process.
Then run the code, it will ask for some input arguments, depending on the script you are running. You can either add them one by one, when asked for them, or faster method is to add them prior to running as parameters.
When all the arguments are added, it might take a while to do all the calculations but after a while all output files should appear in the output folder.

To actually run any of those scripts: If an input file is needed, change the working directory (writing: cd C:complete\path\) to that of the input file.
If no input files are needed the working directory is where the results will be stored. 
Then a following command will run the script:

python TOECalc.py argument1 argument2 argument3

Be sure to include the name of the script to run with the .py extension. If you don't know which arguments to include, don't put any and you will be ask for them.
All words separate with space, no commas or dots. For decimal numbers use dot.
In ReadMe for each specific script you will find specific information on the arguments for that code.
If the arguments are not what the code expects, it would ask you to input them manually.

There are two main folders: 
Cubic_symmetry - Contains scripts for calculating Elastic constants for cubic crystals only. Additionally contains further processing of these results. Mostly based on https://arxiv.org/abs/1111.2737
Arbitrary_symmetry - Contains scripts for calculating elastic constants for crystals of arbitrary symmetry. Does not contain any further processing as such calculations are typically structure dependent. 



Additional informations are in folders with each script. The code is also heavily annotated so a help() command might give you some other insight.

author: Pavel Kroupa
TOECalc was created during an internship at Czech Academy of Science, department of Materials, under the supervision of Mgr. Martin Friák PhD.
