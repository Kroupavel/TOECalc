Dear User,
thank you for using TOECalc. If you have any problems running the code, feel free to contact me at kroupavel(at)gmail(dot)com. I am happy to help.
If you are using this script for publishing, please acknowledge the use of this software by citing the article TBC.

To run this code, you need to have a working installation of Python (version 2.7.6). 
A simple solution to this can be to download an interpreter like Canopy (which was used in making of TOECalc).
The code then can be run through the interpreter or in a command line(this option is recommended).

To run it in the interpreter change the working directory to the folder that contains the input files with the numerical data you would like to process.
Then run the code, it will ask for some input arguments, depending on the script run. You can either add them one by one, or faster method is possible through the command line
When all the arguments are added, it might take a while to do all the calculations but after a while all output files should appear in the output folder.

When running through the command line, few problems might occur at first. Typical ones are missing library or cmd not finding any installation of python. 
There are plenty of guides on the internet that help you solve that problem. If you keep running into problem and no solution on internet works, contact me.
To actually run any of those scripts: If an input file is needed, change the working directory (writing: cd C:complete\path\) to that of the input file.
If no input files are needed the working directory is where the results will be stored. 
Then a following command will run the script:

python C:complete\path\to\the\script\TOECalc.py argument1 argument2 argument3

Be sure to include the name of the script to run with the .py extension. If you don't know which arguments to include, don't put any and you will be ask for them.
All words separate with space, no commas or dots. For decimal numbers use dot.
In ReadMe for each specific script you will find specific information on the arguments for that code.
If the arguments are not what the code expects, it would ask you to input them manually.

Scripts:
Deformation matrices generator - when given a strain value will calculate the necessary deformation matrices to use in the quantum calculations.
TOECalc-complete - This is the main script that has all the functions of other scripts. From the ab initio results, it will calculate SOEC, TOEC and other properties. 
SOEC,TOEC calculator - This will calculate the SOEC and TOEC from ab initio data only.
Postprocessing - This takes the input of SOEC and TOEC and do the post processing (Young moduli, elasticity...) for this given data. Good for analysis of experimental measurements of SOEC and TOEC.
GUI Constant calculator + matrix generator - This will create a simple GUI in which the SOEC and TOEC can be calculated. Also the deformation matrices can be generated in the GUI.
			- The main advantage is that the data will be sorted and you can manually pick which sets you want to use 
			- e.g. compare the result if lower 5 strain values are used to the results gained when the higher 5 strain values are used.

Additional informations are in the folders with the script. The code is also heavily annotated so a help() command might give you some other insight.

author: Pavel Kroupa
TOECalc was created during an internship at Czech Academy of Science, department of Materials, under the supervision of Mgr. Martin Friák PhD.
