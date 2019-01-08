TOECalc is a collection of Python scripts for calculation of Second and Third order elastic constants from quantum Energy-Strain calculations and other useful properties of crystals. The implemented equations are mostly based on https://arxiv.org/abs/1111.2737 . Main advantage of these scripts is that they can combine several quantum calculations into a single prediction of elastic constants through a least square method. This greatly increases the accuracy and reduces fluctuation of results as compared to results obtained from single data set, or obtained by averaging the results of individually calculated constants.

The code is designed to be run from the command line. You can either directly specify the required parameters, or if not sure which are needed running the code without any parameter will cause the code to ask you for all relevant parametrs.

Most of the different functions require an input file, and all of them create an output file. These files will be looked for or stored in the working directory. When running from command line, the working directory is the directory you are currently in.

In general the scripts are called in the following manner:

python C:/path_to_TOECalc/.../TOECalc.py argument1 argument2 argument3

In ReadMe for each specific script you will find specific information on the arguments for that code.
If the arguments are not what the code expects, it would ask you to input them manually.

Additional informations are in folders with each script. The code is also heavily annotated so a help() command might give you some other insight. If still unsure, feel free to contact me on Github. I am happy to help.

author: Pavel Kroupa
TOECalc was created during an internship at Czech Academy of Science, department of Materials, under the supervision of Mgr. Martin Friák PhD.

The author takes no responsibility for the results produced by this code. The calculations were tested and verified and the code was used in published work, but it was not possible to test for all possible edge cases.
