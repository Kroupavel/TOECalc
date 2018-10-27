This program is responsible for the post-processing of the second and third order elastic constants (SOEC, TOEC).
The input of the constants must be made through a file (see example) - The file must contain all the SOEC and TOEC in specific order, separated with a space or a new line.
The input arguments are:  Input file, Output file, Pressure change, number of iteration steps, number of integration steps.
The input file must be in the working directory. The name must be inputted without the .txt suffix.
The output file will create a separate folder in the working directory of that name.
	If a folder of that name already exists it will use that folder and might overwrite files in it!
The hydrostatic pressure applied in GPa. So to simulate 1 GPa, write 1.
The number of iteration steps will determine the density of points for "surface integral" when calculating the directional Young moduli. 
        It also sets the smoothness of the final plots.
        This number affects the computing time the most so it is important to set it accordingly to desired result!:
        For normal plot 200-300 is a recommended value.
        For publication level graphs a 400-500 might be useful at a cost of higher computing time.
        If only numerical results are needed it can be set to 1 which will greatly decrease the computing time.
The number of integration steps is for calculating the pressure derivatives. The effect on computing time is small.
        For normal pressure changes like 1,10,100 GPa the effect of n. of iteration steps is also small. However a 1000 is used as a standard value.   
        
The output are three graphs (Young moduli (E) at zero pressure, the difference in E for given difference in pressure and the same graph normalized, all as a function of direction)
Also a text file is created with numerical data (more details in store method). 

cd C:\Users\Pavel\Documents\Python\data\constants
python C:\Users\Pavel\Documents\Python\TOECalc\Post_processing\Postprocessing.py Cuconstants Curesults 1 200 1000

if the python script is in the same folder as the input file (the working directory) the whole path can be omitted when calling the script, e.g:
python Postprocessing.py Cuconstants Curesults 1 200 1000