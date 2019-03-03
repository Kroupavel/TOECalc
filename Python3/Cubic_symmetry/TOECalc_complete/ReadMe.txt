The is the complete version of TOECalc. It calculates second and third order elastic constants (SOEC, TOEC) from 
quantum ab initio energy-strain calculations and further process the results into other useful properties of material, such as Young moduli and pressure derivatives. 
         
The output with the results will be stored in a new folder that will be created in the working directory. 
        
The following input arguments are needed:
Name of the input file, without the .txt extension. An example of input file is  provided. 
	The input file must be in the working directory. 
Name of the folder into which to store the results. The folder will be created into the working directory. 
	If a folder of that name already exists it will store the data into that folder, but any content in the folder with the same name will be overwritten.
Pressure change(in GPa): The code calculates the change of SOEC and other properties when a hydrostatic pressure is applied. 
	The magnitude of the pressure is specified by the user through this parameter. To input 1 GPa write 1
Number of iteration steps will determine the density of points for "surface integral" when calculating the directional Young moduli. 
	It also sets the smoothness of the final plots.
        This number affects the computing time the most so it is important to set it accordingly to desired result!:
        For normal plot 200-300 is a recommended value.
        For publication level graphs a 400-500 might be useful at a cost of higher computing time(up to a minute).
        If only numerical results are needed it can be set to 1 which will greatly decrease the computing time.
The number of integration steps is for calculating the pressure derivatives. The effect on computing time is small.
        For normal pressure changes like 1,10,100 GPa the effect on the accuracy is small. 1000 is used as a standard value. 

for example: python C:/path_to_TOECalc/TOECalc.py example Output 1 200 1000

assuming the working directory is changed to where the input file is eg: cd C:\...\TOECalc\data
If the python script is not in the same folder as the input file, then the whole path to the python script needs to be specified as shown above.
If however the script is in the same folder, only the name of the script can be written, e.g: python TOECalc.py example Output 1 200 1000

The output are three files and three graphs: Constants.txt contains the SOEC and TOEC calculated from all the available data.
LeaveOneOut.txt stores the results of calculating the constants while omitting a different set of data each time. This serves as a consistency check for outliers and invalid data.
Additional data.txt contains other properties of material such as: anisotropy, directional Young moduli, SOEC, all at zero and user-defined pressure. It also contains some polycrystalline properties like shear.
Three graphs are produced: 1) Young moduli (E) at zero pressure, 2) the difference in E for given difference in pressure and 3) the same graph normalized by the Young moduli at zero pressure, all as a function of direction
