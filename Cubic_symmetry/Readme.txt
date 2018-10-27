This folder contains the following:

Scripts:
Deformation matrices generator - when given a strain value will calculate the necessary deformation matrices to use in the quantum calculations.
TOECalc-complete - This is the main script that has all the functions of other scripts. From the ab initio results, it will calculate SOEC, TOEC and other properties. 
SOEC,TOEC calculator - This will calculate the SOEC and TOEC from ab initio data only.
Postprocessing - This takes the input of SOEC and TOEC and do the post processing (Young moduli, elasticity...) for this given data. Good for analysis of experimental measurements of SOEC and TOEC.
GUI Constant calculator + matrix generator - This will create a simple GUI in which the SOEC and TOEC can be calculated. Also the deformation matrices can be generated in the GUI.
			- The main advantage is that the data will be sorted and you can manually pick which sets you want to use 
			- e.g. compare the result if lower 5 strain values are used to the results gained when the higher 5 strain values are used.
