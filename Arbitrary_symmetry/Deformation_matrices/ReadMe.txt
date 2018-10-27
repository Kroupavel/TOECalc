This program create the necessary deformation matrices that is needed for ab initio calculations. 
The format is adjusted to serve as input for VASP.
12 matrices are generated based on six matrices A1-A6 and a positive and negative delta.
The input of this script is a value of delta and a name of the output file without .txt .
calling the script:
python Matrices.py 0.015 Matrices1

If the python script lies in the working directory, the whole path to the script can be omitted
