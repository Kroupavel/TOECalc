This program create the necessary deformation matrices that is needed for ab initio calculations. 
The format is adjusted to serve as input for VASP.
12 matrices are generated based on six matrices A1-A6 and a positive and negative delta.
The input of this script is a value of delta and a name of the output file without .txt .
The script should be called from the folder which contains the training data
calling the script:

python C:/path_to_Matrices/Matrices.py 0.015 Matrices1
