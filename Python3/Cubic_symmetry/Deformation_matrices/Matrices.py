import numpy as np
from cmd import Cmd
import sys

class Main(Cmd):
    """
    This program create the necessary deformation matrices that is needed for ab initio calculations. 
    The format is adjusted to serve as input for VASP.
    12 matrices are generated based on six matrices A1-A6 and a positive and negative delta.
    The input of this script is a value of delta and a name of the output file without .txt .
    e.g: a=Main(0.021, "matrices")
    """
    def __init__(self):
        if (len(sys.argv)!= 3): 
           
            self.delta = float(input("Strain:"))
            self.filename = input("Output file:") 

        else:
            self.delta = float(sys.argv[1])
            self.filename = str(sys.argv[2])
            
        self.f = open(self.filename+".txt", 'w')
        self.calculate()
        
    def calculate(self):
        """
        This method contains the formulas of the matrices which are added to a Identity matix,
        to produce the defromation matrices for a given delta."""
        self.base = np.array([[1.0, 0., 0.],
                                [0., 1.0, 0.],
                                [0., 0., 1.0]])
                                
        self.A1 = np.array([[self.delta, 0., 0.],
                                [0., 0., 0.],
                                [0., 0., 0.]])
                                
        self.A2 = np.array([[self.delta, 0., 0.],
                                [0., self.delta, 0.],
                                [0., 0., 0.]])
                                
        self.A3 = np.array([[self.delta, 0., 0.],
                                [0., self.delta, 0.],
                                [0., 0., self.delta]])
                                
        self.A4 = np.array([[self.delta, 0., 0.],
                                [0., 0., self.delta],
                                [0., self.delta, 0.]])
                                
        self.A5 = np.array([[self.delta, self.delta, 0.],
                                [self.delta, 0., 0.],
                                [0., 0., 0.]])
            
        self.A6 = np.array([[0., self.delta, self.delta],
                                [self.delta, 0., self.delta],
                                [self.delta, self.delta, 0.]])    
                             
        i=1
        self.list = (self.A1, self.A2, self.A3, self.A4, self.A5, self.A6)
        for matrix in self.list:
            self.matrix=self.base+matrix
            self.negmatrix = self.base-matrix
            self.f.write("#A%s for delta =%s: \n" %(i, self.delta))
            self.f.write('\n'.join('  '.join(str(cell) for cell in row) for row in self.matrix))   
            self.f.write(" \n#A%s for delta = -%s: \n" %(i, self.delta))
            self.f.write('\n'.join('  '.join(str(cell) for cell in row) for row in self.negmatrix))
            self.f.write(" \n\n")
            i+=1

        self.f.close() 
        
a=Main()