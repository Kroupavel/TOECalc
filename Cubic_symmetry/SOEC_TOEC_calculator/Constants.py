from cmd import Cmd
import numpy as np
import sys
import os.path


class Main(Cmd):
    def __init__(self):    
        """
        This program is responsible for afterprocessing of the second and third order elastic constants (SOEC, TOEC).
        The input of the constants must be made through a file (see readme or example for more information)
        The input arguments are:  Input file, Output file
        The input file must be in the working directory. A new file with the Output file name will be created in the working directory.
        The results (two text files) will be stored in the new file.   
        If outputfilename is already used, the file will be overwritten.           
        """              
        if (len(sys.argv)!= 3): 
            self.inp = raw_input("Input file:")
            self.outfile = raw_input("Output file(Select or create a folder to which data is stored):") 

        else:
            self.inp = str(sys.argv[1])
            self.outfile = str(sys.argv[2])

            
        if self.outfile!="":
            if not os.path.exists(self.outfile):
                os.makedirs(self.outfile)
        
        a = Calculate()
        a.initialize(self.inp, self.outfile)
        
class Calculate(Cmd):
    """This program needs to be called with arguments:inputfilename, outputfilename\
    where the names of the files should be without the .txt extension. \
    This program has the following methods:\
    `initialize` prepares the necessary variables and calls the methods
    `calc`  calculates the SOEC and TOEC using all the data available. Producing one set of results.\
    `leaveone` which prints a series of SOEC and TOEC. For each set of results one set of data is omitted.\
    This is a good check of consistency of results and thus the reliability of data and script."""


    def initialize(self,  inp, outpath):

        self.input = inp
        self.outpath = outpath
        self.data =  np.genfromtxt(self.input+".txt", dtype=float, comments="#")
        self.calc()
        self.ndeltas = self.data[0]
        if self.ndeltas>1:
            self.leaveone() 
                         

    def calc(self):
        """calculates the elastic constants from the input values of energy.\
        Uses the Script class which uses the least square method to find a solution to set of equations with 9 unknowns.\
        It creates a file with the constants and the residual, which is a good measure of reliability"""
        self.ndeltas = self.data[0]                         
        self.volume = self.data[1]
        self.eqenergy = self.data[2]
        n=int(self.data[0])
        # Rework the value of strain into needed form.
        i = 0        
        self.delt = [0.0]*n
        for i in range (0, n):                 
            self.delt[i] = self.data[i+3]
                 
        # Energies are prepared for the input format.    
        self.values = [0.0]*((n)*12) #the list to which energies will be stored         
        # inports all the energy values and save them into self.values:      
        for j in range(0, len(self.data) -3-n ):          
            self.values[j] = self.data[j +3+n]
                
        a = Constants(self.ndeltas, self.volume, self.eqenergy,self.delt, self.values )        
        self.results, self.residual = a.solve()
        
        #Storing algorithm
        self.completeName = os.path.join(self.outpath, "Constants.txt")
        f = open(self.completeName, 'w')     
        constants = ["  #C11 in GPa","  #C12 in GPa", "  #C44 in GPa", "  #C111 in GPa", "  #C112 in GPa", "  #C123 in GPa", "  #C144 in GPa","  #C166 in GPa", "  #C456 in GPa"]
        for i in range (0,9):
            f.write(str(self.results[i]))
            f.write(str(constants[i]))
            f.write("\n")
            
            
            
        f.write(str(self.residual[0]))
        f.write(" # residuals")
        f.close()
                
    def leaveone(self):
        """Calculates the constants using Leave one out method. It omitts one set of data out of the calculation.\
        It gradually omits all sets, one at a time and print the results into specified file\
         the function should be called: leaveone inputname outputname   without .txt 
         This serves as a check of reliablity of sets of data. If there is a large deviation, the data might be suspicious of error"""
          
        
        self.completeName2 = os.path.join(self.outpath, "LeaveOneOut.txt")  
        f2 = open(self.completeName2, 'w') # creates an output file
        # Prepares all the data for input into solving script
        self.n = self.data[0]
        self.volume = self.data[1]
        self.eqenergy = self.data[2]
        self.deltas = self.data[3:self.n+3]
        self.energy = self.data[self.n+3:]
        self.values = np.array([0.0]*(int(self.n-1)*12))
        self.delt = np.array([0.0]*(int(self.n-1)))   
        
        # THe main "i" loop is responsible for omitting one set of data each time.
        #This is done through integers being compared to the intiger on the values being sorted  
        # It also stores the values each time they are calculated.
        #The "k" loop is for sorting value of strains 
        #The "j" loop is for sorting the energy values
        for i in range(0, int(self.n)):
            j=0
            help2=0
            help=0
            # These help intigers are used to get rid of a line full of zeros.
            #E.g. after omitting one set of data the next is safed into the place of the previous one, instead of leaving a blank space
  
            for k in range(0, int(self.n)):
                if k!= i:
                    self.delt[k-help] = self.deltas[k]
                else:
                    help=1 
            for j in range (0,  int(self.n*12)):
                if int(j/12)!= i:
                    self.values[j-help2] = self.energy[j]
                else:
                    help2=12
               
            a = Constants((self.n-1), self.volume, self.eqenergy,self.delt, self.values )
            self.results, residual = a.solve()
               
            f2.write("when the set number %d was omitted the following data was obtained\n" %i)
            constants = ["  #C11 in GPa","  #C12 in GPa", "  #C44 in GPa", "  #C111 in GPa", "  #C112 in GPa", "  #C123 in GPa", "  #C144 in GPa","  #C166 in GPa", "  #C456 in GPa"]
            for i in range (0,9):
                f2.write(str(self.results[i]))
                f2.write(str(constants[i]))
                f2.write("\n")
            f2.write(str(residual))
            f2.write(" # residuals\n")
            
            f2.write("\n")
        f2.close() 
  
class Constants:
    """This code calculates the SOECs and TOECs when given input of the form:"
    "Number of measurements, Volume, Eq. Energy, all strain values in a list and all the energies in another list"
    "The energies are in the order A1+, A1-, A2+, A2- for first strain and than repeated for the other values of strain"
    "The solution is found using the least square method numpy.linalg.lstsq"""
          
    def __init__(self, ndeltas = 1, volume = 47.11, eqenergy = 0, delta = [0.021], energy = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])):
        # the energies are stored with equilibirum energy as an independent variable,
       
      self.delta = delta
      self.volume = volume
      self.eqenergy = eqenergy
      self.energy = energy
      self.ndeltas = ndeltas
    
    def solve(self):
        """
        This is the method responsible for the calculation of SOEC and TOEC. When feeded appropriate infromation by its __init__ it
        creates arrays and matrix necessary for least square method solution.
        Mathematically: A*X=B Where X is an 1-D array of the 9 elastic constants
        A is a n*9 matrix where n is equal to 12*number of delta values used. So it is always overdefined
        B is a 1-D array of energy density difference.
        The values in the X array are being calculated for given A and B
        """
        i=0
        a=[] # coefficient matrix
        b = [] # right side matrix (inhomogenous part)
        
        for i in range (0,int(self.ndeltas)):
            self.ddelta = self.delta[i]
            Constants.matrix(self)
            deltasquared = self.delta[i]**2
            deltacubed = abs(self.delta[i]**3)
            
            x=0
            for x in range (0,12):
                count= x+(i*12)  
                             
                Pos = ((self.energy[count]-self.eqenergy)*160.217662000/(self.volume*self.determinant[x])) # changing the energy to joules                   
                b.append(Pos)
                
                    
        # matrix A with coefficients for matrix x with: C11, C111, C12, C112, C123, C44, C144, C166, C456        
            
            a.extend([[0.5*deltasquared, deltacubed/6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.5*deltasquared, -1*deltacubed/6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [1*deltasquared, deltacubed/3, 1*deltasquared, 1*deltacubed, 0.0, 0.0, 0.0, 0.0, 0.0],
                [deltasquared, -1*deltacubed/3, deltasquared, -1*deltacubed, 0.0, 0.0, 0.0, 0.0, 0.0],
                [deltasquared*3/2, deltacubed/2, 3*deltasquared, 3*deltacubed, deltacubed, 0.0, 0.0, 0.0, 0.0],
                [deltasquared*3/2, -deltacubed/2, 3*deltasquared, -3*deltacubed, -deltacubed, 0.0, 0.0, 0.0, 0.0],
                [deltasquared/2, deltacubed/6, 0.0, 0.0, 0.0, 2*deltasquared, 2*deltacubed, 0.0, 0.0],
                [deltasquared/2, -deltacubed/6, 0.0, 0.0, 0.0, 2*deltasquared, -2*deltacubed, 0.0, 0.0],
                [deltasquared/2, deltacubed/6, 0.0, 0.0, 0.0, 2*deltasquared, 0.0, 2*deltacubed, 0.0],
                [deltasquared/2, -deltacubed/6, 0.0, 0.0, 0.0, 2*deltasquared, 0.0, -2*deltacubed, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 6*deltasquared, 0.0, 0.0, 8*deltacubed],
                [0.0, 0.0, 0.0, 0.0, 0.0, 6*deltasquared, 0.0, 0.0, -8*deltacubed]])
                     
        newb = np.array(b)    
        amatrix = np.matrix(a)
        x, residues, rank, s = np.linalg.lstsq(amatrix, newb)

        #rearrange the output. This is as the ideal output is in different order than the equations
        dummy = x[1]
        x[1] = x[2]
        x[2] = x[5]
        x[5] = x[4]
        x[4] = x[3]
        x[3] = dummy

        return x, residues
        
        
    def matrix(self):
        """
        This method generates the deformation matirces.
        The determinant is calculated and used as a correction to the volume of the crystal.
        (As the volume of the crystal change from the equilibirium value when deformed and it is important to incorporate this change in the script)"""
        self.mdelta = self.ddelta
        self.base = np.array([[1.0, 0., 0.],
                                [0., 1.0, 0.],
                                [0., 0., 1.0]])
                                
        self.A1 = np.array([[self.mdelta, 0., 0.],
                                [0., 0., 0.],
                                [0., 0., 0.]])
                                
        self.A2 = np.array([[self.mdelta, 0., 0.],
                                [0., self.mdelta, 0.],
                                [0., 0., 0.]])
                                
        self.A3 = np.array([[self.mdelta, 0., 0.],
                                [0., self.mdelta, 0.],
                                [0., 0., self.mdelta]])
                                
        self.A4 = np.array([[self.mdelta, 0., 0.],
                                [0., 0., self.mdelta],
                                [0., self.mdelta, 0.]])
                                
        self.A5 = np.array([[self.mdelta, self.mdelta, 0.],
                                [self.mdelta, 0., 0.],
                                [0., 0., 0.]])
            
        self.A6 = np.array([[0., self.mdelta, self.mdelta],
                                [self.mdelta, 0., self.mdelta],
                                [self.mdelta, self.mdelta, 0.]])    
                              
        i=0
        self.list = (self.A1, self.A2, self.A3, self.A4, self.A5, self.A6)
        self.determinant = np.array([0.0]*12)
        for matrix in self.list:
            self.matrix = self.base+matrix
            self.negmatrix = self.base-matrix
            self.determinant[2*i] = abs(np.linalg.det(self.matrix))
            self.determinant[2*i+1] = abs(np.linalg.det(self.negmatrix))                        
            i+=1
   
a = Main()
