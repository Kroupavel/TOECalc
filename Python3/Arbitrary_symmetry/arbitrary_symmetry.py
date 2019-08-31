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
            self.inp = input("Input file:")
            self.outfile = input("Output file(Select or create a folder to which data is stored):") 

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
        self.ndeltas = self.data[0] # number of data sets (e.g. one data set for each value of strain)                   
        self.volume = self.data[1] # volume of the crystal
        self.eqenergy = self.data[2] # the zero stress energy
        n=int(self.data[0])
        # Rework the value of strain into needed form.
        i = 0        
        self.delt = [0.0]*n # here we collect the strain values into one array
        for i in range (0, n):                 
            self.delt[i] = self.data[i+3]
                 
        # Energies are prepared for the input format.    
        self.values = [0.0]*((n)*28) #the list to which energies will be stored 
        # pk4: here we change 12 to 28 as there are 14 deformation matrices and for each we have positive and negative stress        
        # inports all the energy values and save them into self.values:      
        for j in range(0, len(self.data) -3-n ):          
            self.values[j] = self.data[j +3+n]
        #print(" the energy values are %s", self.values)        
        a = Constants(self.ndeltas, self.volume, self.eqenergy,self.delt, self.values )        
        self.results, self.residual = a.solve()
        
        #Storing algorithm
        self.completeName = os.path.join(self.outpath, "Constants.txt")
        f = open(self.completeName, 'w') 
        #pk4: here need to change the constants name and their number, there will be 20 of them    
        
        constants = ["  #C11 in GPa","  #C12 in GPa","  #C13 in GPa","  #C14 in GPa","  #C33 in GPa", "  #C44 in GPa", "  #C111 in GPa", "  #C112 in GPa","  #C113 in GPa","  #C114 in GPa", "  #C123 in GPa","  #C124 in GPa","  #C133 in GPa","  #C134 in GPa", "  #C144 in GPa","  #C155 in GPa","  #C222 in GPa","  #C333 in GPa","  #C344 in GPa", "  #C444 in GPa"]
        for i in range (0,20):
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
        self.n = int(self.data[0])
        self.volume = self.data[1]
        self.eqenergy = self.data[2]
        self.deltas = self.data[3:self.n+3]
        self.energy = self.data[self.n+3:]
        self.values = np.array([0.0]*(int(self.n-1)*28))
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
            for j in range (0,  int(self.n*28)):
                if int(j/28)!= i:
                    self.values[j-help2] = self.energy[j]
                else:
                    help2=28
               
            a = Constants((self.n-1), self.volume, self.eqenergy,self.delt, self.values )
            self.results, residual = a.solve()
            
            #pk4: here again we will have to change the constants      
            f2.write(" When the set number %d was omitted the following data was obtained \n" %i)
            constants = ["  #C11 in GPa","  #C12 in GPa","  #C13 in GPa","  #C14 in GPa","  #C33 in GPa", "  #C44 in GPa", "  #C111 in GPa", "  #C112 in GPa","  #C113 in GPa","  #C114 in GPa", "  #C123 in GPa","  #C124 in GPa","  #C133 in GPa","  #C134 in GPa", "  #C144 in GPa","  #C155 in GPa","  #C222 in GPa","  #C333 in GPa","  #C344 in GPa", "  #C444 in GPa"]
            for i in range (0,20):
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
      #print("Delta values", self.delta)
      #print("volume", self.volume)
      #print("eq energy", self.eqenergy)
      #print("energy", self.energy)
      #print("number of deltas", self.ndeltas)
    
    def solve(self):
        """
        This is the method responsible for the calculation of SOEC and TOEC. When feeded appropriate infromation by its __init__ it
        creates arrays and matrix necessary for least square method solution.
        Mathematically: A*X=B Where X is an 1-D array of the 20 elastic constants
        A is a n*20 matrix where n is equal to 28*number of delta values used. So it is always overdefined
        B is a 1-D array of energy density difference.
        The values in the X array are being calculated for given A and B
        """
        i=0
        a=[] # coefficient matrix
        b = [] # right side matrix (inhomogenous part)
        
        #Iterates over every strain value
        for i in range (0,int(self.ndeltas)):
            self.ddelta = self.delta[i]
            Constants.matrix(self)
            dltsquared = (self.delta[i]**2)/2
            dltcubed = abs(self.delta[i]**3)/6
            
            #calculating the energy density, which is the result vector B
            x=0
            for x in range (0,28):
                count= x+(i*28)  
                # pk4: make sure the energy calc is corret, but should be             
                b.append((self.energy[count]-self.eqenergy)*160.217662000/(self.volume*self.determinant[x])) # changing the energy to joules                   
    
                
                    
        # matrix A with coefficients for matrix x with: C11, C111, C12, C112, C123, C44, C144, C166, C456        
            a.extend([[dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, dltcubed, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, -dltcubed, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [2*dltsquared, 2*dltsquared, 0.0, 0.0, 0.0, 0.0, 4*dltcubed, 6*dltcubed, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2*dltcubed, 0.0, 0.0, 0.0],
                [2*dltsquared, 2*dltsquared, 0.0, 0.0, 0.0, 0.0, -4*dltcubed, -6*dltcubed, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, +2*dltcubed, 0.0, 0.0, 0.0],
                # c11, c12, c13, c14, c33, c44, c111, c112, c113, c114, c123, c124, c133, c134, c144, c155, c222, c333, c344, c444
                [2*dltsquared, 2*dltsquared, 4*dltsquared, 0.0, dltsquared, 0.0, 4*dltcubed, 6*dltcubed, 6*dltcubed, 0.0, 6*dltcubed, 0.0, 6*dltcubed, 0.0, 0.0, 0.0, -2*dltcubed, dltcubed, 0.0, 0.0], 
                [2*dltsquared, 2*dltsquared, 4*dltsquared, 0.0, dltsquared, 0.0, -4*dltcubed, -6*dltcubed, -6*dltcubed, 0.0, -6*dltcubed, 0.0, -6*dltcubed, 0.0, 0.0, 0.0, 2*dltcubed, -dltcubed, 0.0, 0.0],                 
                [dltsquared, 0.0, 0.0, 4*dltsquared, 0.0, 4*dltsquared, dltcubed, 0.0, 0.0, 6*dltcubed, 0.0, 0.0, 0.0, 0.0, 6*dltcubed, 0.0, 0.0, 0.0, 0.0, 8*dltcubed],
                [dltsquared, 0.0, 0.0, 4*dltsquared, 0.0, 4*dltsquared, -dltcubed, 0.0, 0.0, -6*dltcubed, 0.0, 0.0, 0.0, 0.0, -6*dltcubed, 0.0, 0.0, 0.0, 0.0, -8*dltcubed],                
                [dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dltcubed, 0.0, 0.0, 0.0],
                [dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -dltcubed, 0.0, 0.0, 0.0],                
                [0.0, 0.0, 0.0, 0.0, dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dltcubed, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -dltcubed, 0.0, 0.0],                
                [dltsquared, 0.0, 2*dltsquared, 0.0, dltsquared, 0.0, 0.0, 0.0, 3*dltcubed, 0.0, 0.0, 0.0, 3*dltcubed, 0.0, 0.0, 0.0, dltcubed, dltcubed, 0.0, 0.0],                
                [dltsquared, 0.0, 2*dltsquared, 0.0, dltsquared, 0.0, 0.0, 0.0, -3*dltcubed, 0.0, 0.0, 0.0, -3*dltcubed, 0.0, 0.0, 0.0, -dltcubed, -dltcubed, 0.0, 0.0],                                
                [0.0, 0.0, 0.0, 0.0, 0.0, 4*dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8*dltcubed],
                [0.0, 0.0, 0.0, 0.0, 0.0, 4*dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8*dltcubed],
                [2*dltsquared, -2*dltsquared, 0.0, 0.0, dltsquared, 0.0, 0.0, 0.0, 6*dltcubed, 0.0, -6*dltcubed, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dltcubed, 0.0, 0.0],
                [2*dltsquared, -2*dltsquared, 0.0, 0.0, dltsquared, 0.0, 0.0, 0.0, -6*dltcubed, 0.0, 6*dltcubed, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -dltcubed, 0.0, 0.0],
                [dltsquared, 0.0, 0.0, -4*dltsquared, 0.0, 4*dltsquared, 0.0, 0.0, 0.0, -6*dltcubed, 0.0, -12*dltcubed, 0.0, 0.0, 0.0, 12*dltcubed, dltcubed, 0.0, 0.0, 8*dltcubed],
                [dltsquared, 0.0, 0.0, -4*dltsquared, 0.0, 4*dltsquared, 0.0, 0.0, 0.0, +6*dltcubed, 0.0, +12*dltcubed, 0.0, 0.0, 0.0, -12*dltcubed, -dltcubed, 0.0, 0.0, -8*dltcubed],                
                [dltsquared, 0.0, 0.0, 0.0, 0.0, 4*dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12*dltcubed, 0.0, dltcubed, 0.0, 0.0, 0.0],                 
                [dltsquared, 0.0, 0.0, 0.0, 0.0, 4*dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12*dltcubed, 0.0, -dltcubed, 0.0, 0.0, 0.0],                 
                [0.0, 0.0, 0.0, 0.0, dltsquared, 4*dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dltcubed, 12*dltcubed, 0.0],                                 
                [0.0, 0.0, 0.0, 0.0, dltsquared, 4*dltsquared, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -dltcubed, -12*dltcubed, 0.0],
                [dltsquared, 0.0, 2*dltsquared, 4*dltsquared, dltsquared, 4*dltsquared, dltcubed, 0.0, 3*dltcubed, 6*dltcubed, 0.0, 0.0, 3*dltcubed, 12*dltcubed, 12*dltcubed, 0.0, 0.0, dltcubed, 12*dltcubed, 8*dltcubed],                                                 
                [dltsquared, 0.0, 2*dltsquared, 4*dltsquared, dltsquared, 4*dltsquared, dltcubed, 0.0, -3*dltcubed, -6*dltcubed, 0.0, 0.0, -3*dltcubed, -12*dltcubed, -12*dltcubed, 0.0, 0.0, -dltcubed, -12*dltcubed, -8*dltcubed],               
                [2*dltsquared, 2*dltsquared, 0.0, 0.0, 0.0, 4*dltsquared, 4*dltcubed, 6*dltcubed, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12*dltcubed, 12*dltcubed, -2*dltcubed, 0.0, 0.0, 0.0],                
                [2*dltsquared, 2*dltsquared, 0.0, 0.0, 0.0, 4*dltsquared, -4*dltcubed, -6*dltcubed, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12*dltcubed, -12*dltcubed, 2*dltcubed, 0.0, 0.0, 0.0]])
                     
        newb = np.array(b)    
        amatrix = np.matrix(a)
        x, residues, rank, s = np.linalg.lstsq(amatrix, newb, rcond=None)

        return x, residues
        
        
    def matrix(self):
        """
        This method generates the deformation matirces.
        The determinant is calculated and used as a correction to the volume of the crystal.
        (As the volume of the crystal change from the equilibirium value when deformed and it is important to incorporate this change in the script)"""
        #will need more and different matrices
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
                                
        self.A7 = np.array([[0., 0., 0.],
                                [0., self.mdelta, 0.],
                                [0., 0., 0.]])  
                                
        self.A8 = np.array([[0., 0., 0.],
                                [0., 0., 0.],
                                [0., 0., self.mdelta]]) 
                                
        self.A9 = np.array([[0., 0., 0.],
                                [0., self.mdelta, 0.],
                                [0., 0., self.mdelta]]) 
                                
        self.A10 = np.array([[0., 0., 0.],
                                [0., 0., self.mdelta],
                                [0., self.mdelta, 0.]])
                                
        self.A11 = np.array([[0., self.mdelta, 0.],
                                [self.mdelta, 0., 0.],
                                [0., 0., self.mdelta]])
                                
        self.A12 = np.array([[0., 0., 0.],
                                [0., self.mdelta, self.mdelta],
                                [0., self.mdelta, 0.]])   
                                
        self.A13 = np.array([[0., 0., self.mdelta],
                                [0., self.mdelta, 0.],
                                [self.mdelta, 0., 0.]]) 
                                
        self.A14 = np.array([[0., 0., self.mdelta],
                                [0., 0., 0.],
                                [self.mdelta, 0., self.mdelta]]) 
                                
        self.A15 = np.array([[self.mdelta, 0., 0.],
                                [0., 0., self.mdelta],
                                [0., self.mdelta, self.mdelta]]) 
                                
        self.A16 = np.array([[self.mdelta, 0., self.mdelta],
                                [0., self.mdelta, 0.],
                                [self.mdelta, 0., 0.]])                         
                              
        i=0
        self.list = (self.A1, self.A2, self.A3, self.A4, self.A7, self.A8, self.A9, self.A10, self.A11, self.A12, self.A13, self.A14, self.A15, self.A16)
        self.determinant = np.array([0.0]*28)
        for matrix in self.list:
            self.matrix = self.base+matrix
            self.negmatrix = self.base-matrix
            self.determinant[2*i] = abs(np.linalg.det(self.matrix))
            self.determinant[2*i+1] = abs(np.linalg.det(self.negmatrix))                        
            i+=1
   
a = Main()