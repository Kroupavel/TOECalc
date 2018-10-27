# -*- coding: utf-8 -*-
import Tkinter as tk
import numpy as np
import sys

class GUI(tk.Frame):
    """
    This class creates a GUI for producing the deformation matrices and calculating the elastic constants.
    Important information is contained within the README file that should be located in the same file as this script.
    Input information is contained in the input text field and an example input is in the example file. However main points will be repeated:
        Don't mix the numerical data with any string characters or semicolon as it will break the code. Data is separated by space, comma or a new line.
        Dot is used as a decimal point. 
    """
    def __init__(self, root):
    
        # Canvas creator
        tk.Frame.__init__(self)
        self.canvas = tk.Canvas(root)
        self.canvas.pack
        
        # scroller
        self.yscroll = tk.Scrollbar(root, orient = "vertical", command = self.canvas.yview)
        self.canvas.configure(yscrollcommand = self.yscroll.set)
        self.yscroll.pack(side = "right", fill="y")
        
        self.xscroll = tk.Scrollbar(root, orient = "horizontal", command = self.canvas.xview)
        self.canvas.configure(xscrollcommand = self.xscroll.set)
        self.xscroll.pack(side = "bottom", fill = "x")
        
        self.canvas.pack(side ="left", fill="both", expand=True)
        
        # Frame 1 creator
        self.frame1 = tk.Frame(self.canvas, width=750, height=500)
        self.canvas.create_window((0,0), window = self.frame1, anchor="nw",tags="self.frame1")                                  
        self.frame1.bind("<Configure>", self.onFrameConfigure)
        
        # Frame two is created separately due to different function
        # Frame 3 (results) creator
        self.frame3 = tk.Frame(self.canvas, width=750, height=450)        
        self.canvas.create_window((750,0), window = self.frame3, anchor="nw",tags="self.frame3")
        self.frame3.bind("<Configure>", self.onFrameConfigure)
        
        # Frame 4 (matrix) creator
        self.frame4 = tk.Frame(self.canvas, width=750, height=500)        
        self.canvas.create_window((750,450), window = self.frame4, anchor="nw",tags = "self.frame3")
        self.frame4.bind("<Configure>", self.onFrameConfigure)
        
        #Adding widgets to frame 1
        tk.Label(self.frame1, text = "Input:").grid(row=0, column=0, pady=10, sticky="W")
        self.input = tk.DoubleVar()
        tk.Label(self.frame1, text = "Strain for generating matrices:").grid(row=1, column=0, sticky="W")
        self.mdelta = tk.DoubleVar()
        tk.Entry(self.frame1, textvariable=self.mdelta).grid(row=1, column=1, sticky="W")
        tk.Button(self.frame1,text = "Generate matrices", command=self.gen_matrices).grid(row=1, column=3, stick="W")
        
        self.Txt = tk.Text(self.frame1, height=20, width=82, font=15)
        self.Txt.grid(row=2, column=0, columnspan=4, pady=10)
        self.Txt.insert(tk.END, "Please input your data in the following form:(numbers to be separated by Enter or Space)\
                    \nNumber of different strain values (e.g. number of set of measurements)\nEquilibrium volume in angstroem cubed\
                    \nEquilibrium Energy in eV\nAll strain values\nAll Energy values in the order:\
                    \nA1(strain>0), A1(strain<0),A2(strain>0), A2(strain<0),... all twelve values for first strain\
                    \nRepeat for each value of strain\
                    \nPlease only input numbers and decimal dots as any letter will break the number sorting algorithm") 
        
        tk.Button(self.frame1,text = "Check Data", command = self.data_show).grid(row=3, column=0, stick="W")
        
        # adding widgets to frame 3
        tk.Label(self.frame3, text="Output:").grid(row=0, column = 0, pady=5, sticky="W")
        self.Txt2 = tk.Text(self.frame3, height=22, width=60, font=15)
        self.Txt2.insert(tk.END, "Here the results will be displayed.")
        self.Txt2.grid( row=1, column=0, columnspan=4, pady=10)
      
        # adding widgets to frame 4
        tk.Label(self.frame4, text="Matrices:").grid(row=0, column = 0,pady=5, sticky="W")
        self.Txt4 = tk.Text(self.frame4, height=30, width=60, font=15)
        self.Txt4.insert(tk.END, "Here the matrices will be displayed.")
        self.Txt4.grid( row=1, column=0, columnspan=4)
        
        self.fr2=False # Marker whether a frame2 is created
        
      
    def data_crunch(self):
        """
        This method will collect the data and sort through them, preparing them for manual check and further processing.
        The sorting algorithm is as follows:
         The loop searches through the list of numbers adding each number to a list, until it finds a space or enter.\
         When space or enter is found it checks that the list is non empty(problem where multiple spaces or rows)\
         It takes the list and joins all the symbols in the list into single float variable.\
         The variable is added to the self.results list which contains the complete numbers\
         The 'searched' symbols can be easily changed or added.  
            """
        
        self.input = self.Txt.get(1.0, tk.END)
        imax = len(self.input)
        i=0
       
        self.list=[]
        self.results = []
        
    
        for i in range(0, imax):
            if (self.input[i]== " ") or (self.input[i]== "\n") or (self.input[i]==","):
                if len(self.list) !=0:
                    self.list = [x.encode('ascii') for x in self.list]
                    s = "";
                    seq = self.list
                    joined = float(s.join(seq))
                    self.results.append(joined)
                    self.list = []
                    
                else:
                    None
                
            else:                
                self.list.append(self.input[i])
                      
        return self.results
        sys.stdout.flush()
       
    def data_show(self):
        """This method displays the processed data for manual check and prepares it for the calculation
         If any value is changed within the displayed data the calculator will use the alternated value
        Changing ndelta wont change the number of collumns but will calculate with that value"""
        
        self.data = self.data_crunch() # Importing processed data
        tk.Button(self.frame1,text = "Calculate TOECs", command=self.data_calc).grid(row=3, column=1, stick="W")
        
        #Recreating the frame every time new values are inputed. Easier to destroy the whole frame than 
        # to check what neeeds to be recreated
        if self.fr2==True:
            self.frame2.destroy()
            self.fr2 = False 
            
        self.frame2 = tk.Frame(self.canvas,width=750, height=500)
        self.canvas.create_window((0,500), window=self.frame2, anchor="nw", 
                                  tags="self.frame2")
        self.frame2.bind("<Configure>", self.onFrameConfigure)
        self.fr2 = True
        
        tk.Label(self.frame2, text = "The data was interpreted in the following way:"). grid(row=0, column = 0, columnspan=2, sticky="W")
        
        #number of delta values
        tk.Label(self.frame2, text = "Number of strain values used = ").grid(row=1, column=0, sticky="W")
        self.ndeltas = tk.DoubleVar()
        self.ndeltas.set(self.data[0])
        tk.Entry(self.frame2, textvariable = self.ndeltas).grid(row=1, column=1, sticky="W")
        
        # Volume       
        tk.Label(self.frame2, text = "Volume of the crystal in angstroem cubed = ").grid(row=1, column=2, sticky="W" )
        self.volume = tk.DoubleVar()
        self.volume.set(self.data[1])
        tk.Entry(self.frame2, textvariable = self.volume).grid(row = 1, column=3, sticky="W") 
        
        # Eq Energy
        tk.Label(self.frame2, text = "Equilibrium energy in eV = ").grid(row=2, column=0, sticky="W" )
        self.eqenergy = tk.DoubleVar()
        self.eqenergy.set(self.data[2])
        tk.Entry(self.frame2, textvariable = self.eqenergy).grid(row = 2, column=1, sticky="W")
                
        # Strain
        tk.Label(self.frame2, text = "The values of strain = ").grid(row=3, column=0, sticky="W", pady=10)
        n = int(self.data[0])
        
        i = 0        
        self.delt = [0.0]*n
        for i in range (0, n):
            self.delt[i] = tk.DoubleVar()
            self.delt[i].set(self.data[i+3])
            tk.Entry(self.frame2, textvariable=self.delt[i]).grid(row = 3+i, column=1, sticky="W")
     

        # Energies
        tk.Label(self.frame2, text="δ>0").grid(row=4+n, column=1)
        tk.Label(self.frame2, text="δ<0").grid(row=4+n, column=2)
        self.values = [0.0]*((n)*12) #the list to which energies will be stored
         
        # These loops scan over all 6 sets of energies (A1-A6 for both positive and negative strains)
        # for each of the n sets of measurements. It stores all the enrgies in one list and displays them           
        j=0
        self.leaveoo = [0.0]*n 
        for j in range(0, n):   
            i=0                  
            for i in range(0,6):
                
                tk.Label(self.frame2, text="E(A%d)= " %(i+1)).grid(row=6+n+i+(7*j), column=0)
                self.values[(2*i)+(j*12)] = tk.DoubleVar()
                self.values[(2*i)+(j*12)].set(self.data[3+n+(j*12)+2*i])            
                self.entry = tk.Entry(self.frame2, textvariable=self.values[(2*i)+(j*12)]).grid(row=6+n+i+(7*j), column=1)  
                            
                self.values[(2*i)+1+(12*j)] = tk.DoubleVar()
                self.values[(2*i)+1+(j*12)].set(self.data[3+n+(j*12)+(2*i)+1])
                self.entry = tk.Entry(self.frame2, textvariable=self.values[(2*i)+1+(12*j)]).grid(row=6+n+i+(7*j), column=2)                
            
                              
            #if j!=0: #separates the sets of results for better arrangment                
            self.leaveoo[j] = tk.IntVar()
            tk.Checkbutton(self.frame2, text="Omit following data", variable=self.leaveoo[j]).grid(row=5+n+(7*j), column=0, sticky="W")
                          
    def data_calc(self): 
        """This method takes values from the entry fields created in data_show, which can be altered from the original values
         It then feed them and calls Script method which returns the results which are displayed by this method in frame3        
         Create a list of used deltas. It checkes the value of checkbox and if the checkbox is TRUE,
         it doesn't append the delta to a list that will be used for calculation. The same happens for the set of 12 energies corresponding to that delta
        """
        i=0 
        strains = []
        self.omitted = 0
        for i in range (0, int(self.ndeltas.get())):
            if self.leaveoo[i].get()!=1:
                strains.append(self.delt[i].get())
            
            else:
                self.omitted +=1 # We need to have an effective ndeltas to use when calling the Script
       
        j=0 # Getting values of energies
        ener=[]
        for j in range (0, int(12*self.ndeltas.get())):
              if self.leaveoo[int(j/12)].get()!=1: # need to omit the block of energies which are to be omitted. Use the fact that int() round always down
                  ener.append(self.values[j].get())
                       
        energies=np.array(ener)
        
       # ndelta, volume, eqenergy, strains, energy
        a = Script(self.ndeltas.get()-self.omitted, self.volume.get(), self.eqenergy.get(),strains, energies )
        
        self.d, residual = a.solve()
       
         #Final results display        
        constants = np.array(["C11   = ", "C12   = ", "C44   = ", "C111 = ", "C112 = ", "C123 = ", "C144 = ", "C166 = ", "C456 = "])
        self.Txt2.delete(0.0, tk.END)
        for i in range (0,9):
            self.Txt2.insert(tk.END, constants[i])
            self.Txt2.insert(tk.END, self.d[i])      
            self.Txt2.insert(tk.END, " #(GPa)\n") 
        self.Txt2.insert(tk.END, residual[0])    
        self.Txt2.insert(tk.END, " #The residuals")
        
                
        if self.d[2]>0 and self.d[0]> abs(self.d[1]) and self.d[0]+2*self.d[1]>0:
            self.Txt2.insert(tk.END, "\n#This material fulfills the stability criteria for cubic crystal!" ) 
            
        else:
            self.Txt2.insert(tk.END, "\n#This material is unstable.")
        
        sys.stdout.flush()
        
    def gen_matrices(self):
        """
        This method uses Matrix class to generate the deformation matrices for a given value of delta(strain) and then displays them.
        There is a special class for generating the matrices as there is another part of the code that uses the matrices but for a different
        purpose so a basic generator can be for both of these functions. 
        """
        self.Txt4.delete(0.0, tk.END)
        self.delta = self.mdelta.get()
        a = Matrix(self.delta)
        self.base, self.A1, self.A2, self.A3, self.A4, self.A5, self.A6 = a.results() 
                               
        i=1
        self.list2 = (self.A1, self.A2, self.A3, self.A4, self.A5, self.A6)
        for matrix in self.list2:
            self.matrix = self.base+matrix
            self.negmatrix = self.base-matrix
            #The printing format is picked such that it fits the VASP input format.
            self.Txt4.insert(tk.END, "#A%s for delta =%s: \n" %(i, self.mdelta.get()))
            self.Txt4.insert(tk.END,  '\n'.join('  '.join(str(cell) for cell in row) for row in self.matrix))
            self.Txt4.insert(tk.END, "\n#A%s for delta = -%s: \n" %(i, self.mdelta.get()))
            self.Txt4.insert(tk.END, '\n'.join('  '.join(str(cell) for cell in row) for row in self.matrix))
            self.Txt4.insert(tk.END, "\n")                      
            i+=1
               
        
    def onFrameConfigure(self, event):
        """Adjusts the scroll acordingly to changes int the size of the canvas"""
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))


class Script:
    """This method calculates the SOECs and TOECs when given input of the form:
    Number of strain values(delta), Volume, Eq. Energy, all strain values in a list and all the energies in an array
    The energies are in the order A1+, A1-, A2+, A2- for first strain and than repeated for the other values of strain
    The solution is found using the least square method numpy.linalg.lstsq
    Units: energies in eV, volume in Angstroem cubed"""   
    
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
            Script.matr(self, self.ddelta) #The fractional change to the volume of the crystal is produced in the determinant array
            deltasquared = self.delta[i]**2
            deltacubed = abs(self.delta[i]**3)
            
            x=0
            for x in range (0,12):
                count= x+(i*12)                               
                Pos = ((self.energy[count]-self.eqenergy)*160.217662000/(self.volume*self.determinant[x])) # The number coeficient comes from changing the energy to joules                   
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
                     
        newb=np.array(b)    
        amatrix=np.matrix(a)
        x, residues, rank, s = np.linalg.lstsq(amatrix, newb)

        #rearrange the output. This is as the ideal output is in different order than the equations
        dummy=x[1]
        x[1]=x[2]
        x[2]=x[5]
        x[5]=x[4]
        x[4]=x[3]
        x[3]= dummy
        return x, residues
        
    def matr(self, delta): 
        """
        This method calls the Matrix class to generate the apropriate deformation matrices.
        It calculates the determinant of the matrices which is the coefficient by which the volume increase/decrease upon deformation.
        It then creates a list of these "adjustions" of the volume that are used in the main calculation for more accurate results.
        """                     
        self.ddelta = delta
                
        a = Matrix(self.ddelta)
        self.base, self.A1, self.A2, self.A3, self.A4, self.A5, self.A6 = a.results()
        
        i=0
        self.list = (self.A1, self.A2, self.A3, self.A4, self.A5, self.A6)
        self.determinant = np.array([0.0]*12)
        for matrix in self.list:
            self.matrix = self.base+matrix
            self.negmatrix = self.base-matrix
            self.determinant[2*i] = abs(np.linalg.det(self.matrix))
            self.determinant[2*i+1] = abs(np.linalg.det(self.negmatrix))            
            
            i+=1
        
class Matrix():
    """
    This class is an assistant class that generates the deformation matrices for a given value of delta to be used in other parts of the code.    
    This class must be first initialized by A=Matrix(delta) and then matrices are obtained through: a.results()
    """    
    def __init__(self, delta):
        self.mdelta = delta
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
                                
    def results(self):                                                        
        return self.base, self.A1, self.A2, self.A3, self.A4, self.A5, self.A6
                                             

if __name__ == "__main__":
    root=tk.Tk()
    GUI(root).pack(expand=True)
    root.mainloop()