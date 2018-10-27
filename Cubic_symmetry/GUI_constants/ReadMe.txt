This script creates GUI for manual calculations of SOEC and TOEC. It is also capable of generating the deformation matrices. 
However it is not capable of any post processing like calculating Young moduli and pressure derivatives.
To run just write:
python C:\Users\Pavel\Documents\Python\TOECalc\GUI_constants\GUI.py

A window will open. On the top left there is a single small entry field with a "Generate matrices" button next to it.
Write a value of strain into the entry field and press the button. In the bottom right field the matrices will appear in the form suitable for Vasp calculation.

The large text field on the left contains information on how to input data form the main calculation. See example files. 
Example-annotated.txt will not work when copied into the field as it contains non numerical characters. Example-working.txt will work.
Be careful not to have any non-numerical characters, even commas in the input text. Separate values by spaces or new lines.

When you written your appropriate data, press the "Check data" button. A form with entry fields filled with your data will appear. 
You can check that the data have been recognised and digested correctly. If you make any change to the values in the entry fields, the changes will be taken into consideration in calculations.
A new button will appear "Calculate TOECs" If you press it, SOEC and TOEC will appear in the top right field.