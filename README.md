# NZShearWallDesign
A script which takes output tables relating to shear walls from ETABS summarised within an Excel Spreadsheet and runs NZS 3101:2006 checks for in-plane bending, shear, confinement and anti-buckling. Contact Ray Laxmidas (ray.laxmidas@mottmac.com) for any queries or to report any bugs.
The script takes ETABS Outputs summarised in a spreadsheet originally developed by Mitchell Mulvey (mitchell.mulvey@mottmac.com) for Australian shear wall design. Note the associated Python Macro for Australian design is not required for this script to run. The scipt runs NZS 3101:2006 checks based on the spreadsheet developed within Mott MacDonald NZ by Darrin Liddell.
These script have been written to assit in the automation of the design and documentation of the shear walls in the CUMA in Christchurch. Primarily in the West, South and Lower Tiers.

## MainControl.py
The file MainControl.py controls 3 key functions.

(1) Extraction of key design parameters from the "Shear Wall Script Design Spreadsheet" and passing these parameters into the NZShearWallChecks.py files which contains the individual checks. The MainControl.py exports the outputs of the checks into either a console or a .txt file as well summarising the results into an Excel File. 

(2) Controlling some specific parameters which are passed into the NZShearWallChecks.py which are not defined in the Shear Wall Script Design Spreadsheet these include safety factor for Shear Checks which is set at 0.75, the size of stirrup reinforcement for anti-buckling and confinement which is set at 10mm.  

(3) Calculating design forces. The peak tension was combined with the peak moment (of all load cases without compression). 

The MainControl.py was developed specifically for the CUMA project but may be revised for future projects. A key assumption in the MainControl.py which is specific to the CUMA project is that moment capacities are based on the dowel/starter bars (single layered) and then the main reinforcement is sized with a double layer of equivalent area.

## NZShearWallChecks.py
The file NZShearWallChecks.py contains all shear wall design checks NZS 3101:2006. The checks are primarily written for elastic design.

## Notes:
These scripts have been validated through shear walls calculated on Mott MacDonald spreadsheets for the CUMA East Stand carried out by Andrew Wei (andrew.wei@mottmac.com).

## Running of Scripts:

(1) Download the files 

(2) Extract ETABS tables and place into the "Shear Wall Script Design Spreadsheet" excel file.

(3) Open CMD (see example below)

(4) Change directory to file path.

(5) Activate virtual environment.

(6) Run Script

(7) An output file summarising the check results will be outputs. Detailed calculations for Peer Review will be printed to the console or can be saved to DesignCheckOuptuts.txt if line 18 of MainControl.py is not commented out.

Example:

chdir C:\Downloads\NZShearWallDesign

.venv\Scripts\Activate.bat

python.exe MainControl.py
