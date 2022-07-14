"""
Shear Wall Design to NZ Standard
Created by: Ray Laxmidas

The following script extracts the design inputs from the Script Design Spreadsheet
and passes them into the NZ Shear Wall Checks functions. The moments design is based 
on the starter bars / dowels. However this can be adjust to run checks for the main 
reo in the wall.
"""
import sys
import pandas as pd
import NZShearWallChecks as NZSWC
import warnings
warnings.filterwarnings("ignore")

#Define File Path for Ouptuts 
file_path = 'DesignCheckOuptuts.txt'
sys.stdout = open(file_path, "w")

#Importing Data
shear_walls_parameters = pd.read_excel (r'Shear Wall Script Design Spreadsheet.xlsm', 
                    sheet_name='Master', 
                    skiprows=[0,2],
                    usecols = 'C,G:U')
shear_walls_parameters.set_index('Pier')

print("Shear Wall Parameters are shown below:")
print(shear_walls_parameters)

shear_walls_parameters['P'] = ''
shear_walls_parameters['M3'] = ''
shear_walls_parameters['V2'] = ''

shear_wall_forces = pd.read_excel (r'Shear Wall Script Design Spreadsheet.xlsm', 
                    sheet_name='Pier Forces', 
                    skiprows=[0,2],
                    usecols = 'C,F:K')

shear_wall_forces['M3'] = shear_wall_forces['M3'].abs()
shear_wall_forces.set_index('Pier')

#Find Maximum Moments, Shears and Axial Forces for Each Pier
print("Maximum Axial Forces are shown below:")
P_max = shear_wall_forces.loc[shear_wall_forces.groupby('Pier').idxmax()['P']]
P_max.set_index('Pier')
print(P_max)

shear_wall_forces['V2'] = shear_wall_forces['V2'].abs()
print("Maximum Shear Forces are shown below:")
V2_max = shear_wall_forces.loc[shear_wall_forces.groupby('Pier').idxmax()['V2']]
V2_max.set_index('Pier')
print(V2_max)

print("Maximum Moments are shown below:")
#Drop all load cases with compression (which is beneifical to the design) and then find the maximum moment.
#shear_wall_forces.drop(shear_wall_forces[shear_wall_forces['P'] < 0].index, inplace = True)
M3_max = shear_wall_forces.loc[shear_wall_forces.groupby('Pier').idxmax()['M3']]
M3_max.set_index('Pier')
print(M3_max)

i = 0
while i < len(shear_walls_parameters):
    shear_walls_parameters['P'][i] = P_max.iat[i,1]
    shear_walls_parameters['M3'][i] = M3_max.iat[i,6]
    shear_walls_parameters['V2'][i] = V2_max.iat[i,2]
    i = i + 1

#Setup Output Columns:
#Flexural Outputs:
shear_walls_parameters['Maximum longitudinal bar diameter'] = ''
shear_walls_parameters['Starter Dia Check'] = ''
shear_walls_parameters['Maximum spacing of vertical reinforcement'] = ''
shear_walls_parameters['Starter Spacing Check'] = ''
shear_walls_parameters['Vertical reinforcement ratio'] = ''
shear_walls_parameters['Minimum vertical reinforcement ratio'] = ''
shear_walls_parameters['Ratio Check'] = '' 
shear_walls_parameters['Neutral axis depth'] = ''
shear_walls_parameters['M Uti of Starter Reo'] = ''  

#Shear Outputs:
shear_walls_parameters['v_n'] = ''
shear_walls_parameters['Shear Stress Check'] = ''
shear_walls_parameters['V_c'] = ''
shear_walls_parameters['Shear Bar Spacing'] = ''
shear_walls_parameters['Max spacing of horizontal reinforcement'] = ''
shear_walls_parameters['Horizontal Spacing Check for Shear'] = ''
shear_walls_parameters['A_v'] = ''
shear_walls_parameters['A_vmin'] = ''

#Shear Friction of Dowels Ouputs:
shear_walls_parameters['Shear friction uti of dowels'] = ''

#Antibuckling Outputs:
shear_walls_parameters['Area of longitudinal bars restrained by one horizontal stirrup'] = ''
shear_walls_parameters['Vertical spacing of horizontal stirrups'] = ''
shear_walls_parameters['Max vertical spacing of horizontal stirrups'] = ''
shear_walls_parameters['Stirrup Spacing check'] = ''
shear_walls_parameters['A_te'] = ''
shear_walls_parameters['Selected aiameter of stirrup tie for antibuckling'] = ''
shear_walls_parameters['Stirrup dia check'] = ''

#Confinement Ouptus:
shear_walls_parameters['Ash = '] = ''
shear_walls_parameters['Stirrup diameter for confinement'] = ''
shear_walls_parameters['No of stirrup legs for confinement'] = ''
shear_walls_parameters['Stirrup spacing for confinement zones'] = ''

#Run Shear Wall Checks
i = 0
while i < len(shear_walls_parameters):
    print('------------------------------------------------------------------')
    print('-------------------------------',shear_walls_parameters['Pier'][i],'-------------------------------')
    print('-------------------------------------------------------------------')
    #Setup variables
    wall_types=1
    t=shear_walls_parameters['Width'][i]
    l_w=shear_walls_parameters['Depth'][i]
    fc=shear_walls_parameters['fc'][i]
    fy=shear_walls_parameters['fsy'][i]
    fyt=shear_walls_parameters['fsy'][i]
    dbl=shear_walls_parameters['Verti.Dia.'][i]
    sv=shear_walls_parameters['Verti. Spacing'][i]      
    nL=1
    N_max=-1*shear_walls_parameters['P'][i]
    M_max=shear_walls_parameters['M3'][i]
    V_max=shear_walls_parameters['V2'][i]
    cover=shear_walls_parameters['Cover'][i]
    dross=1
    d_bh=shear_walls_parameters['Horiz.Dia.'][i]
    s_h=shear_walls_parameters['Horiz. Spacing'][i]
    
    #Run starter bars for moment checks:   
    flexural_design = NZSWC.flexural_design(wall_types,t,l_w,fc,fy,fyt,dbl,sv,nL,N_max,M_max,cover,dross,d_bh)
    shear_walls_parameters['M Uti of Starter Reo'][i] = shear_walls_parameters['M3'][i]/flexural_design[0]
    shear_walls_parameters['Starter Dia Check'][i] = flexural_design[1]
    shear_walls_parameters['Starter Spacing Check'][i] = flexural_design[2]
    shear_walls_parameters['Ratio Check'][i] = flexural_design[3]
    shear_walls_parameters['Neutral axis depth'][i] = flexural_design[4]
    shear_walls_parameters['Maximum spacing of vertical reinforcement'][i] = flexural_design[5]
    shear_walls_parameters['Maximum longitudinal bar diameter'][i] = flexural_design[6]
    shear_walls_parameters['Vertical reinforcement ratio'][i] = flexural_design[7]
    shear_walls_parameters['Minimum vertical reinforcement ratio'][i] = flexural_design[8]
    
    #Run checks for shear:
    phiV=0.75
    shear_design = NZSWC.shear_design(phiV,V_max,N_max,M_max,d_bh,t,l_w,fc,fyt)
    shear_walls_parameters['Shear Bar Spacing'][i] = shear_design[0]
    shear_walls_parameters['Max spacing of horizontal reinforcement'][i] = shear_design[1]
    shear_walls_parameters['Horizontal Spacing Check for Shear'][i] = shear_design[3]
    shear_walls_parameters['v_n'][i] = shear_design[4]
    shear_walls_parameters['Shear Stress Check'][i] = shear_design[2]
    shear_walls_parameters['V_c'][i] = shear_design[5]
    shear_walls_parameters['A_v'][i] = shear_design[6]
    shear_walls_parameters['A_vmin'][i] = shear_design[7]
    
    #Run checks for Shear Friction of Dowels:
    phi=0.7
    u_sf=1.0
    shear_friction = NZSWC.Shear_friction_of_dowels(phi,u_sf,dbl,l_w,sv,cover,d_bh,fy,N_max,V_max)
    shear_walls_parameters['Shear friction uti of dowels'][i] = shear_walls_parameters['V2'][i]/shear_friction
    
    #Run checks for AntiBuckling
    ds = 10
    anti_buckling = NZSWC.AntiBuckling(s_h,wall_types,t,fy,fyt,dbl,ds)
    shear_walls_parameters['Area of longitudinal bars restrained by one horizontal stirrup'][i] = anti_buckling[0]
    shear_walls_parameters['Vertical spacing of horizontal stirrups'][i] = anti_buckling[1]
    shear_walls_parameters['Max vertical spacing of horizontal stirrups'][i] = anti_buckling[2]
    shear_walls_parameters['Stirrup Spacing check'][i] = anti_buckling[3]
    shear_walls_parameters['A_te'][i] = anti_buckling[4]
    shear_walls_parameters['Selected aiameter of stirrup tie for antibuckling'][i] = anti_buckling[5]
    shear_walls_parameters['Stirrup dia check'][i] = anti_buckling[6]
    
    #Run checks for Confinement
    direction = 1
    c = flexural_design[4]
    s_vmax = flexural_design[5]
    no_legs=2
    confinement = NZSWC.confinement(wall_types,cover,direction,sv,ds,t,fc,c,l_w,s_vmax,fyt,no_legs,s_h)
    shear_walls_parameters['Ash = '][i] = confinement[0]
    shear_walls_parameters['Stirrup diameter for confinement'][i] = confinement[1]
    shear_walls_parameters['No of stirrup legs for confinement'][i] = confinement[2]
    shear_walls_parameters['Stirrup spacing for confinement zones'][i] = confinement[3]
    
    i = i+1
    
shear_walls_parameters.to_excel("ShearWallDesignOuput.xlsx")