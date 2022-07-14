# -*- coding: utf-8 -*-
"""
The following are the functions to check shear walls for in - plane bending and
shear against NZ code. Confinement and anti-buckling are also calculated as to 
NZS 3101:2006.

Contact:
ray.laxmidas@mottmac.com

Inputs Variable Names:
    wall_types =
        1 = Elastic
        2 = Limited Ductile
        3 = Fully Elastic
    t = wall thickness
    l_w = Horizontal lenght of the wall
    fc = Concrete compressive strenght
    fy = Yield strength of longitudinal reinforcement
    fyt = Yield strength of horizontal reinforcement
    dbl = Diameter of longitudinal bars
    sv = Spacing of vertical reinforcement
    nL = number of layers of reinforcement 
    N_max = maximum axial force (where Compression = (-), Tension = (+), kN)
    M_max = maximum major axis bending (kNm)
    dross =
        1 = for dross back calculation
        2 = for main wall reo
    d_bh = diameter of horizontal reinforcement
"""
import math
import numpy as np
import pandas as pd

#Function is used to find the tensile force in each of the bars in the shear wall.
def check_bars(c,c_end,l_w,fy,dbl,nL,b1,a1,t,fc,N_max,sv):
    column_names = ["bar", "bar status", "es", "fs", "Fs", "Fs(d - a/2)"]
    df = pd.DataFrame(np.zeros((100,6)),columns = column_names)
    
    for ind in df.index:
        #Bar Positions
        if ind == 0:
            df.loc[ind]['bar'] = c_end 
        else:
            df.loc[ind]['bar'] = sv + df.loc[ind-1]['bar']
        
        #Set Bar Status
        if df.loc[ind]['bar'] > l_w - c_end:
            df.loc[ind]['bar status'] = -1
        elif df.loc[ind]['bar'] < c:
            df.loc[ind]['bar status'] = 0 #Bar in compression
        else: 
            df.loc[ind]['bar status'] = 1 #Bar in tension
    
        #Setting e.s
        if df.loc[ind]['bar status'] == 0:
           df.loc[ind]['es'] = 0
        elif df.loc[ind]['bar status'] == -1:
           df.loc[ind]['es'] = 0
        else:
           df.loc[ind]['es'] = 0.003*(df.loc[ind]['bar'] - c)/c
       
        #Setting fs
        df.loc[ind]['fs'] = min(df.loc[ind]['es']*200000,fy)
        
        #Setting Fs
        df.loc[ind]['Fs'] = (df.loc[ind]['fs']*(dbl**2)*math.pi*nL)/4000
        
        #Setting Fs(d-a/2)
        df.loc[ind]['Fs(d - a/2)'] = df.loc[ind]['Fs']*(df.loc[ind]['bar']-0.5*(b1*c))/1000
    
    #Calculate total tension in bars
    Ts1 = df['Fs'].sum(skipna=True)   
    Ts2 = df['Fs(d - a/2)'].sum(skipna=True)
    
    #Cclaulate Concrete Compression
    Cc = a1*b1*c*t*fc/1000
    N_max2 = N_max*(0.5*l_w-(0.5*b1*c))/1000
    
    #Target for goal seek
    delta = Ts1 + N_max - Cc
    
    #Moment Capacity
    M_capacity = 0.85*(Ts2 + N_max2)
    
    return [delta,M_capacity] 

"""
The following function calculatures the flexural capacity of the shear walls.
Outputs:
    Moment Capacity
    Bar Diameter Check (Ture = OKAY, False = Not OKAY)
    Spacing Check (True = OKAY, False = Not OKAY)
    Ratio Check (True = OKAY, False = Not OKAY)     
"""
def flexural_design(wall_types,t,l_w,fc,fy,fyt,dbl,sv,nL,N_max,M_max,cover,dross,d_bh): 
    print('-------------------------Flexural Checks--------------------------')
    print('-------------------------Input parameters-------------------------')
    print('Wall thickness = ', t, 'mm')
    print('Horizontal lenght of the wall = ', l_w, 'mm')
    print('fc = ', fc, 'MPa')
    print('fy = ', fy, 'MPa')
    print('Diameter of longitudinal bars = ', dbl, 'mm')
    print('Diameter of horizontal reinforcement = ', d_bh, 'mm')
    print('Spacing of vertical reinforcement = ', sv, 'mm')
    print('Number of layers of reinforcement =', nL, 'Layers')
    print('N*', N_max, 'kN')
    print('M*', M_max, 'kNm')
    
    if dross == 1:
        print('Calculations are for starter bars')
    
    print('----------------------Spacing and Reo Checks-----------------------')
    #Calculate the maixmum longituindal bar diameter
    if wall_types == 1: 
        dbl_max = t/7
    elif wall_types == 2:
        dbl_max = t/8
    else:
        dbl_max = t/10
    print('Maximum bar diameter =' ,round(dbl_max,2),'mm')
    
    #Checking the bar diameter
    if dbl_max > dbl:
        bar_dia_check = 1 #Bar diameter OKAY
        print('Selected Bar Diameter is OKAY')
    else:
        bar_dia_check = 0 #Bar diameter NOT OKAY
        print('Selected Bar Diameter is NOT OKAY')
    
    #Calculate the maximum spacing
    sv_max = min(450,3*t)
    print('Maximum bar spacing =' ,round(sv_max,0),'mm')
    
    #Checking the spacing
    if sv > sv_max:
        spacing_check = 0 #Spacing NOT OKAY
        print('Selected Bar Spacing is NOT OKAY')
    else: 
        spacing_check = 1 #Spacing OKAY
        print('Selected Bar Spacing is OKAY')
    
    #Calculate the cover to end bar, each end
    if dross == 1:
        c_end = cover + d_bh + 0.5*dbl + 0.5*sv
    else:
        c_end = cover + d_bh #Main reo
    
    #Calculate wall reinforcement area
    As = nL*dbl**2*math.pi*l_w/(4*sv)
    print('Area of wall reo =',  round(As,0), 'mm^2')
    
    #Calculate vertical reinforcement ratio
    pn = As / (t*l_w)
    print('Reo Rate =',round(pn,6))
    
    #Calculate minimum reinforcment ratio
    pn_min = max(math.sqrt(fc)/(4*fy),0.7/fy,0.0014)
    print('Minimum Reo Rate =',round(pn_min,6))
    
    #Checking reinforcement ratio
    if pn > pn_min:
        ratio_check = 1 #Ratio OKAY
        print('Reo Rate OKAY')
    else:
        ratio_check = 0 #Ratio NOT OKAY
        print('Reo Rate NOT OKAY')
        
    #Moment capacity
    if fc < 55:
        a1 = 0.85
    elif fc < 80:
        a1 = 0.85 - 0.004*(fc - 55)
    else:
        a1 = 0.75
    
    if fc < 30:
        b1 = 0.85
    elif fc <= 55:
        b1 = 0.85 - 0.008*(fc - 30)
    else:
        b1 = 0.65
  
    c = 1
    result = check_bars(c,c_end,l_w,fy,dbl,nL,b1,a1,t,fc,N_max,sv)
    lower = 0
    upper = 10000
    
    #Finding the value of the compression block using binary search.
    while abs(result[0])>= 10:
        if result[0] < 0:     
            upper = c
            c = (lower + upper)/2
            result = check_bars(c,c_end,l_w,fy,dbl,nL,b1,a1,t,fc,N_max,sv)
        if result[0] > 0:
            lower = c
            c = (lower + upper)/2
            result = check_bars(c,c_end,l_w,fy,dbl,nL,b1,a1,t,fc,N_max,sv)
    
    print('-------------------------------Results-------------------------------')
    print('Depth to Neutral Axis:',round(c,1), 'mm')    
    result = check_bars(c,c_end,l_w,fy,dbl,nL,b1,a1,t,fc,N_max,sv)
    
    print('The value In-Plane Moment Capacity is:', round(result[1],0), 'kNm')
    M_capacity = result[1]       
    
    print('Moment Utilisation:', round((M_max/result[1])*100,0),'%')
          
    return M_capacity, bar_dia_check, spacing_check, ratio_check,c,sv_max,dbl_max,pn,pn_min 

"""
The following function calculatures the flexural capacity of the shear walls.
Outputs:
    Spacing of horizontal reinforcement
    Shear Stress
    Spacing Check (True = OKAY, False = Not OKAY)    
"""
def shear_design(phi,V_max,N_max,M_max,dbh,t,l_w,fc,f_yt):
    print('----------------------------Shear Checks----------------------------')
    print('--------------------------Input parameters--------------------------')
    print("Maximum Shear = ", V_max,"kN")
    print("Horizontal Bar Diameter = ", dbh,"mm")
    print('-------------------------------Outputs------------------------------')
    vn = 1000*V_max/(phi*t*0.8*l_w)
    print('v.n =', round(vn,2), 'MPa')  
    
    #Concrete Shear Capacity:
    vc = M_max/V_max - (l_w/2000)
    
    
    #Equation 11-14
    EQ11_14 = 0.27*math.sqrt(fc)+(1000*N_max)/(4*l_w*t)
    
    #Equation 11-15
    EQ11_15 = 0.05*math.sqrt(fc)+(l_w/1000*(0.1*math.sqrt(fc)+200*N_max/(l_w*t))/vc) 
    
    if(vc > 0):
        vc = min(EQ11_14,EQ11_15)
    else:
        vc = EQ11_14
    
    V_c = vc*t*0.80*l_w/1000
    print('V.c =', round(V_c,2), 'kN')  
    
    if(vn<min(0.2*fc,8)):
        shear_stress_check = 1
        print("Shear Stress OKAY")
    else:
        shear_stress_check = 0
        print("Stress Stress NOT OKAY")
    
    A_v = (V_max/phi-V_c)*1000000/(0.8*l_w*f_yt)
    A_vmin =0.7*t*1000/f_yt
    
    s_2 = 2*(dbh**2)*math.pi*1000/(max(A_v,A_vmin)*4)
    print('Spacing of horizontal reinforcement =', round(s_2,0), 'mm') 
    s_2max = min(l_w/5,3*t,450)
    
    
    if (s_2 < s_2max):
        spacing_check = 1
        print("Shear Bar Spacing OKAY")
    else:
        spacing_check = 0
        print("Max spacing governs")
    
    return s_2, s_2max, shear_stress_check, spacing_check, vn, V_c, A_v, A_vmin

#The following function checks the shear friction for dowels.
def Shear_friction_of_dowels (phi,u_sf,dbl,l_w,sv,cover,d_bh,fy,N_max,V_max):  
    print("-----------------Checking Shear Friction of Dowel-----------------")
    
    c_end = cover + d_bh + 0.5*dbl + 0.5*sv
    A_s =(dbl**2)/4*math.pi*(((l_w-2*c_end)/sv)+1)
    phiV_sf = phi*u_sf*(A_s*fy+N_max)/1000
    
    print("PhiV_sf = ", round(phiV_sf,0), "kN")
    
    if phiV_sf > V_max:
        print("OKAY")
    else:
        print("NOT OKAY")
        
    return phiV_sf

#The following function checks the antibuckling stirrups.
def AntiBuckling (s_h,wall_type,t,fy,fyt,dbl,ds):
    print('------------------------Antibuckling Checks------------------------')
    
    print("Selected diameter of stirrup tie: ", ds, "mm")
    print("Selected stirrup spacing: ", s_h, "mm")
    
    A_b = 1 * ((dbl**2)/(4))*math.pi
    print("Area of longitudinal bars restrained by one horizontal stirrup:", round(A_b,0), "mm^2")
    
    if wall_type == 3:
        sh_max = min(6*dbl,t/2)
    else: 
        sh_max = min(10*dbl,300)
    
    print("Maximum stirrup spacing: ", sh_max, "mm")    
    
    A_te = (A_b*s_h*fy)/(96*fyt*dbl)
    
    print("Area of one stirrup tie: ", round(A_te,1), "mm^2")
    
    if (s_h <= sh_max):
        stirrup_spacing_check = 1
        print("Stirrup spacing is OKAY")
    else: 
        stirrup_spacing_check = 0
        print("Stirrup spacing is NOT OKAY")
        
    if(ds**2*math.pi/4 > A_te): 
        stirrup_diameter_check = 1    
        print("Stirrup diameter OK")
    else:    
        stirrup_diameter_check = 0
        ("Revise anti-buckling design")
    return A_b,s_h,sh_max,stirrup_spacing_check,A_te,ds,stirrup_diameter_check 

#The following function checks the confinement at the end zones. 
def confinement(wall_type,cover,direction,sv,ds,t,fc,c,l_w,s_vmax,fyt,no_legs,sh):
    print('------------------------Confinement Checks-------------------------')
    
    print("Selected diameter of stirrup tie: ", ds, "mm")
    if direction == 1: #if in major direction
        h_dash=t-2*cover
    else: #if in minor direction
        h_dash=max(t/2,sv+2*ds)
    
    print("h'' =",h_dash,"mm")
    
    if wall_type == 3:
        A_sh = 0.25*sh*h_dash*t*fc*(c/l_w-0.07)/((t-2*cover)*fyt)
    else:
        A_sh = 0.175*sh*h_dash*t*fc*(c/l_w-0.07)/((t-2*cover)*fyt)
                
    print("A_sh = ", round(A_sh,2), 'mm^2') 
    
    if A_sh < 0:
        stirrup_spacing = sh
    else:
        stirrup_spacing = min(sh,sh*(ds**2)*no_legs*math.pi/(4*A_sh))
    
    print("Stirrup spacing must be <= ", stirrup_spacing, 'mm for confinement zones')  
    return A_sh, ds, no_legs, stirrup_spacing
'''
#Testing outputs for EPO1 hand checked on East Stand (All Ouptuts Correct)
flexural_design = flexural_design(
    wall_types=1,
    t=300,
    l_w=8000,
    fc=40,
    fy=500,
    fyt=500,
    dbl=25,
    sv=200,         
    nL=1,
    N_max=-1544,
    M_max=21575,
    cover=35,
    dross=1,
    d_bh=16)

shear_design(phi=0.75,
             V_max=6403,
             N_max=-1544,
             M_max=21575,
             dbh=16,
             t=300,
             l_w=8000,
             fc=40,
             f_yt=500)

Shear_friction_of_dowels (phi=0.7,
                          u_sf=1.0,
                          dbl=25,
                          l_w=8000,
                          sv=200,
                          cover=35,
                          d_bh=16,
                          fy=500,
                          N_max=-1544,
                          V_max=6403)

AntiBuckling (s_h = 250,
              wall_type = 1,
              t = 300,
              fy = 500,
              fyt = 500,
              dbl = 25,
              ds = 10)

confinement(wall_type = 1,
            cover= 35,
            direction = 1,
            sv = 200,
            ds = 10,
            t = 300,
            fc = 40,
            c = 845.68,
            l_w =8000,
            s_vmax = 450,
            fyt = 500,
            no_legs = 2,
            sh = 250)
'''