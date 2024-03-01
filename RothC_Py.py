######################################################################################################################
# 
#  RothC python version    
#
#  This python version was translated from the Fortran code by Alice Milne, Jonah Prout and Kevin Coleman 29/02/2024
#
#  The Rothamsted Carbon Model: RothC
#  Developed by David Jenkinson and Kevin Coleman
#
#  INPUTS: 
#
#  clay:  clay content of the soil (units: %)
#  depth: depth of topsoil (units: cm)
#  IOM: inert organic matter (t C /ha)
#  nsteps: number of timesteps 
#
#  year:    year
#  month:   month (1-12)
#  modern:   %modern 
#  TMP:      Air temperature (C)
#  Rain:     Rainfall (mm)
#  Evap:     open pan evaporation (mm)
#  C_inp:    carbon input to the soil each month (units: t C /ha)
#  FYM:      Farmyard manure input to the soil each month (units: t C /ha)
#  PC:       Plant cover (0 = no cover, 1 = covered by a crop)
#  DPM/RPM:  Ratio of DPM to RPM for carbon additions to the soil (units: none)
#
#  OUTPUTS:
#
#  All pools are carbon and not organic matter 
#
#  DPM:   Decomposable Plant Material (units: t C /ha)
#  RPM:   Resistant Plant Material    (units: t C /ha)
#  Bio:   Microbial Biomass           (units: t C /ha)
#  Hum:   Humified Organic Matter     (units: t C /ha)
#  IOM:   Inert Organic Matter        (units: t C /ha)
#  SOC:   Soil Organic Matter / Total organic Matter (units: t C / ha)
#
#  DPM_Rage:   radiocarbon age of DPM
#  RPM_Rage:   radiocarbon age of RPM
#  Bio_Rage:   radiocarbon age of Bio
#  HUM_Rage:   radiocarbon age of Hum
#  Total_Rage: radiocarbon age of SOC (/ TOC)
#
#  SWC:       soil moisture deficit (mm per soil depth)
#  RM_TMP:    rate modifying fator for temperature (0.0 - ~5.0)
#  RM_Moist:  rate modifying fator for moisture (0.0 - 1.0)
#  RM_PC:     rate modifying fator for plant retainment (0.6 or 1.0)

######################################################################################################################

import pandas as pd
import numpy as np

# Calculates the rate modifying factor for temperature (RMF_Tmp)
def RMF_Tmp (TEMP):
         
    if(TEMP<-5.0):
        RM_TMP=0.0
    else:
        RM_TMP=47.91/(np.exp(106.06/(TEMP+18.27))+1.0)
   
    return RM_TMP

# Calculates the rate modifying factor for moisture (RMF_Moist)
def RMF_Moist (RAIN, PEVAP, clay, depth, PC, SWC):
        
    RMFMax = 1.0
    RMFMin = 0.2

# calc soil water functions properties
    SMDMax=-(20+1.3*clay-0.01*(clay*clay))
    SMDMaxAdj = SMDMax * depth / 23.0
    SMD1bar = 0.444 * SMDMaxAdj
    SMDBare = 0.556 * SMDMaxAdj
      
    DF = RAIN - 0.75 * PEVAP

    
    minSWCDF=np.min (np.array([0.0, SWC[0]+DF]))
    minSMDBareSWC=np.min (np.array([SMDBare, SWC[0]]))
      
    if(PC==1):
        SWC[0] = np.max(np.array([SMDMaxAdj, minSWCDF]))
    else:
        SWC[0] = np.max(np.array([minSMDBareSWC, minSWCDF]))
       
      
    if(SWC[0]>SMD1bar): 
        RM_Moist = 1.0
    else:
        RM_Moist = (RMFMin + (RMFMax - RMFMin) * (SMDMaxAdj - SWC[0]) / (SMDMaxAdj - SMD1bar) )   
    
    return RM_Moist

# Calculates the plant retainment modifying factor (RMF_PC)
def RMF_PC (PC):
     
    if (PC==0):
        RM_PC = 1.0
    else:
        RM_PC = 0.6

    return RM_PC

# Calculates the decomposition and radiocarbon
def decomp(timeFact, DPM,RPM,BIO,HUM, IOM, SOC, DPM_Rage, RPM_Rage, \
           BIO_Rage, HUM_Rage, Total_Rage, modernC, RateM, clay, C_Inp, FYM_Inp, DPM_RPM):

    zero = 0e-8
# rate constant are params so don't need to be passed
    DPM_k = 10.0
    RPM_k = 0.3
    BIO_k = 0.66
    HUM_k = 0.02 
    
    conr = np.log(2.0) / 5568.0

    tstep = 1.0/timeFact    # monthly 1/12, or daily 1/365  
      
    exc = np.exp(-conr*tstep) 
      
 
# decomposition
    DPM1 = DPM[0] * np.exp(-RateM*DPM_k*tstep)
    RPM1 = RPM[0] * np.exp(-RateM*RPM_k*tstep)      
    BIO1 = BIO[0] * np.exp(-RateM*BIO_k*tstep)      
    HUM1 = HUM[0] * np.exp(-RateM*HUM_k*tstep) 

      
    DPM_d = DPM[0] - DPM1
    RPM_d = RPM[0] - RPM1      
    BIO_d = BIO[0] - BIO1
    HUM_d = HUM[0] - HUM1 
      
    
    x=1.67*(1.85+1.60*np.exp(-0.0786*clay))
                    
# proportion C from each pool into CO2, BIO and HUM      
    DPM_co2 = DPM_d * (x / (x+1))
    DPM_BIO = DPM_d * (0.46 / (x+1))
    DPM_HUM = DPM_d * (0.54 / (x+1))
      
    RPM_co2 = RPM_d * (x / (x+1))
    RPM_BIO = RPM_d * (0.46 / (x+1))
    RPM_HUM = RPM_d * (0.54 / (x+1))    
      
    BIO_co2 = BIO_d * (x / (x+1))
    BIO_BIO = BIO_d* (0.46 / (x+1))
    BIO_HUM = BIO_d * (0.54 / (x+1))
      
    HUM_co2 = HUM_d * (x / (x+1))
    HUM_BIO = HUM_d * (0.46 / (x+1))
    HUM_HUM = HUM_d * (0.54 / (x+1))  
           
      
# update C pools  
    DPM[0] = DPM1
    RPM[0] = RPM1
    BIO[0] = BIO1 + DPM_BIO + RPM_BIO + BIO_BIO + HUM_BIO
    HUM[0] = HUM1 + DPM_HUM + RPM_HUM + BIO_HUM + HUM_HUM    
      
# split plant C to DPM and RPM 
    PI_C_DPM = DPM_RPM / (DPM_RPM + 1.0) * C_Inp
    PI_C_RPM =     1.0 / (DPM_RPM + 1.0) * C_Inp

# split FYM C to DPM, RPM and HUM 
    FYM_C_DPM = 0.49*FYM_Inp
    FYM_C_RPM = 0.49*FYM_Inp      
    FYM_C_HUM = 0.02*FYM_Inp   
      
# add Plant C and FYM_C to DPM, RPM and HUM   
    DPM[0] = DPM[0] + PI_C_DPM + FYM_C_DPM
    RPM[0] = RPM[0] + PI_C_RPM + FYM_C_RPM  
    HUM[0] = HUM[0] + FYM_C_HUM
      
# calc new ract of each pool      
    DPM_Ract = DPM1 *np.exp(-conr*DPM_Rage[0])
    RPM_Ract = RPM1 *np.exp(-conr*RPM_Rage[0]) 
      
    BIO_Ract = BIO1 *np.exp(-conr*BIO_Rage[0])
    DPM_BIO_Ract = DPM_BIO * np.exp(-conr*DPM_Rage[0])
    RPM_BIO_Ract = RPM_BIO * np.exp(-conr*RPM_Rage[0])
    BIO_BIO_Ract = BIO_BIO * np.exp(-conr*BIO_Rage[0])
    HUM_BIO_Ract = HUM_BIO * np.exp(-conr*HUM_Rage[0])
      
    HUM_Ract = HUM1 *np.exp(-conr*HUM_Rage[0])   
    DPM_HUM_Ract = DPM_HUM * np.exp(-conr*DPM_Rage[0])
    RPM_HUM_Ract = RPM_HUM * np.exp(-conr*RPM_Rage[0])
    BIO_HUM_Ract = BIO_HUM * np.exp(-conr*BIO_Rage[0])
    HUM_HUM_Ract = HUM_HUM * np.exp(-conr*HUM_Rage[0])
      
    IOM_Ract = IOM[0] *np.exp(-conr*IOM_Rage[0]) 
      
# assign new C from plant and FYM the correct age   
    PI_DPM_Ract = modernC * PI_C_DPM
    PI_RPM_Ract = modernC * PI_C_RPM
      
    FYM_DPM_Ract = modernC * FYM_C_DPM
    FYM_RPM_Ract = modernC * FYM_C_RPM
    FYM_HUM_Ract = modernC * FYM_C_HUM          
      
# update ract for each pool        
    DPM_Ract_new = FYM_DPM_Ract + PI_DPM_Ract + DPM_Ract*exc
    RPM_Ract_new = FYM_RPM_Ract + PI_RPM_Ract + RPM_Ract*exc    
    
    BIO_Ract_new = (BIO_Ract + DPM_BIO_Ract + RPM_BIO_Ract + 
                    BIO_BIO_Ract + HUM_BIO_Ract )*exc
      
    HUM_Ract_new = FYM_HUM_Ract + (HUM_Ract + DPM_HUM_Ract + 
                                   RPM_HUM_Ract + BIO_HUM_Ract + HUM_HUM_Ract )*exc  
      
      
    SOC[0] = DPM[0] + RPM[0] + BIO[0] + HUM[0] + IOM[0]     
      
    Total_Ract = DPM_Ract_new + RPM_Ract_new + BIO_Ract_new + HUM_Ract_new + IOM_Ract
      

# calculate rage of each pool.      
    if(DPM[0]<=zero): 
        DPM_Rage[0] = zero
    else:
        DPM_Rage[0] = ( np.log(DPM[0]/DPM_Ract_new) ) / conr

    
    if(RPM[0]<=zero):
        RPM_Rage[0] = zero
    else:
        RPM_Rage[0] = (np.log(RPM/RPM_Ract_new) ) / conr 
    
    if(BIO[0]<=zero):
        BIO_Rage[0] = zero
    else:
        BIO_Rage[0] = ( np.log(BIO/BIO_Ract_new) ) / conr
    
    
    if(HUM[0]<=zero):
        HUM_Rage[0] = zero
    else:
        HUM_Rage[0] = ( np.log(HUM/HUM_Ract_new) ) / conr
    
    
    if(SOC[0]<=zero):
        Total_Rage[0] = zero
    else:
        Total_Rage[0] = ( np.log(SOC[0]/Total_Ract) ) / conr    
    
    return

def RothC(timeFact, DPM,RPM,BIO,HUM,IOM, SOC, DPM_Rage, RPM_Rage, BIO_Rage, HUM_Rage, \
               Total_Rage, modernC, clay, depth,TEMP,RAIN,PEVAP,PC,DPM_RPM, \
                C_Inp, FYM_Inp, SWC): 
            
     
# Calculate RMFs     
    RM_TMP = RMF_Tmp(TEMP)
    RM_Moist = RMF_Moist(RAIN, PEVAP, clay, depth, PC, SWC)
    RM_PC = RMF_PC(PC)

# Combine RMF's into one.      
    RateM = RM_TMP*RM_Moist*RM_PC
      
    decomp(timeFact, DPM,RPM,BIO,HUM, IOM, SOC, DPM_Rage, RPM_Rage, BIO_Rage, HUM_Rage, Total_Rage, modernC, RateM, clay, C_Inp, FYM_Inp, DPM_RPM)
    
   
    return


######################################################################################################
# program RothC_Python
import os
print(os.getcwd())
os.chdir(["INPUT DIRECTORY PATH") # Change to path of RothC_input.dat
print(os.getcwd())
   
# set initial pool values   
DPM = [0.0]
RPM = [0.0]
BIO = [0.0]
HUM = [0.0]
SOC = [0.0]

DPM_Rage = [0.0]
RPM_Rage = [0.0]
BIO_Rage = [0.0]
HUM_Rage = [0.0]
IOM_Rage = [50000.0]  

# set initial soil water content (deficit) 
SWC = [0.0]
TOC1=0.0
   
# read in RothC input data file
df_head = df_head = pd.read_csv('RothC_input.dat', skiprows = 3, header = 0, nrows = 1, index_col=None, delim_whitespace=True) 
clay = df_head.loc[0,"clay"]
depth = df_head.loc[0,"depth"]
IOM = [df_head.loc[0,"iom"]]
nsteps = df_head.loc[0,"nsteps"]
df = pd.read_csv('RothC_input.dat', skiprows = 6, header = 0, index_col=None, delim_whitespace=True)
print (df)
df.columns =['t_year', 't_month', 't_mod', 't_tmp','t_rain','t_evap', 't_C_Inp', 't_FYM_Inp', 't_PC', 't_DPM_RPM']
  
       
# run RothC to equilibrium     
k = -1
j = -1
      
SOC[0] = DPM[0]+RPM[0]+BIO[0]+HUM[0]+IOM[0]

print (j, DPM[0], RPM[0], BIO[0], HUM[0], IOM[0], SOC[0])

timeFact = 12
    
test = 100.0
while (test > 1E-6):
    k = k + 1
    j = j + 1 

    if( k == timeFact):
        k = 0
    
    TEMP = df.t_tmp[k]
    RAIN = df.t_rain[k]
    PEVAP =df.t_evap[k]
    
    PC = df.t_PC[k]
    DPM_RPM = df.t_DPM_RPM[k]
    
    C_Inp = df.t_C_Inp[k]
    FYM_Inp = df.t_FYM_Inp[k]
    
    modernC = df.t_mod[k] / 100.0   

    Total_Rage=[0.0]

    RothC(timeFact, DPM,RPM,BIO,HUM,IOM, SOC, DPM_Rage, RPM_Rage, BIO_Rage, HUM_Rage, Total_Rage, \
        modernC, clay, depth,TEMP,RAIN,PEVAP,PC,DPM_RPM, C_Inp, FYM_Inp, SWC)  
        
    #each a year calculates the difference between previous year and current year (counter =12 monthly model)
    if (np.mod(k+1, timeFact)== 0):
        TOC0 = TOC1
        TOC1 =DPM[0]+RPM[0]+BIO[0]+HUM[0]
        test = abs(TOC1-TOC0)            
         
Total_Delta = (np.exp(-Total_Rage[0]/8035.0) - 1.0) * 1000.0     
      
print( j, DPM[0], RPM[0], BIO[0], HUM[0], IOM[0], SOC[0], Total_Delta)

year_list = [[1, j+1, DPM[0], RPM[0], BIO[0], HUM[0], IOM[0], SOC[0], Total_Delta[0]]]

month_list = []        

for  i in range(timeFact, nsteps):
  
    TEMP = df.t_tmp[i]
    RAIN = df.t_rain[i]
    PEVAP =df.t_evap[i]
    
    PC = df.t_PC[i]
    DPM_RPM = df.t_DPM_RPM[i]
    
    C_Inp = df.t_C_Inp[i]
    FYM_Inp = df.t_FYM_Inp[i]
    
    modernC = df.t_mod[i] / 100.0
    
    RothC(timeFact, DPM,RPM,BIO,HUM,IOM, SOC, DPM_Rage, RPM_Rage, BIO_Rage, HUM_Rage, Total_Rage, \
          modernC, clay, depth,TEMP,RAIN,PEVAP,PC,DPM_RPM, C_Inp, FYM_Inp, SWC)  
         
    Total_Delta = (np.exp(-Total_Rage[0]/8035.0) - 1.0) * 1000.0
    
    print(C_Inp, FYM_Inp, TEMP, RAIN, PEVAP, SWC[0],  PC,  DPM[0],RPM[0],BIO[0],HUM[0], IOM[0], SOC[0])
    
    month_list.insert(i-timeFact, [df.loc[i,"t_year"],df.loc[i,"t_month"], DPM[0],RPM[0],BIO[0],HUM[0], IOM[0], SOC[0], Total_Delta[0]])
        
    if(df.t_month[i] == timeFact):
        timeFact_index = int(i/timeFact)   
        year_list.insert(timeFact_index, [df.loc[i,"t_year"],df.loc[i,"t_month"], DPM[0],RPM[0],BIO[0],HUM[0], IOM[0], SOC[0], Total_Delta[0]])
        print( i, DPM, RPM, BIO, HUM, IOM, SOC, Total_Delta)

output_years = pd.DataFrame(year_list, columns=["Year","Month","DPM_t_C_ha","RPM_t_C_ha","BIO_t_C_ha","HUM_t_C_ha","IOM_t_C_ha","SOC_t_C_ha","deltaC"])     
output_months = pd.DataFrame(month_list, columns=["Year","Month","DPM_t_C_ha","RPM_t_C_ha","BIO_t_C_ha","HUM_t_C_ha","IOM_t_C_ha","SOC_t_C_ha","deltaC"])

output_years.to_csv("year_results.csv", index = False)
output_months.to_csv("month_results.csv", index = False)

