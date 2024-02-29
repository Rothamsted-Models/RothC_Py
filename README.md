# The Rothamsted carbon model (RothC)

## Purpose

RothC models the turnover of organic carbon in non-waterlogged top-soil.  It accounts for the effects of soil texture, temperature, moisture content and plant cover on the turnover process. It uses a monthly time step to calculate total organic carbon (t ha<sup>-1</sup>), microbial biomass carbon (t ha<sup>-1</sup>) and Δ<sup>14</sup>C (from which the equivalent radiocarbon age of the soil can be calculated). 

## Development history

The first version of RothC created by David Jenkinson and James Rayner in 1977 (Jenkinson and Rayner, 1977).

In 1987 an updated version was published, see Jenkinson et al. (1987).  This version included the prediction of the radiocarbon age of the soil, the pools POM (physically stabilized organic matter) and COM (chemically stabilized organic matter) were replaced with Hum (humified organic matter) and IOM (inert organic matter), and the microbial biomass pool was split into BioA (autochthonous biomass) and BioZ (zymogenous biomass).  

**In 1990, the two biomass pools were combined into a single pool (Jenkinson, 1990) this version is the standard version of the model, that this code refers to.**

Other published developments of the model include:

Farina et al. (2013) modified the soil water dynamics for semi-arid regions.

Giongo et al. (2020) created a daily version and modified the soil water dynamics, for Caatinga shrublands, in the semiarid region, North-East Brazil.

 
## Description of files included

### RothC_description.pdf
This file contains the description of the model.


### RothC_Py.py
This file contains the RothC code in Python language. Details of the inputs required, pools modelled, and units are in the code.


### RothC_input.dat  
This file contains input variables for the model.  

At the start of the file values for **clay** (%), **soil depth** (cm), **inert organic matter** (IOM, t C ha<sup>-1</sup>) and **number of steps** (nsteps) are recorded.  
Following that there is a table which records monthly data on **year**, **month**, **percentage of modern carbon**  (%), **mean air temperature** (Tmp, °C), **total monthly rainfall** (Rain, mm), **total monthly open-pan evaporation** (Evap, mm), **all carbon input entering the soil** (from plants, roots, root exudates) (C_inp, t C ha<sup>-1</sup>), **carbon input from farmyard manure** (FYM, t C ha<sup>-1</sup>), **plant cover** (PC, 0 for no plants e.g. bare or post-harvest, 1 for plants e.g. crop or grass), and the **DPM/RPM ratio** (DPM_RPM) of the carbon inputs from plants.

### year_results.csv
This file contains the yearly values of the SOC (both the pools and Total) and the delta 14-carbon.

The pools are:  
**Year**  
**Month** 	- Always December for the yearly output  
**DPM** 	- Decomposable plant material (t C ha<sup>-1</sup>)  
**RPM** 	- Resistant plant material (t C ha<sup>-1</sup>)  
**BIO** 	- Microbial biomass (t C ha<sup>-1</sup>)  
**HUM**	- Humified organic matter (t C ha<sup>-1</sup>)  
**IOM** 	- Inert organic matter (t C ha<sup>-1</sup>)  
**SOC**	- Total soil organic carbon (t C ha<sup>-1</sup>)  
**deltaC** 	- delta <sup>14</sup>C (‰)  


The total organic carbon (soil organic carbon) is equal to the sum of the 5 pools. 

TOC or SOC = DRM + RPM + BIO + HUM + IOM 
     
### month_results.csv
This file contains the monthly inputs, rate modifying factors, SOC pools.

**Year**  
**Month**  
**DPM_t_C_ha**		- Decomposable plant material (t C ha<sup>-1</sup>)  
**RPM_t_C_ha**		- Resistant plant material (t C ha<sup>-1</sup>)  
**BIO_t_C_ha**		- Microbial biomass (t C ha<sup>-1</sup>)  
**HUM_t_C_ha**		- Humified organic matter (t C ha<sup>-1</sup>)  
**IOM_t_C_ha**		- Inert organic matter (t C ha<sup>-1</sup>)  
**SOC_t_C_ha**		- Total soil organic carbon (t C ha<sup>-1</sup>)  

## Requirements
The code was written in Python 3.9.7.

## Installation/set-up

A directory path will need to be provided as indicated in the code ([“INPUT DIRECTORY PATH”]), to read in RothC_input.dat. 

**Example of how to run the model**  
The file RothC_input.dat contains all the inputs data needed to run the model. The month results (month_results.csv) and year results (year_results.csv) files correspond to this input file as an example. 
The model is normally run to equilibrium using average temperature, rainfall, open pan evaporation, an average carbon input to the soil, the equilibrium run is to initialise the soil carbon pools. Once the soil carbon pools have been initialised, the model is run for the period of interest. The met data (temperature, rainfall and evaporation) can be average or actual weather data. The carbon input to the soil can be: 1) adjusted so the modelled output matches the measured data, or 2) can be estimated from yield data (Bolinder et al., 2007), or NPP data.  


## References

Bolinder MA, Janzen HH, Gregorich EG, Angers DA, VandenBygaart AJ. An approach for estimating net primary productivity and annual carbon inputs to soil for common agricultural crops in Canada. Agriculture, Ecosystems & Environment 2007; 118: 29-42.  
Farina R, Coleman K, Whitmore AP. Modification of the RothC model for simulations of soil organic C dynamics in dryland regions. Geoderma 2013; 200: 18-30.  
Giongo V, Coleman K, Santana MD, Salviano AM, Olszveski N, Silva DJ, et al. Optimizing multifunctional agroecosystems in irrigated dryland agriculture to restore soil carbon - Experiments and modelling. Science of the Total Environment 2020; 725.  
Jenkinson DS. The Turnover of Organic-Carbon and Nitrogen in Soil. Philosophical Transactions of the Royal Society of London, Series B: Biological Sciences 1990; 329: 361-368.  
Jenkinson DS, Hart PBS, Rayner JH, Parry LC. Modelling the turnover of organic matter in long-term experiments at Rothamsted. INTECOL Bulletin 1987; 15: 1-8.  
Jenkinson DS, Rayner JH. Turnover of soil organic matter in some of the Rothamsted classical experiments. Soil Science 1977; 123: 298-305.  

