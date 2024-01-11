# Kinetic modeling of leucine-mediated signaling and protein metabolism in human skeletal muscle

## Overview

The muscle protein synthesis kinetic model allows for the simulation of leucine-mediated mTOR signaling via insulin- and leucine-independent pathways and protein synthesis in skeletal muscle of young adults. The model can be used to simulate leucine feeding of different doses and timing to visualize the downstream signaling cascade. 

This tool is being developed by the [Clarke Laboratory for Quantitative Exercise Biology](https://www.sfu.ca/clarkelab-bpk.html). Further information on the muscle protein synthesis kinetic model can be found in the following study.

> McColl, T. J. & Clarke, D. C. Kinetic modeling of leucine-mediated signaling and protein metabolism in human skeletal muscle. iScience 108634 (2023) doi:10.1016/j.isci.2023.108634.
<picture>
  <img src="https://media.github.sfu.ca/user/1053/files/4532768b-11d6-48e7-bad4-34be4f65177e">
</picture>

## Getting Started
The 'McColl_2023_MuscleProteinSynthesisKineticModel_231024' folder contains all required files (scripts, experimental data, input data) to run the model and replicate the manuscript figures. 

### Installation
Download the package to a local folder (e.g., '~/MuscleProteinSynthesisKineticModel/') by extracting the ZIP file or by running the following terminal command: 
```
git clone https://github.com/tjmccoll/MuscleProteinSynthesisKineticModel.git
```

### Running the model

To simulate the model:
1. Open the ‘McColl_2023_model execution_230922.mlx’ live script from the ‘Code’ folder in MATLAB 2022a or newer. 
2. Update the ‘FolderPath’ variable on line 4 to the file path that contains the locally saved ‘McColl_2023_Muscle Protein Synthesis Kinetic Model_231024’ folder. I.e., update the fileparts(‘…’) function.
3. Run each section of the script.
4. Each manuscript figure provides the prompt: ‘Save figure? Y/N’. Enter 'Y' in the command window to save the figure as an .eps file. The figure will then be saved as an .eps file in the ‘Output Plots’ folder. Enter ’N’ in the command window to only create figure within the Matlab live script output. 

### File list
#### 'Code' folder:
* ‘McColl_2023_model execution_230922.mlx’: This script contains the code to simulate the model and create all the plots included in the manuscript and supplementary information. 
#### 'Experimental Data' folder:
* ‘230221_experimental data.xlsx’: This excel file contains all experimental data collected for model calibration and validation.
* ‘200832_Biolo, 3-pool parameters.xlsx’: 3-pool parameter values at baseline. These experimental values were used to calibrate the baseline k-values that control the 3-pool parameters. The comparison between the simulated and experimental baseline 3-pool parameter values are calculated at lines 42-43 in the ‘McColl_2023_model execution_230922.mlx’ script.
* '230308_Unit conversions.xslx': An excel file that provides the calculations for the conversion of units from experimental data sets to the input data in the model.
#### 'Functions' folder:
* Contains all required functions to run the model execution script. Information pertaining to each function is contained within the function script.
#### 'Input Data' folder:
* ’230922_Initial Values, IOM.xlsx’: Contains the initial values for each species that is inputted into the ODE function. 
* ’231023_K-values, IOM.xlsx’: Contains the kinetic rate parameters that are inputted into the ODE function.
#### 'Meta-analysis' folder:
* Contains the R code and experimental data needed replicate the p-Akt and p-p70S6K meta-analysis and spline regressions.
#### 'Optimizer' folder:
* 	•	Contains the ‘runOptimizer_230923’ code to run the parameter optimizer and the ‘objectiveFunction_230923’ function that the optimizer seeks to minimize.
#### 'Output Plots' folder:
* Empty folder
* The model execution script will create a folder corresponding to the date and time that line 9 is run where saved model figures are stored.

## Contact
tmccoll@sfu.ca
