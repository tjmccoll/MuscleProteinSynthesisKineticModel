# Kinetic modeling of leucine-mediated signaling and protein metabolism in human skeletal muscle

## Overview

The muscle protein synthesis kinetic model allows for the simulation of leucine-mediated mTOR signaling via insulin- and leucine-independent pathways and protein synthesis in skeletal muscle of young adults. The model can be used to simulate leucine feeding of different doses and timing to visualize the downstream signaling cascade. 

This tool is being developed by the Clarke Laboratory for Quantitative Exercise Biology. Further information on the muscle protein synthesis kinetic model can be found in the following study.

> Insert link to BioRxiv or published paper
<img src="https://media.github.sfu.ca/user/1053/files/4532768b-11d6-48e7-bad4-34be4f65177e">

## Getting Started
The 'McColl_2023_MuscleProteinSynthesisKineticModel_230221' folder contains all required files (scripts, experimental data, input data) to run the model and replicate the manuscript figures. 

### Installation
1. Download the 

### Running the model

To simulate the model:
1. Open the ‘McColl_2023_model execution_230221.mlx’ live script in the ‘Code’ folder in MATLAB 2022a or newer. 
2. Update the ‘FolderPath’ on line 4 to the file path that contains the locally saved ‘McColl_2023_Muscle Protein Synthesis Kinetic Model_230221’ folder. I.e., update the fileparts(‘…’) function.
3. Run each section of the script.
4. Each manuscript figure provides the prompt: ‘Save figure? Y/N’. Enter 'Y' in the command window to save the figure as an .eps file. The figure will then be saved as an .eps file in the ‘Output Plots’ folder. Enter ’N’ in the command window to only create figure within the Matlab live script output. 
