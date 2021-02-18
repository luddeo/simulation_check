# Example file

A simple model containing three species (Aa, Bb and Cc) that diffuse in a cell geometry. the species Aa has a diffusion coefficient of 0.05 um^2/s, Bb of 0.5 um^2/s and Cc of 5 um^2/s. 

A particle of Aa type turns into a Bb type with a rate of 1 1/s (giving a dwell timer of 1 s), Bb type turn into Cc type with a rate of 0.2 1/s (giving a dwell time of 5 s) and Cc type turn into Aa type with a rate of 2 1/s (giving a dwell time of 0.5 s).

In the beginning there are 1000 particles of the type Bb and the simulation run for 10.5 seconds (simulated time) to equilibrate the system, and then continue to run until 40 seconds (simulated time).

The geometry of the cell is divided into small cubes with sides of a length 0.01 um. The scale used for the coordinates are in cubes.

The model is executed with MesoRD:
> mesord "simple_model.xml" -i 1 -I 50 -c 0.02 -C -1 -E -p -g -t 40 -q 0.010 um -K -x 10.5 -w 0.0057

The previous script calculated the diffusion coefficient (D), occupancy (pOcc) and dwell time (DT) of the different species. 


|    | Aa | Bb | Cc |
| :--- | :---: | :---: | :---: |
| D: | 0.04889986  +/-  0.04191294 | 0.4702457  +/-  0.3915139 | 4.045476  +/-  3.401735 |
| pOcc: | 0.1532398  +/-  0.0105523 | 0.7705715  +/-  0.01312642 | 0.07618876  +/-  0.00770415 |
| DT: | 0.9556315  +/-  0.9279875 | 4.011812  +/-  3.804945 | 0.4871705  +/-  0.4905858 |
