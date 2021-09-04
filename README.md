# Silicon cluster

This programm tries to find the energedical lowest constellation of a silicon cluster.
The silicon atoms interact with a Bazant force field and for movement a velocity verlet algorithm is used. The optimiztion is dome with the methode of of simulated annealing and furthermore every 100 steps a steepest descent optimiztion with energy feedback is applied.

The project consists of several parts. :

- A initial constellation of silicon atoms is created on a grid
- Using the method of simulted annealing the system envolves and is cooled down step by step. 
- Each 100 steps the the energetically best constellation is searched using a steepest descend optimization method with energy feedback.
- The ten optimal constellations are stored and can be visualized with V_Sim

## To Execute the programm

To execute the .f90 programm a frotran 90 compiler is nedded that creates an .exe file or similar.

Furthermore the file "minimal_energies_silicon.exe" can be executed directly and the simultion is performed.

## Visualization

To visualize the final configuration files for showing with V_Sim are generated.

# Contact

Philipp Bronner: [philippbronner-at-t-online.de](http://philippbronner-at-t-online.de)

Project Link: 

