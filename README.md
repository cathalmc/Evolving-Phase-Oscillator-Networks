# Evolving-Phase-Oscillator-Networks
Train a network of phase oscillators to climb a hill using a basic genetic algorithm. Should work on most machines if you have a new-ish version of the gcc compiler. Launch with "python runSim.py". The .sh file is for running on a cluster 
Install the gcc compiler for windows from: https://gcc.gnu.org/releases.html


# Quick overview
The "neurons" used here are linear phase oscillators. 
At each timestep their "phase" variable x increases by 1: $x_{n+1} = x_{n} + 1$. 
When $x$ reaches its period, i.e. $x_{n} = P$, it spikes, sending a signal to the neurons its connected to, 
then resets back to 0. The signal can be positive or negative and will increases or decrease the phase $x$ of the "neurons" its connected to.
The code for this in lPRC class in the cpp file.  
![image](https://user-images.githubusercontent.com/22299424/218301520-79c825ab-64c3-400c-8523-71df88eab2a0.png)
  
  

These neurons are connected together in random networks. 
The first 5 neurons are considered input neurons, and they receive a signal (their current height on the hill). 
The next 5 neurons force the system to move to the left when they are active, and the next five push the system to the right. 
The remaining neurons interact with themselves and the intputs and ouputs in order to move the system up or down the hill. 
![image](https://user-images.githubusercontent.com/22299424/218301555-8659d5cc-020a-4091-9b72-fcf47af31f67.png)  

The genetic algorithm initialises a population of eg 100 of these neural systems, and tests how good they are at climbing the hill.
The top 10 are kept into the next generation and the rest deleted. 
These 10 networks are mutated slightly and duplicated (reproduction) to give a new population of 100 better systems.
The process is repeated for many generations and eventually the best hill climbers should emerge. 
