# stable_lds
This simple repository provides MATLAB code to estimate a linear dynamical system from data solving a convex program. The development branch <a href="http://github.com/epfl-lasa/stable_lds/tree/inverse_lds_">inverse_lds_</a> presents a variant to estimate both the attractor and the nonlinear stable dynamics (assuming a linear parameter varying model) from data. This code relies on YAMLIP and the sedumi solver to solve the convex optimization problems that arise in the estimation. 

To run a simple demo of the code first init and update the respective submodules. In the terminal, go to yout stable_lds folder
```
$ cd your_stable_lds_folder
```
then
```
$ git submodule update --init --recursive
```
To run the demo, in the MATLAB command window run
```
>> demo
```
A figure will pop up where you can draw as many trajectories as you want. Once you are done click 'stop recording' and you will see the streamlines of the resulting dynamical system
![Exemplary linear DS](plot/lds.jpg)

The code can estimate the attractor and dynamics of the provided trajectories or constrain the attractor to a state provided a priori.
