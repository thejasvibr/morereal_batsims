# More Real BatSims
A module that generate somewhat more realistic bat echolocation simulations. 

IMPORTANT : this module only generates SINGLE bat flight & call paths

Absolutely required inputs are:

* `-seed` : the seed number is the only thing that is required, aside from the other arguments that can be changed (Type `python batsimulation.py -h` for more info)

## An example run

```
python batsimulation.py -seed 78464
```


For now there are a bunch of hard-coded things:

* Array geometry : an 8 channel array with mics on/very close to the wall
* Room geometry : it's a shoe-box room. The dimensions can be changed
* Max/min speed : set to be within 3-6 m/s
* sampling rate: set to 192 kHz
* max reflection order: set to 1st order (source-wall-mic)

Other things that are randomly chosen - and dependent on the starting seed:

* trajectory : created by choosing a bunch of random points within the x,y,z limits, creating a spline and checking if the speed is within the min/max limits
* call duration: chosen between 5-7 ms
* min/max call frequency : 15-22 & 88-92 kHz
* call type : one of logarithmic, linear or hyperbolic



Baseline simulations
--------------------
The baseline simulations were created with the following seed numbers: 1234, 3456, 5678 
To replicate the runs enter the following code:


```
python batsimulation.py -seed <seednumberhere>
```

Package requirements
--------------------
pyroomacoustics
matplotlib 
numpy
scipy
pandas
