A program to generate cyclic voltammograms under different circumstances. Examples of generated CVs are shown in examples/

The program takes an input file and command line options. An example input file is shown below;

```CONTENTS
Ox = 100
Red = 0
S = 0
P = 0
ECHEM
peake:7
scanrate:0.1
cells:50
output:1000
diffusion:0
REACTIONS
S + Red > P + Ox (1)
```

In CONTENTS all chemical species should be defined. Ox and Red are required for the electrochemistry to work, however for a purely chemical simulation these can be set to 0. 

In ECHEM there are 5 options. Peake defines the maximum variation in the potential; the potential will scan from -peake to peake. Scanrate determines how quickly the model moves from one to the other (this does also depend on the timestep as detailed below). cells controls how many cells are produced; These have spatial diffusion which allows for segregation of bulk and electrode. Changing line 258 in cell.py can allow you to print out the concentrations for all cells; analysing this output should allow concentration gradients to be plotted. 

output determines how frequently output is written to the file. For small time steps larger values may be useful to prevent the output files being massive. diffusion determines how freely diffusion occurs. This is based on a very simple model where the new concentration in position x is given by (diffusion*(x-1) + (1-2*diffusion)x + diffusion*(x+1)) where x-1 and x+1 refer to the neighbouring cells.

REACTIONS defines all non electrochemical reactions occurring. These are defined as;

`REACTANTS > PRODUCTS (RATE CONSTANT)`

In the example above there is no substrate and no diffusion; this would therefore be expected to give a centrally symmetric cyclic voltammogram. One thing I would like to change is to have independent diffusion coefficients for different species, which would allow for things like surface mounted species to be used too.

The file is run by

`python file FILENAME LENGTH TIMESTEP OUTPUT`

So if I wanted to run a file called CELL1 I might do this by;

`python file CELL1 2 0.001 CELL1.out `

The length and timestep may need to be tuned to prevent overflows.
