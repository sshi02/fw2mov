# fw2mov.py
## Description
fw2mov.py produces figures and stitches them to produce a .mp4 by pulling from LOG.txt and output produced by FUNWAVE-TVD. 
There is a decent amount of flexibility in configuring input and output, see flags below. 
#### Core features:
- Parses from LOG.txt for model parameters, and auto handles output and bathy paths
- Simply place the file into the working directory, and execute via command line.
- Auto handles 1-D and 2-D, using output to automatically set figure and colorbar limits.
- Auto handles binary and ASCII output IO type.

## Prerequisites 
The script depends on matplotlib and built to expect Python 3.4+ environments

## Quick Start
Drop fw2mov.py into your working folder. 
This folder should contain LOG.txt produced by FUNWAVE-TVD.
Run the following:

```
cd ./<workingdirectory>
python3 fw2mov.py -i LOG.txt
```
Observe that the only necessary argument is the LOG.txt path as input file, which is the only command line argument fw2mov.py techincally needs to produce full movies.
- *1-D cases:* <br />
  fw2mov.py produce a 1-D figure for any model with Nglob = 3. To produce a 1-D figure for any model, specify the `--1d` flag when calling fw2mov.py. fw2mov.py will read the row closest to half of Nglob.

- *2-D cases:* <br />
  To produce a 2-D figure for any model, specify the `--2d` flag when calling fw2mov.py.

## Flags
All flags arguments take precedent over parameters parsed from LOG.txt
- `-h`: help/usage
- `-i, --ifile`: input file
- `-o, --odir`: FUNWAVE output/result folder
- `--1d`: toggle 1-D figure output
- `--2d`: toggle 2-D figure output
- `--wid`: figure width
- `--len`: figure length
- `--arrow`: toggle velocity arrows
- `--dpi`: figure dpi
- `--fps`: movie frames per second
- `--xmin`: figure xlim, minimum argument
- `--xmax`: figure xlim, maximum argument
- `--ymin`: figure ylim, minimum argument
- `--ymax`: figure ylim, maximum argument
- `--mint`: minimum time interval evaluated
- `--maxt`: maximum time interval evaluated

Fully Flagged Example<br />
`$ python3 fw2mov.py -i LOG.txt -odir ./exampleoutput --2d --wid 20 --len 6 --fps 10 --xmin 10 --xmax 300 --ymin 10 --ymax 100 --mint 10 --maxt 20`

## Default Behavior/Values
- dpi = 300
- width = 18
- length = 5
- xmin = 0
- xmax = Mglob
- ymin = 0
- ymax = Nglob
- minimum time = 0
- max time = number of eta files

## Notes
- dep.out in output is tremendously useful and best practice, since fw2py checks for it when constructing bathy
