# fw2mov.py
## Description
fw2mov.py produces figures by pulling from FUNWAVE-TVD LOG.txt and output. 
The advantage to pulling from LOG.txt is that most data needed to produce figures is neatly standardized by FUNWAVE-TVD output.
There is a decent amount of flexibility in configuring input and output, see flags below. 
Simply place the file into the working directory, and execute via command line.

## Prerequisites 
The script depends on matplotlib.

## Quick Start
Drop fw2mov.py into your working folder. 
This folder should contain LOG.txt produced by FUNWAVE-TVD.
Run the following:

```
cd ./<workingdirectory>
python3 fw2mov.py -i LOG.txt
```

- *1-D cases:* <br />
  fw2mov.py produce a 1-D figure for any model with Nglob = 3. To produce a 1-D figure for any model, specify the `--1d` flag when calling fw2mov.py. fw2mov.py will read the row closest to half of Nglob.

- *2-D cases:* <br />
  To produce a 2-D figure for any model, specify the `--2d` flag when calling fw2mov.py.

## Flags
All flags arguments take precedent over parameters parsed from LOG.txt
- `-h`: help/usage
- `-i, --ifile`: input file
- `-o, --odir`: FUNWAVE output/result folder
- `--1d`: 1-D figure output
- `--2d`: 2-D figure output
- `--wid`: figure width
- `--len`: figure length
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

