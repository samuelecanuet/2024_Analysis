# MCP Reconstruction Script

This Python script is designed to read and indicate corners of a grid with the mouse, in order to generate a guess txt file to fit the data.

## Prerequisites

- **GSL** (https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz)
- **fasterac** (https://faster.in2p3.fr/index.php/download/category/2-software?download=20:fasterac-2-17-0-tar-gz)
- **ROOT**
- run commands :

```bash
./configure
```

```bash
source ~/.cshrc
```
## Usage

To run the script, use the following command:

```bash
python MCP_gui.py -f <input_filename> [-y]
```

- `-f` or `--filename`: Specify the input file (either .root or .fast format). 
- `-y` OR `--year`: Specify the file year


## Version 
### v1.0 : 
- GUI for calibration and reconstruction
- Manual calibration with corner point
- Polynomial fit
- Gaussian 1D
