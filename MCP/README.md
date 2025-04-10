////// README FOR PLOTTING ////////

A) Use the GUI to generate the guess corner txt file
B) Fitting:

- analysis explain in google slides
- use Calibration.xx (writting true)
- running with Calibration 1 2 3 4 5 6 7 8 9 (par of fit if not fixed and thread where chi2 written)
    ---> All the analysis is ON 
    1) fitting grid with grid convoluted from guess_centers.txt to out_centers.txt
    2) Fitting 2D polynomial from x, y to X and Y (calibration of the grid)
    3) Using fitted fucntion on measurement data and fit [grid x gauss] * Res
- plots in Calibration.root



- Details:
-- calibration of the grid at 4T with MCP_008_4T
-- measurement of the grid at 4T with MCP_010_4T
