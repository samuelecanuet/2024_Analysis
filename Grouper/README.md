# Analysis program for WISArD experiment

This code contains some routines to analyse WISArD experimental data from raw root files. It can be use to group data in event between detectors, clean them and merge all the run after a calibration.

## Prerequisites
- **ROOT** 
- **GSL**
- **root2mpl** for plotting some spectrums in python

## Usage
To compile the code : 

```bash
make clean; make all
```

### To group the data in events
You can personnalize time window coincidence between detectors to built group and then make event. 
*Rear* silicon detectors signal is taken as trigger and all the others as followers. Additionaly double strip signal or double rear signal are discarded and not included in the saved tree. Finnaly a cleaning is performed to avoid remove background and strange energy deposited with respect to rear/strip signals.

```bash
Grouper <Run Number>
```

### To gain match runs 
The goal is to match a run with a reference run. For all the silicon detector the reference run is 115, only on 32Ar with the thin catcher.

```bash
Matcher <Run Number>
```
or 

```bash
Matcher all
```

### To calibrate detectors
For silicon detectors 

## TODO
- Gain match SiPMs
- Calibrate SiPMs

## Sequence executing
- **Grouper**           : Input *.root                              /           Output *_grouped.root
- **Matcher**           : Input *_grouped.root                      /           Output matched.root
- **Merger**            : Input *_grouped.root & matched.root       /           Output **_merged.root
- **Calibration**       : Input **_merged.root                      /           Output Calibrated.root
- **SiPM_Matching**     : Input **_merged.root                      /           Output SiPM_Matched.root
- **SiPM_Calibration**  : Input **_merged.root & SiPM_Matched.root  /           Output SiPM_Calibrated.root
- **Analysis**          : Input **_merged.root 
                            & Calibrated.root                       /           Output **_analysed.root
                            & SiPM_Calibrated.root

