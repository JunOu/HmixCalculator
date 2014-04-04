HmixCalculator
==============

Heat of mixing calculator using Miedema's scheme

Save the code as *.py and *.dat in the same directory 

Then use the command python *.py to run.

Dependencies: numpy and matplotlib

Note: hmix.decideA() sets a to be 0.14 when the element is alkali and 0.0 for the rest as the paper says only for alkali hmix.decideA() is important. 

Also note: check the original paper to make sure the data is right. There might be some errors since the data is captured from the pdf file and the software seems can not tell '3' and 8 or '5' and '6' very precisely.

Data format: | Element | Phi | Nws1/3 | Vm 2/3 | R/P | Transition Metal (T) or Not (N) 

Units are consistant with those used in the paper.
