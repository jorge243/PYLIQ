# PYLIQ
State of the art free soil liquefaction analysis software from SPT data
## Installation
```bash
$ pip install PYLIQ
```
## Usage example:
```python
from PYLIQ import analyze

input_file = 'PYLIQ Input.xlsx'

analyze(input_file)
```
 ### Features:
 - Fast and easy data input from excel file
 - Batch processing from multiple SPT borings at once
 - Automated generation of results data in excel format
 - Automated generation of results graphics in pdf or png format
 ### Analyses options:
 - Deterministic Standard Penetration Test (SPT) liquefaction triggering analysis according to state of the art method proposed by Idriss & Boulanger (2014)
 - Estimation of liquefaction effects / consequences:
    - Liquefaction potential index (LPI)
        - Iwasaki et al. (1978)
    - Liquefaction severity number (LSN)
        - van Ballegooy et al. (2014)
    - Liquefaction-induced damage
        - Ishihara (1985)
    - Lateral spreading / displacements
        - Idriss & Boulanger (2008)
        - Youd et al. (2002)
    - Post-liquefaction vertical reconsolidation settlements
        - Idriss & Boulanger (2008) *based on Ishihara & Yoshimine (1992)*
    - Post-liquefaction residual strength
        - Olson & Stark (2002)
        - Idriss & Boulanger (2008)
        - Weber (2015)
