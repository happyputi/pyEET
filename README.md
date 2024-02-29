# General Information
- Photosystem I (PSI) network analysis.
- Python code to cumpute excitation energy transfer between chlorophylls in PSI LHC supercomplex.

# General notes
- This repository hosts python code to compute FRET rate and the network effect of the excitation energy transfer of PSI LHC supercomplex.
- The python code `Forster.py` analyzes the Förster resonance energy transfer rate between chlorophylls. 
- The python code `k2t.py` analyzes the Förster resonance energy transfer time between chlorophylls. 
- The python code `GF.py` analyzes the generalized Förster energy transfer rate between pigment clusters. 

# Requirements
Numpy

# Usage
1 Computational analysis of FRET rates between Chlorophylls  
Excuting `Forster.py`  by python can cumpute FRET rates between Chlorophylls:
```bash
python3 Forster.py H_matrix_file spec_tmppath
```
2 Calculate the energy transfer time  
Excuting `k2t.py`  by python can cumpute energy transfer time between Chlorophylls based on the energy transfer rate matrix obtained in the previous step:
```bash
python3 k2t.py k_filename t_filename filelistname
```
3 Computational analysis of generalized Förster energy transfer rate between pigment clusters  
Excuting `GF.py`  by python can cumpute energy transfer rate between pigment clusters:
```bash
python3 FGF.py H_matrixfile Donor_atart Doner_end Accept_start Accept_end spectmp outdir
```

# Contact
Created by [@gaojun](gaojun@mail.hzau.edu.cn) and [@guojianping](guojianping@webmail.hzau.edu.cn) feel free to contact us!

