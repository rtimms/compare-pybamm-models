# Compare models in PyBaMM

This repository shows an example of how to compare models for speed/accuracy in PyBaMM. It contains two parameter sets from the literature:

1. Parameters for an LG M50 cell taken from the paper 

    > Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W. Dhammika Widanage, and Emma Kendrick. ["Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models."](https://iopscience.iop.org/article/10.1149/1945-7111/ab9050) Journal of the Electrochemical Society 167 (2020): 080534

    with additional thermal parameters from
    
    > Ferran Brosa Planella, Muhammad Sheikh, and W. Dhammika Widanage, ["Systematic derivation and validation of a reduced thermal-electrochemical model for lithium-ion batteries using asymptotic methods."](https://www.sciencedirect.com/science/article/pii/S0013468621008148?via%3Dihub) Electrochimica Acta Volume 388, 2021, 138524
    
    The parameter set is stored in `Chen2020_params.py` and provides functional forms for parameters that depend on state (e.g., open circuit potential as a function of concentration) and scalar values for the remaining parameters. The associated measured data for the functions are also provided as `.csv` files in the `data` folder with the same name as the function.

2. Parameters for a Kokam SLPB 75106100 cell, from the papers

    > Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of a lithium-ion battery I. determination of parameters." Journal of the Electrochemical Society 162.9 (2015): A1836-A1848.

    >Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of a lithium-ion battery II. Model validation." Journal of The Electrochemical Society 162.9 (2015): A1849-A1857.

    The thermal material properties are for a 5 Ah power pouch cell by Kokam. The data are extracted from

    > Zhao, Y., et al. "Modeling the effects of thermal gradients induced by tab and surface cooling on lithium ion cell performance."" Journal of The Electrochemical Society, 165.13 (2018): A3169-A3178.

    The parameter set is stored in `Ecker2015_params.py` and provides functional forms for parameters that depend on state (e.g., open circuit potential as a function of concentration) and scalar values for the remaining parameters.    
    

## 🚀 Installation
In order to run the models and load in any data you will need to install `pybamm`. The notebooks have been tested on PyBaMM Version 22.2. To install the required python packages on Linux/Mac OS use the following terminal commands:

1. Clone the repository
```bash
https://github.com/rtimms/compare-pybamm-models
```
2. Change into the `compare-pybamm-models` directory 
```bash
cd compare-pybamm-models
```
3. Create a virtual environment (optional)
```bash
virtualenv env
```
4. Activate the virtual environment (optional)
```bash
source env/bin/activate
```
5. Install the required packages
```bash 
pip install -r requirements.txt
```

PyBaMM is available on GNU/Linux, MacOS and Windows. For more detailed instructions on how to install PyBaMM, see [the PyBaMM documentation](https://pybamm.readthedocs.io/en/latest/install/GNU-linux.html#user-install).
