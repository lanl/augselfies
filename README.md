## augSELFIES
Augmented implementation and modified versions of Self-Referencing Embedded Strings (SELFIES).

The core distinction here is that while the original SELFIES implementation focuses on interconverting 
between SMILES and SELFIES strings with molecular graphs predominantly used as an afterthought, here
graphs and operations on them are at the forefront. Secondarily, native SELFIES re-uses chemical symbol tokens
for indexing, whereas this implementation supports alterative indexing schemes. 

For the original implmentation of SELFIES and references, see https://github.com/aspuru-guzik-group/selfies.


===

# Installation 

First install pixi onto your machine as well as a mirror of group-selfies (https://github.com/aspuru-guzik-group/group-selfies). Next, while in a directory for augselfies, run

```bash 
pixi install 
```
to install default required dependencies and
```bash 
pixi shell 
```
from the root directory of this project to install the dependencies and to activate the envioronment.

PyPI/pip Instructions:

It is *not* recommended to install from pip, which loses some of the precision that pixi provides for general dependency management. However, for
ease of installation this project may be installed from source via 

```bash
pip install .["pip"]
```
or from PyPI via 

```bash
pip install augselfies["pip"]
```

This is provided as a convenience but may cause unintended behavior. 

12.02.2025: GroupSELFIES must be separately sourced at this time.

# Usage 

augSELFIES should be used as a library, predominantly for data processing. 

It contains methods for converting SELFIES to numeric SELFIES (numSELFIES), implementations of 
which are in `augselfies.numeralization`. 

augSELFIES also contains processes for data augmentation, creating multiple equivalent SELFIES/numSELFIES for the some underlying
molcular graph. See `augselfies.augmentation` for details and implmentation. 

# Testing 

This repository is desgined to be tested via pytest. Run 
```bash
pixi shell --environment test
pytest --cov src 
```
to run all unit tests and determine current code coverage. 
# API Documentation

See the HTML files in `/docs`. To regenerate documentation, run 
```bash
pixi shell --environment docs 
sphinx-build -M html docs docs/_build 
```

# Copyright Notice 

This project is MIT-licensed (see `LICENSE.txt`). 


© 2025. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos

National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.

Department of Energy/National Nuclear Security Administration. All rights in the program are

reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear

Security Administration. The Government is granted for itself and others acting on its behalf a

nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare

derivative works, distribute copies to the public, perform publicly and display publicly, and to permit

others to do so.


This project has been approved for open-source release under number O#: O4990.
