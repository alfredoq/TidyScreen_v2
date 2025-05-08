<img title="a title" alt="Alt text" src="./images/tidy_screen_title.png">

---
### **Project Overview** 
---

TidyScreen is a package developed at [MedChemLab](https://unitefa.conicet.unc.edu.ar/linea-investigacion-quimica-medicinal-y-diseno-de-farmacos/) with the aim of providing a structured framework for the **design**, **execution** and **documentation** of a virtual drug screening campaign.

Overall, the package provides functionalities capable of creating and organizing a given chemical space (*ChemSpace module*) intended to be explored through synthetically feasible chemical reactions. In a nutshell, the corresponding working chemical space is efficiently created from commercially available reactants a first stage of the virtual campaign execution. In this way, synthetic pathways are enumerated.

As part of further screening stages, the execution of molecular docking and molecular dynamics studies are also managed within TidyScreen.

A core feature common to TidyScreen philosophy is the use of [SQL databases](https://www.w3schools.com/sql/sql_intro.asp) to store in an organized fashion all the information relevant the campaign progress, including simulation conditions, results and raw data. In this respect, sharing of the whole project in the context of collaborative work and/or reproducing reported studies is straightforward.

---
#### Installation

```bash
$ conda create -n tidyscreen python=3.10 chemicalite
$ conda activate tidyscreen
$ pip install git+https://github.com/alfredoq/TidyScreen_v2
```
&nbsp;

#### Updating the package

We are constantly adding functionalities to TidyScreen and/or applying identified bugfixes. In order the make the available to the package installed in the corresponding conda environment (without loosing the project database :-P ), run the following command from the corresponding TidyScreen environment:

```bash
pip install --upgrade git+https://github.com/alfredoq/TidyScreen_v2
```

**Additional requirements not installed by CONDA**:


- [Meeko development version](https://github.com/forlilab/Meeko@develop ): this packages was developed by [ForliLab](https://forlilab.org/) at Scripps Research Institute, and is used in TidyScreen workflows to prepare .pdbqt file intended for docking assays. Currently, TidyScreen needs a specific development version of Meeko capable of reading .mol2 atomic charges. This development version should be installed in the corresponding TidyScreen enviroment as follows:

```bash
$ pip install git+https://github.com/forlilab/Meeko@develop  
```

- [*AutoDock-GPU*](https://github.com/ccsb-scripps/AutoDock-GPU): TidyScreen has been prepared to work in conjunction with AutoDock-GPU, which has been developed in the [ForliLab](https://forlilab.org/) at Scripps Research Institute. We acknowledge Stefano Forli, Diogo Santos-Martins and Andreas Tillack for the kind feedback during TidyScreen development.
&nbsp;
- [Amber](https://ambermd.org/) MD engine: this software package is required to confer TidyScreen the capability to prepare and document molecular dynamics simulations of docked poses.

---

In order to use TidyScreen, users can access the [documentation/WIP](null) describing the project and specific working examples.
