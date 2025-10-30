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

**Suggested installation procedure:** A simplified installation procedure is to download the [bash installation script](https://github.com/alfredoq/TidyScreen_v2/blob/main/tidyscreen_installation.sh) and execute it as follows:

```bash
$ chmod 777 tidyscreen_installation.sh
$ ./tidyscreen_installation.sh

## The creation con dedicated environment named 'tidyscreen' will be accomplished

$ conda activate tidyscreen

# Start working
```

**Step-by-step installation** : Below you will find instruction on how to manually perform all the actions included in the [bash installation script](https://github.com/alfredoq/TidyScreen_v2/blob/main/tidyscreen_installation.sh) in case any specific parameter is required to be modified:

```bash
$ conda create -n tidyscreen python=3.12 
$ conda activate tidyscreen
$ pip install git+https://github.com/alfredoq/TidyScreen_v2

# Next, we will install several accesory packages required by TidyScreen
$ conda install -c conda-forge ambertools==23.6 espaloma espaloma_charge chemicalite visidata vmd-python

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
- [Autodock4](https://forlilab.org/code/): is required to compute docking grids.
```bash
$ conda install -c bioconda autodock autogrid
```

- [*AutoDock-GPU*](https://github.com/ccsb-scripps/AutoDock-GPU): TidyScreen has been prepared to work in conjunction with AutoDock-GPU, which has been developed in the [ForliLab](https://forlilab.org/) at Scripps Research Institute. We acknowledge Stefano Forli, Diogo Santos-Martins and Andreas Tillack for the kind feedback during TidyScreen development.
&nbsp;

- [Ersilia Models Hub](https://ersilia.gitbook.io/ersilia-book): this package is required if the user intends to use make use of the [Ersilia Open Source Initiative](https://www.ersilia.io/) prediction models as part of chemical space prioritization. To install the hub run the following code within the TidyScreen environment:

```bash
# With the TisyScreen environment activated, clone the source repo and install it:
$ git clone https://github.com/ersilia-os/ersilia.git
$ cd ersilia
$ pip install -e .
$ bentoml # this action will install the BentoML version required by the Ersilia Model Hub
```

- [Docker](https://www.docker.com/products/docker-desktop/): Ersilia Models Hub requires a system-wide installation of Docker (container management system). Installation instructions for different operating systems can be found [here](https://docs.docker.com/desktop/setup/install/linux/), while instructions specific for Ubuntu can be found [here](https://docs.docker.com/desktop/setup/install/linux/ubuntu/)

- [Redis server](https://redis.io/open-source/):

```bash
$ conda install -c conda-forge redis-server
```

-[AutodockTools] ligand preparation scripts: Some actions (i.e. parameterization of ligands for dockings) uses custom scripts that needs to be installed into a separate conda environment named 'adt', and that uses an older Python version (2.7). To create that environment and the corresponding tools execute the following commands:

```bash
$ conda create -n adt python=2.7 
$ conda activate adt
$ conda install -c insilichem autodocktools-prepare
```

In order to use TidyScreen, users can access the [documentation](https://alfredoq.github.io/TidyScreen_v2_docs_new/) describing the project and specific working examples.
