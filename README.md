# Directory structure & content
___

This application is to visualize HLA and epitope on the in 3D. The application, can do the visualisation per transplant number in PROCARE 1 or per given set of epitopes. In the latter case, the minimum number of HLA (which are of highest expression) are depicted that could accommodate all the given epitopes. Filtering can be done via

- Donor type
- HLA class: HLA should be passed in as string
- HLA Molecule
- Transplant Number
- ElliPro Score
- By Epitopes: Epitopes should all be based in as a comma separated string like '45GV','45EV', '46VY' or  '45GV,45EV, 46VY',
...


# Directories:
---

- app: main application content with the front-end and all the relevant files for the backend

- Notebook: Contains data preparation for the consumed data files in data directory and some

- data: This repository includes all the required data by the application. Please MAKE SURE these are not tracked by cloud based data or code version controls
    - HLA pdb files
    - HLA expression data base
    - [Epitope registry data base](https://www.epregistry.com.br/index/databases/database/ABC): This also includes some feature absent in the registry for example epitope  distance to cell membrane. Link to
    - Log: The log file created during the running of application
    - DESA and Transplant data

# Running the application:
---

The application can be launched via running the app.py file.

# Environments for Productionising:
---

- Pipfile: for making local virtual environments for installing python packages. To produce the virtual environment, below steps are required:

    * Install the pipenv package `pip install --user pipenv`
    * Create the virtual environment `pipenv install --python 3.8`. (Pipfile packages will automatically install)
    * Activate the virtual environment `pipenv shell`

- There are some docker-compose and Dockerfile available to productionize the application


# Adjustments to Epitope visualisation:
---

In the `min_hlavsep` and `_get_min_hlavsep` method, HLA frequencies are included to only visualize the epitopes belonging to [high frequency HLAs](http://www.allelefrequencies.net). Thus, pay a careful attention to most_freq_hla flag in applying the method.
