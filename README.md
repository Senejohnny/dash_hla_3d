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
    - Epitope registry data base: This also includes some feature absent in the registry for example epitope  distance to cell membrane.
    - Log: The log file created during the running of application
    - DESA and Transplant data

# Environments for Productionising:
---

- Pipfile: for making local virtual environments for installing python packages

- There are some docker-compose and Dockerfile available to productionize the application
