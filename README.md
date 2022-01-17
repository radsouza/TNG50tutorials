# TNG50tutorials


Setup:

Download the illustris_python library from https://github.com/illustristng/

Essential files to download from the TNG50-1 respository.
- Sublink merger trees (~613.5 GB across 11 files)
- Snapshot-99 subhalo and groupcatalog (~5.9 GB), also offsets (~1 GB).
Organize these files according to the propoer data structure

Optional files needed to choose MW-mass galaxies.
- (Selection on rotation) stellar Circularities, Angular Momenta, and Axis Ratios (100 Mb in 1 file)
- (Selection on accreted stellar fraction) Stellar Assembly Catalog (~34 GB)
- (Section on aperture masses) To choose MW-mass galaxies Aperture masses (~50 MB)


There are four main notebooks in this directory:
1. ChoosingMWhaloes.ipynb: A notebook to choose MW-mass galaxies based on mass, isolation and circularity. Outputs a file MWM31SubfindID.npy which stores the data.
2. Get_merger_tree-TNG-50.ipynb: A notebook to find the massive progenitors and massive satellites of the chosen MW-mass galaxies. Outputs a file ans_TNG50.npy which store information about the progenitors/satellites of the chosen galaxies.
3. AccretedStars-Single.ipynb: Taking a MW-mass galaxy, orienting it and dividing into insitu and accreted stars
4. InfallingSubhaloes-Single.ipyb: Find all the subhaloes that fell into the galaxy. Demonstrate that there is a scattering of binding energy in the surving subhaoes, when there is a massive accretion.




### Dividing stellar particles into insitu/accreted

- Separting stars into insitu and accreted is done with the use of a file which details the birth snapshot, subhalo and group of every stellar particle found in Snapshot99: origins_99.hdf5 (11 Gb).

- One calculates the tags for the stellar particles using the code findtags.py. Given a subhalo at snapshot 99, it searches the merger trees and divides the particles into three categories: i) those born on the main progenitor branch (mpb), ii) those born within the main virial halo (but not on the mpb), iii) those born outside the virial halo.

- For now, to speed up calculations for a large number of galaxies, the code uses a helper data file (data_99.hdf5 ~9.8 GB), which is derived from the snaphot data, and helps find the particleIDs for a given subhalo or a given group.

- You need to get origins_99.hdf5 and data_99.hdf5 from me.


