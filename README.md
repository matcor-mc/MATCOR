# MATCOR
One of the challenges in data analysics and data mining is the verification of the quality of the input data. Especially for large datasets it is not feasible to check each datum and automatization of this task is desirable. MATCOR (MATerials CORrelation) is a software to automate inter-database property validation/verification by finding the closest matching materials in two TestAPI supporting materials databases. The software in this repository is adjusted for materials verification between Materials Project (https://materialsproject.org/) and AFLOW (http://www.aflowlib.org/). 

If you use this software, please cite: 
Marquez Chavez, J. and Kiefer, B., "MATCOR, a program for the cross-validation of material properties between databases", Computational Materials Science (2020), https://doi.org/10.1016/j.commatsci.2020.110103. 

The implemented search strategy avoids the usage of ICSD tags. Instead, MATCOR searches for consistent compositions and space group to reduce the number of structure comparisons. This creates a pool of materials to the refinment of the comparison. After sorting this pool by magnitude difference between the proprty between reference and target database, MATCOR first attempts in a first pass a Hubbard-U consistent verification. If this fails in-consistent Hubbard-U verification is permitted. As a last step to avoid false positive materials verification, MATCOR uses the structure matcher utility available through pymatgen (https://pymatgen.org/).

In this repository you find the source code (MATCOR.py) and examples that we used to test MATCOR. 


Quickstart
To execute the program all you need is a list of valid materials identifiers ("reference_id_list") in AFLOW or Materials Project. Change the header of the python 3.8+ "MATCOR.py" script to match your needs, and specify a outputfile name "result_file". Simply execute within the pymatgen environment:

python3 MATCOR.py reference_id_list result_file

Files and user edits:


densitygit: base_property_MP_1 = 'density', base_property_MP_2 = '', base_property_AFLOW = 'density'
bandgapgit: base_property_MP_1 = 'band_gap', base_property_MP_2 = '', base_property_AFLOW = 'Egap'


