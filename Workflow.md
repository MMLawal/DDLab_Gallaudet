A workflow tutorial to identify small-molecule drug candidates.
1. Target Selection and Preparation:
•	Identify a well-defined disease target: Choose a protein or pathway with a clear role in the disease process. We are working on anticancer drug design and polo-like kinase 1 as a target protein.
•	Obtain a 3D structure of the target: The polo-box domain (PBD) portion of PLK1 is available from different X-ray crystallographic studies. The selected PBD structure from the protein data bank (PDB) is 4HCO (www.rcsb.org/structure/4HCO). 
•	Prepare the target structure: Clean and optimize the structure using software tools. Add missing residues and hydrogens, assign atom types, and remove unnecessary components, including water molecules and co-crystalized ligands. We have selected chain A of 4HCO with 4 missing residues in its loop and several buried terminal residues. Missing residue information is accessible in the raw X-ray crystal. Download and view in WordPad or TextEdit
Protein preparation for docking or screening can be done with various molecular modeling software, including MODELLER (salilab.org/modeller/), ChimeraX (www.cgl.ucsf.edu/chimerax/download.html) via integrated MODELLER, and online CHARMMGUI (www.charmm-gui.org/). Due to OS variants, here is an illustration involving CHARMMGUI for protein preparation after registering with your academic email.
a)	Click Input Generator 
b)	Click PDB Reader  
c)	Scroll down and input the PDB ID 4HCO 
d)	Click Next Step to see the details of the crystal structure and select PROA – protein chain A
 
e)	Click Next Step, it auto-selects the missing loop  
f)	Click Next Step and download step1_pdbreader.pdb 
Error: step1_pdbreader.pdb structure looks like , added missing residues failed. Redo the steps. 
To avoid complications from residue naming conventions across subsequent software for other steps, rename histidine (HSD) in step1_pdbreader.pdb to HIS.

2. Compound Library Selection:
•	Choose a compound library: This could be a free or commercially available library, a focused library based on known bioactivity, or a newly generated library. Our group is currently considering various natural compound databases like COCONUT (https://coconut.naturalproducts.net/) and SuperNatural 3.0 (https://bioinf-applied.charite.de/supernatural_3/index.php).
•	Filter the library: Consider factors like chemical space coverage, drug-likeness properties (molecular weight, logP, Tanimoto similarity, etc.), and other physicochemical properties based on available data. The “Advanced Search” option in COCONUT has several properties to trim down the screening based on desired chemical components. Our team members have screened COCONUT after carefully selecting five natural products identified through a literature search of experimental studies with characteristics as candidate PLK1 inhibitors. Students performed a structural similarity search in the COCONUT to identify other compounds like these 5 parent natural compounds. The protocol showed small molecules totaling 258 compounds, requiring further mining to remove isomers. The final list includes 25, 25, 25, 25, and 6 compounds as potential PLK1-PBD inhibitors grouped under Quercetin, Xanthatin, Chaetoglobosin E, Luteolin, and Resveratrol derivatives, respectively. An illustration of the search and filtering is available herein.
a)	Click the COCONUT link above, type luteolin  
Dave guys can help with other steps 

3. Virtual Docking:
•	Docking software selection: Choose a suitable docking program based on your needs, such as AutoDock Vina, Glide, or DOCKER. For broad accessibility, CB-Dock2 (https://cadd.labshare.cn/cb-dock2/index.php) is our choice of software because its parameterization is based on AutoDock Vina and requires few parameters to execute. 
•	Define the docking site: Specify the binding pocket or region of interest on the target protein where potential drug candidates are expected to bind. We currently use a blind docking approach to explore all possible binding sites and select the most prominent since various experimental studies have established about 3 promising binding cavities for PBD. 
•	Docking simulation: Perform docking simulations of the filtered compound library against the prepared target structure. This generates scores or binding affinities reflecting the predicted interaction strength between each compound and the target. 
a)	Here is a tutorial on docking with CB-Dock2. 
b)	Go to the link and click Dock 
c)	Upload the prepared Protein and Ligand
d)	Click    
Note: This protocol is slow for multiple ligands, students can install AutoDock Vina locally to do multiple compound docking. For MacBook users, a summary of the steps is available on our GitHub page (https://github.com/MMLawal/DDLab_Gallaudet/blob/main/AutoDock_Vina/install_vina.txt).
•	Analysis and Prioritization:
a)	Score analysis: Analyze the docking scores or binding affinities. A higher score typically indicates a more favorable interaction. Set a threshold or use ranking methods to identify the top-scoring compounds.
b)	Visual inspection: Analyze the docked poses of high-scoring compounds to assess their interactions with the target's binding pocket. CB-Dock2 has a visualization interface, other tools like PyMOL can be used for visualization.

4. ADMET/S prediction
Evaluate selected compounds for ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) properties using ADMETLab 2.0 or SwissADME and stability (S) with density functional theory (DFT) calculations in ORCA.
•	Here is a step-by-step guide on how to run the ADMETlab 2.0 online.
a)	Access the website: Visit admetmesh.scbdd.com to access the ADMETlab 2.0 web server.
b)	Choose prediction mode: Select either "Screening" for batch (multiple compounds) predictions or "Single-molecule evaluation" for individual molecules.
c)	Input molecules: To screen, paste a list of SMILES strings or upload an SDF or TXT file containing molecules (ensure no headers or indexes).
d)	Submit task: Click the "Submit" button to initiate predictions.
e)	View results: After screening, the results table displays SMILES, 2D structures, and a "View" button for each molecule. Click "View" to access detailed predictions for a specific molecule. Download the results as a CSV file for further analysis. For a single-molecule evaluation, an interface displayed the prediction.
•	Running SwissADME
a)	Access the SwissADME website: Open your web browser and go to www.swissadme.ch
b)	Input your molecules: For a single molecule, use the molecule sketcher on the left to draw the structure or paste the molecule's SMILES. Other acceptable molecule file formats include SDF, MOL, and PDB extensions. Click "Transfer to SMILES" to generate the SMILES code. For multiple molecules, paste a list of SMILES codes (one per line) into the text box on the right or include molecule names separated by spaces after each SMILES.
c)	Run the calculations: Click the "Run" button (red) below the SMILES list. SwissADME will process the molecules and display the results on the same page.
d)	View the results: The results are organized into several tabs containing relevant information and visualizations for interpreting the predictions. These headings include physicochemical properties, lipophilicity, water solubility, pharmacokinetics, drug-likeness, and medicinal chemistry.
e)	Export results: Click the CSV icon to download a comma-separated text file or click the clipboard icon to copy the results for pasting elsewhere.
•	Running ORCA
a)	Install (https://orcaforum.kofo.mpg.de/app.php/portal)
b)	Optimize ligands, then run frequency calculation on the optimized molecule to obtain HOMO and LUMO for band gap estimation to depict S (stability)
c)	Build your initial structure in Avogadro and save. 
 create a rough structure or import ligand
   clean it up.
   calculation setup.
 click on Generate…  and save
d)	Run ORCA for all the input files prepared in step 2. 
ORCA execution on a local terminal (MacOS)
 
ORCA execution on a remote server Darwin at UDEL.
 
e)	Analyze your results and enter information in the spreadsheet. 
Use a command prompt (any Terminal) to extract the HOMO and LUMO in ORCA output with the below, whereby the first command reps HOMO and the second is LUMO.
grep "ORBITAL ENERGIES" z8xxx.out -A40 $1 | grep '2.0000' |tail -1 |awk '{print $4}' 
grep "ORBITAL ENERGIES" z8xxx.out -A40 $1 | grep '0.0000' |head -1 |awk '{print $4}'
The likelihood of SOMO/SUMO is high if the molecule has a doublet multiplicity. In such case, we will have the tweak HOMO command to grep '1.0000' or something like that.
You also view the output, search ORBITAL ENERGIES, scroll to find your HOMO and LUMO
  orbital indices 34 and 35 are HOMO and LUMO, respectively.


5. Binding energy calculation with webservers 
Beyond docking scores, the binding energy of ligands to protein could be refined with several approaches, including using freely accessible tools. iSome web servers for ligand-protein binding energy calculations.
a)	fastDRH  -  http://cadd.zju.edu.cn/fastdrh/submit - output is like MMPB/GBSA/ (it looks promising because it uses MMPB/GBSA protocol)
b)	ProDIGY - https://bianca.science.uu.nl/prodigy/lig - output is like docking score 
c)	BAPPL server - http://www.scfbio-iitd.res.in/software/drugdesign/bappl.jsp#anchortag – input file seems complicated
d)	CSM-Lig - https://biosig.lab.uq.edu.au/csm_lig/related - input file requires careful preparation.

•	Running CSM-Lig
a)	Download the docked protein-ligand (PDB format) from CB-Dock2 autodocking.
b)	The ligand (LIG) output from CB-Dock2 no longer has hydrogen atoms and some H atoms are also missing in the protein, add the missing hydrogens.
c)	Open the downloaded docked prot-lig.pdb with ChimeraX (if you don’t have https://www.cgl.ucsf.edu/chimerax/download.html). Then click on ToolsStructure EditingAddH
  
It will pop a window.
d)	Click okay and save in the preferred directory. Doing this for several complexes will be quite a task, but it is possible. 
e)	Convert ligands to SMILES.
f)	Upload in CSM-Lig and run prediction.
 
g)	You can use the docked prot-lig.pdb as downloaded from CB-Dock2 without adding the hydrogen back to note the variation.

6. MD and MMPBSA
Our group often uses NAMD software for molecular dynamics (MD) simulations, but we are currently using AMBER (https://ambermd.org/) because of its robustness. The AMBER tutorial is here: https://ambermd.org/tutorials/ while examples of analysis for various metrics are offered through AmberTools (https://amberhub.chpc.utah.edu/). In most cases, Advanced computing resources or supercomputers facilitate MD simulations. The accessible resources for MD simulations in our lab are Anvil (https://www.rcac.purdue.edu/knowledge/anvil) and NCSA Delta (https://docs.ncsa.illinois.edu/systems/delta/en/latest/index.html). See the summary below to set up and simulate a protein-ligand complex using classical MD (cMD) simulations. Amber is available for all users in Anvil, its installation on NCSA-Delta is still under consideration but an individual installation guide (GPU and MPI-enabled) is available here (https://github.com/MMLawal/DDLab_Gallaudet/blob/main/AMBER/install_amber.txt). 
•	Ligand and protein preparation
a)	Use the docked pose of the desired ligand retrieved from CB-Dock2 docking, other docking software, or co-complexed 3D protein-ligand structure from an experiment.
b)	Add hydrogen back to the ligand in Avogadro by clicking on BuildAdd HydrogensFileSave As… to save in a preferred folder.
 	  
c)	Upload your ligand and protein to your remote work directory.
d)	Load required modules to activate AmberTools/Amber 
e)	Generate forcefield parameter for ligand with antechamber module (https://ambermd.org/tutorials/basic/tutorial4b/index.php) “antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -j 5 -at gaff -c bcc -nc 0” then “parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod.” Change “-nc 0” to whatever charge your ligand has. 
f)	Prepare your protein – the structure uploaded for docking. Remove all the hydrogen atoms to avoid imminent connectivity issues, tleap will add them later and assign appropriate protonation. Type “reduce -Trim step1_protein.pdb > rec.pdb” to remove hydrogen.
g)	Final tleap.in file to prepare a solvated protein-ligand complex for MD simulation. “tleap -s -f tleap.in” The first 5 lines are the forcefield parameters for protein, ligand, water, and ionizing atoms. Delete the comments (!to be used later for MMGBSA) in the file if it pops an error. The total charge of this protein (PBD) is 0, hence the equimolar 0.15M NaCl addition. See this tutorial for clarity https://ambermd.org/tutorials/basic/tutorial8/index.php. 
 
•	MD simulations
a)	Simulating with AMBER usually follows this protocol minimizationheating/equilibration production (https://ambermd.org/tutorials/basic/tutorial0/index.php). The parameter and coordinate files are step3_input.parm7 and step3_input.rst7 required for the simulation. Other extension formats for parameter files are parm, prmtop, and top while other coordinate formats are inpcrd, crd, and rst. Based on a recent study, here are the input files for cMD simulation involving min1.in, min2.in, heat.in, equil.in, and prod.in (200 ns) with a submitting script amber_run.sh. You must learn to work with Linux and understand a few batch scripting concepts. 
Examples:
Submit your simulation with sbatch amber_run.sh.
Cancel with scancel jobid
Check your submitted job progress with squeue -u username
Or all jobs on the infrastructure with squeue

b)	Read/search the Amber manual (https://ambermd.org/doc12/Amber24.pdf) to understand the meaning of each parameter in the input files.
   
   
   

c)	Analyses – there are various metrics to depict information from the simulated trajectories. Usually, water and ions co-simulated with the protein-ligand system will be stripped off to reduce the file size and speed up the post-simulation analyses. Here is an example of a strip.in file to execute using the command cpptraj strip.in
 

•	MMGBSA/PBSA calculation 
a)	The AmberTools has several packages for various analyses, including the cpptraj (https://amberhub.chpc.utah.edu/cpptraj/) and MMPBSA.py for end-point binding free energy calculation (https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.php). See template mmgbsa.in and mmgbsa.sh files. You can create a separate directory for your analyses.
  
b)	Installing AmberTools in your HPC home directory, load Anaconda with "module load anaconda" or whatever module spider anaconda suggests. Anaconda is available in most computing infrastructures. Follow the instructions in this link https://ambermd.org/GetAmber.php. 
 

7. Metadynamics
a)	The variance to this advanced MD simulation is in the sampling approach (https://ambermd.org/tutorials/advanced/tutorial31/index.php). Below are the input files generated with CHARMGUI with some modifications for the center of mass (COM) distance between two residues as the collective variable (CV). The files are step4.0_minimization.mdin, step4.1_equilibration.mdin, cv.dat, and step5_production.mdin with the submission script metad.sh
Uncommenting the details in the script is not required, it will run as it is.
  
 
 
 

b)	Analysis – run this command to obtain the PMF nfe-umbrella-slice 'nfe-bias.nc' > FE.dat 
Plot with Gnuplot or matplotlib in Python/Conda on the remote or your laptop. See the example below to run with python run_FE.py
 
![image](https://github.com/MMLawal/DDLab_Gallaudet/assets/64615834/dc892f43-d646-47a2-8150-75fd4ebee0ea)
