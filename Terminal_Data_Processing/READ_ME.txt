====================================================
ExtractBINOUT.sh
====================================================
Place in main directory.
Extract glstat and matsum files from each impact orientation simulation, give them descriptive names, and then store each within the folders COLLECTED_GLSTAT and COLLECTED_MATSUM, respectively. 

====================================================
ExtractOutput.sh
====================================================
Place in main directory.
Creates the folder "collected_output" which stores separate folders corresponding to each impact orientation. Data on particle fracture is stored within each impact orientation folder. 


================================================
ExtractGLSTAT.py
================================================
Place this script in the COLLECTED_GLSTAT folder created by the executable Extract_BINOUT_files.sh
Results in the csv file glstat_summary.csv. Descriptive labels are added for each column in the first row of the csv.