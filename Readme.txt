Current steps
======================================
Step 1: Directly compare for exact match of novel in reference
help: python ./Codes/novel_filtering.py -h
python ./Codes/novel_filtering.py -i ~/Data_for_validation/pride/peptide_list.txt -o ~/Data_for_validation/pride/peptides.txt

Step 2: If database search result in .tsv format then (A); in .mztab then (B)

A. python ./Codes/Make_list.py -i ~/Data_for_validation/pride/peptides.txt -o ~/Data_for_validation/pride/map_temp.txt -d ~/Data_for_validation/pride/PRIDE_Exp_Complete_Ac_30823.pride.tsv 

B. python ./Codes/Make_list_from_mztab.py -i ~/Data_for_validation/pride/peptides.txt -o ~/Data_for_validation/pride/map_temp.txt -d ~/Data_for_validation/pride/PRIDE_Exp_Complete_Ac_30823.pride.mztab 


Step 3: Exctract the relevant spectra into a single .mgf file [If from single mgf file then give the file name with full path name else just the full path name]

python ./Codes/extract_spec.py -i ~/Data_for_validation/pride/map_temp.txt -s ~/Data_for_validation/pride/PRIDE_Exp_Complete_Ac_30823.pride.mgf -o ~/Data_for_validation/pride/map.txt


Step 4: Get the support sequence from comparing the spectra, b-ions, y-ions

python ./Codes/get_keywords.py -i ~/Data_for_validation/pride/map.txt.mgf.map.txt  -s ~/Data_for_validation/pride/map.txt.mgf -o ~/Data_for_validation/pride/keywords.txt

Step 5: Get the alternative sequences with a report and fasta sequence file

python ./Codes/modification_mutation_search.py -i ~/Data_for_validation/pride/keywords.txt -p ./Data/alternate.fasta.tab  -o ~/Data_for_validation/pride/alternatives.txt

Step 6: Run Msgf+ search on fasta and selected single .mgf.

java -Xmx3500M -jar /Users/anb013/Source_Code_GIT/Proteogenomics_Validation/msgfplus/MSGFPlus.jar  -d  /Users/anb013/Data_for_validation/pride/alternatives.txt.fa -s  /Users/anb013/Data_for_validation/pride/map.txt.mgf  -o /Users/anb013/Data_for_validation/pride/alternative_updated.mzid -t 50ppm -m 0 -inst 0 -e 1 -ti -1,2 -ntt 2 -tda 0 -minLength 8 -maxLength 40 -n 5 -thread 7  -mod  /Users/anb013/Source_Code_GIT/Proteogenomics_Validation/Data/Modifications_msgf.txt 

Step 7: Convert mzid search result to .tsv

python ./Codes/MzidToTsv_withmod.py ~/Data_for_validation/pride/alternative_updated.mzid ~/Data_for_validation/pride/alternative.tsv

Step 8: Run rescoring to get log-odds score tailing true/false proteogenomics call

python ./Codes/Rescore.py -i ~/Data_for_validation/pride/alternative.tsv -s ~/Data_for_validation/pride/map.txt.mgf.map.txt -r ~/Data_for_validation/pride/alternatives.txt.report -o ~/Data_for_validation/pride/output.txt



