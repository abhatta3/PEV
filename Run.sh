#!/bin/bash
echo "Enter Proteogenomics list file name >"
read proteog

echo "Enter Database search file name >"

read search


echo "Enter mgf file name >"

read mgffile

echo "Enter Output file name >"

read outfile

echo "Proteogenomics list file name > $proteog"
echo "Database search file name > $search"
echo "mgf file name > $mgffile"
echo "Output file name > $outfile"

echo "Computing ..."
echo "Running Step 1 of 8: Directly compare for exact match of novel in reference>"

python ./Codes/novel_filtering.py -i $proteog -o $outfile"step1"

echo "Running Step 2 of 8: Make simple PSM file>"



extension=`echo "$search"|awk -F . '{print $NF}'`

if [ $extension == "mztab" ]; then
    python ./Codes/Make_list_from_mztab.py -i $outfile"step1" -o $outfile"step2" -d $search 
else
    python ./Codes/Make_list.py -i $outfile"step1" -o $outfile"step2" -d $search 
fi
     
echo "Running Step 3 of 8: Exctract the relevant spectra into a single .mgf file>"

python ./Codes/extract_spec.py -i $outfile"step2" -s $mgffile -o $outfile"step3"

echo "Running Step 4 of 8: Get the support sequence from comparing the spectra, b-ions, y-ions>"

python ./Codes/get_keywords.py -i $outfile"step3.mgf.map.txt"  -s $outfile"step3.mgf" -o $outfile"step4"

echo "Running Step 5 of 8: Get the alternative sequences with a report and fasta sequence file>"

python ./Codes/modification_mutation_search.py -i $outfile"step4" -p ./Data/alternate.fasta.tab  -o $outfile"step5"

echo "Running Step 6 of 8: Run Msgf+ search on fasta and selected single .mgf>"

java -Xmx3500M -jar ./msgfplus/MSGFPlus.jar  -d  $outfile"step5.fa" -s  $outfile"step3.mgf"  -o $outfile"step6.mzid" -t 50ppm -m 0 -inst 0 -e 1 -ti -1,2 -ntt 2 -tda 0 -minLength 8 -maxLength 40 -n 5 -thread 7  -mod  ./Data/Modifications_msgf.txt 

echo "Running Step 7 of 8: Convert mzid search result to .tsv>"

python ./Codes/MzidToTsv_withmod.py $outfile"step6.mzid" $outfile"step7.tsv"

echo "Running Step 8 of 8: Run rescoring to get log-odds score tailing true/false proteogenomics call>"

python ./Codes/Rescore.py -i $outfile"step7.tsv" -s $outfile"step3.mgf.map.txt" -r $outfile"step5.report" -o $outfile"step8"

echo "Execution complete....Find the final output in "$outfile"step8"
