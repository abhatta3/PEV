#https://pythonhosted.org/pyteomics/data.html
import os
import os.path
import sys, getopt
import time


def print_help():
    print 'python ./Codes/extract_spec.py -i <input PSMs file (From Make_list.py)> -o <outputfile (PSMs)> -s <Spectra files path>'


def extracter(argv):
    inputf = ''
    spec_outfile = ''
    spec_dir = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:s:",["ifile=","ofile=","spath="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print_help()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputf = arg
        elif opt in ("-o", "--ofile"):
            spec_outfile = arg
        elif opt in ("-s", "--spath"):
            spec_dir = arg

    if inputf == '':
        print("Error: Wrong parameter. -i input file is missing")  
        print_help()
        sys.exit()
    
    
    if spec_outfile == '':
        print("Error: Wrong parameter. -o output file name is missing")  
        print_help()
        sys.exit()
    
    
    if spec_dir == '':
        print("Error: Wrong parameter. -s Spectra files path is missing")  
        print_help()
        sys.exit()
        
    if  '.mgf' not in spec_outfile:
        spec_outfile+='.mgf'
        
    

    soutput=open(spec_outfile,"w")
    spec_map_outfile=spec_outfile+'.map.txt' 
    output=open(spec_map_outfile,"w")
 
    print("Output Files:\nCombined spectra file:%s\nPSM File:%s\n" %(spec_outfile,spec_map_outfile))
    
    event_scan={}
    
    with open(inputf,"r") as inputfile:
        for line in inputfile:
            if '#' in line:
                continue
            line=line.strip().split('\t')
            index=int(line[1])
            fname=line[2]
            pev=line[0]
            peptide=line[3]
            etype=line[4]
            if event_scan.has_key(fname):
                sdic=event_scan[fname]
                sdic[index]={'pev':pev,'pseq':peptide,'etype':etype}
                event_scan[fname]=sdic
            else:
                temp={'pev':pev,'pseq':peptide,'etype':etype}
                sdic={index:temp}
                event_scan[fname]=sdic
    
    newscanno=0
    output.write("#Event\tIndex_number\tFile\tPeptide\tEvent_Type\n")            
    for fname in event_scan:
        print fname
        scanlist=event_scan[fname]
        if '.mgf' not in spec_dir:
            PATH=spec_dir+fname
        else:
            PATH=spec_dir
        flag=0
        if os.path.isfile(PATH) and os.access(PATH, os.R_OK):
            with open(PATH,"r") as mgf:
                scanindex=0
                for line in mgf:
                    if len(line)>0:
                        if 'BEGIN IONS' in line:
                            if scanlist.has_key(scanindex):
                                pev=event_scan[fname][scanindex]['pev']
                                output.write("%s\t%d\t%s\t%s\t%s\n" %(pev,newscanno,spec_outfile,event_scan[fname][scanindex]['pseq'],event_scan[fname][scanindex]['etype']))
                                newscanno+=1
                                soutput.write(line)
                                flag=1
                            scanindex+=1
                        elif flag==1 and 'END IONS' not in line:
                            soutput.write(line)
                        elif flag==1 and 'END IONS' in line:
                            soutput.write(line)
                            flag=0
    soutput.close()
    output.close()
                             

if __name__ == "__main__":
    start_time = time.time()
    extracter(sys.argv[1:])
    print("Executes in --- %s seconds ---" %(time.time() - start_time))
