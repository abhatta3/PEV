#https://pythonhosted.org/pyteomics/data.html
import sys, getopt
import time


def print_help():
    print 'python ./Codes/Make_list.py -i <input event file> -o <outputfile (PSMs)> -d <database search file (from ENOSI/MSGF+.../mztab)>\n'
    print 'Optional arguments [For following column numbers the first column number is 0]:'
    print '-a : column number for PSM tag; default 0'
    print '-b : column number for spectrum index number in database search file; default 2'
    print '-c : column number for reported peptide in database search file; default 1'
    print '-x : Event ID column number in inputfile; default 0'
    print '-y : Event type column number in inputfile; default 1'
    print '-z : Event peptide column number in inputfile; default 2'
    

def main(argv):
    
    inputevent = ''
    spec_outfile = ''
    scanf = ''
    d_pep=1
    d_index=2
    d_fname=0
    i_ev=0
    i_seq=2
    i_type=1


    try:
        opts, args = getopt.getopt(argv,"hi:o:d:a:b:c:x:y:z:")
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print_help()
            sys.exit()
        elif opt in "-i":
            inputevent = arg
        elif opt in "-o":
            spec_outfile = arg
        elif opt in "-d":
            scanf = arg

        elif opt in "-c":
            d_pep = int(arg)
        elif opt in "-b":
            d_index = int(arg)
        elif opt in "-a":
            d_fname = int(arg)
        elif opt in "-x":
            i_ev = int(arg)
        elif opt in "-y":
            i_type = int(arg)
        elif opt in "-z":
            i_seq = int(arg)


        
    if inputevent == '':
        print("Error: Wrong parameter. -i event file is missing") 
        print_help()
        sys.exit()
    
    if scanf == '':
        print("Error: Wrong parameter. -d database search file is missing")  
        print_help()
        sys.exit()
    if spec_outfile == '':
        print("Error: Wrong parameter. -o output file name is missing")  
        print_help()
        sys.exit()

    output=open(spec_outfile,"w")
    pep={}
    event_type={}
    novel={}
    with open(inputevent,"r") as inputfile:
        for line in inputfile:
            if '#' in line:
                continue
            line=line.strip().split('\t')
            if len(line)<1:
                continue
            if '.' in line[i_seq]:
                pepseq=line[i_seq].split('.')[1].translate(None, '1234567890._:-*!@#$?')
            else:
                pepseq=line[i_seq].translate(None, '1234567890._:-*!@#$?')
                
            novel[pepseq]=line[i_seq]
            if pep.has_key(pepseq):
                old=pep[pepseq]
                old.append(line[i_ev])
                pep[pepseq]=old
            else:
                pep[pepseq]=[line[i_ev]]
            event_type[pepseq]=line[i_type]
    output.write("#Event\tIndex_number\tFile\tPeptide\tEvent_Type\n") 
           
    with open(scanf,"r") as sf:
        for line in sf:
            if 'PSM' not in line:
                continue
            line=line.strip().split('\t')
            if len(line)<8:
                continue
            if '.' in line[d_pep]:
                pepseq=line[d_pep].split('.')[1].translate(None, '1234567890._:-*!@#$?')
            else:
                pepseq=line[d_pep].translate(None, '1234567890._:-*!@#$?')
            fname=line[d_fname]
            if '=' in line[d_index]:
                index=line[d_index].split('=')[1]
            else:
                index=line[d_index]
                
            if pep.has_key(pepseq):
                pev=pep[pepseq]
                output.write("%s" %(pev[0]))
                if len(pev)>1:
                    for itm in pev[1:]:
                        output.write(",%s" %(itm))
                output.write("\t%s\t%s\t%s\t%s\n" %(index,fname,novel[pepseq],event_type[pepseq]))
    output.close()

if __name__ == "__main__":
    start_time = time.time()
    main(sys.argv[1:])
    print("Executes in --- %s seconds ---" %(time.time() - start_time))
