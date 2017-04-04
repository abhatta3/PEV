import sys, getopt
import time


def print_help():
    print 'python ./Codes/novel_filtering.py -i <input event file> -o <outputfile event file novel peptides>'
    

def main(argv):
    
    inputevent = ''
    outfile = ''


    try:
        opts, args = getopt.getopt(argv,"hi:o:")
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
            outfile = arg

        
    if inputevent == '':
        print("Error: Wrong parameter. -i event file is missing") 
        print_help()
        sys.exit()
    
    if outfile == '':
        print("Error: Wrong parameter. -o output file name is missing")  
        print_help()
        sys.exit()

    novel={}
    with open(inputevent,"r") as inputfile:
        for line in inputfile:
            if '#' in line:
                continue
            line=line.strip().split('\t')
            if len(line)<3:
                continue
            if '.' in line[2]:
                pepseq=line[2].split('.')[1].translate(None, '1234567890._:-*!@#$?')
            else:
                pepseq=line[2].translate(None, '1234567890._:-*!@#$?')
                
            novel[line[0]]={'type':line[1], 'seq':pepseq,'n':True,'id':'','seqold':line[2]}

    with open("./Data/alternate.fasta80.tab","r") as sf:
        for line in sf:
            if '#' in line:
                continue
            line=line.strip().split('\t')
            if len(line)<2:
                continue
            for event in novel:
                if novel[event]['n']==True:
                    if novel[event]['seq'] in line[1]:
                        novel[event]['n']=False
                        novel[event]['id']=line[0]

    output=open(outfile,"w")
    notnovelf=outfile+"_not_novel.txt"
    output_ref=open(notnovelf,"w")

    output.write("#Event\tEvent_Type\tPeptide\n") 
    output_ref.write("#UniProtID\tEvent\tEvent_Type\tPeptide\n") 
                        
    for event in novel:
        if novel[event]['n']==True:
            output.write("%s\t%s\t%s\n" %(event,novel[event]['type'],novel[event]['seqold']))
        else:
            output_ref.write("%s\t%s\t%s\t%s\n" %(novel[event]['id'],event,novel[event]['type'],novel[event]['seqold']))
    output.close()
    output_ref.close()
    
    print ("Novel data are in %s, ready for validation\n" %(outfile))
    print ("Input Novel that mapped to Reference are in %s, discarded from validation\n" %(notnovelf))

if __name__ == "__main__":
    start_time = time.time()
    main(sys.argv[1:])
    print("Executes in --- %s seconds ---" %(time.time() - start_time))
    