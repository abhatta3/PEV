import math
import sys, getopt
import time

def print_help():
    print 'python ./Codes/Rescore.py -i <MSGF+ score .tsv file> -o <outputfile> -s <PSM File from extract_spec.py> -r <Report file from modification_mutation_search.py>'

def rescore(argv):
    ifile = ''
    outputfile = ''
    scan_map_file = ''
    reportfile=''
    try:
        opts, args = getopt.getopt(argv,"hi:o:s:r:",["ifile=","ofile=","sfile=","rfile="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print_help()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-s", "--sfile"):
            scan_map_file = arg
        elif opt in ("-r", "--rfile"):
            reportfile = arg

    if ifile == '':
        print("Error: Wrong parameter. -i input file is missing")  
        print_help()
        sys.exit()
    
    
    if outputfile == '':
        print("Error: Wrong parameter. -o output file name is missing")  
        print_help()
        sys.exit()
    
    
    if scan_map_file == '':
        print("Error: Wrong parameter. -s PSM file is missing")  
        print_help()
        sys.exit()

    if reportfile == '':
        print("Error: Wrong parameter. -r report file is missing")  
        print_help()
        sys.exit()
    
    event_i={}
    with open(reportfile) as statf:
        for line in statf:
            if ',' in line:
                line=line.strip().split(',')
                for w in line:
                    if w!='':
                        event_i[int(w)]=1

    scan_event={}
    event_list={}
    
    with open(scan_map_file,"r") as mapf:
        for line in mapf:
            if '#' in line:
                continue
            line=line.strip().split("\t")
            if len(line)<2:
                continue
            event_list[line[0]]={'type':line[4],'seq':line[3]}
            scan_event[line[1]]={'ev':line[0],'novel':2.0,'modified':2.0}

    max_mod={}
    max_nov={}
    max_pep={}
      
    with open(ifile,"r") as data:
        for line in data:
            if '#' == line[0] or 'XXX_' in line:
                continue
            line=line.strip().split('\t')
            indexno=line[1].split('=')[1]
            evalue=float(line[12])
            tagl=line[9]
            if scan_event.has_key(indexno):
                s_event=scan_event[indexno]['ev']
                '''
                if max_pep.has_key(s_event):
                    if max_pep[s_event]['e']>evalue:
                        max_pep[s_event]['e']=evalue
                        max_pep[s_event]['p']=line[8]
                else:        
                    max_pep[s_event]={'e':evalue,'p':line[8]}
                '''
                scan_tag_novel='|Novel'
                scan_tag_modified='|Modified'
    
                if scan_tag_novel in tagl:
                    if scan_event[indexno]['novel']>evalue:
                        scan_event[indexno]['novel']=evalue
    
                    if max_nov.has_key(s_event):
                        if max_nov[s_event]>evalue:
                            max_nov[s_event]=evalue
                    else:        
                        max_nov[s_event]=evalue
    
                elif scan_tag_modified in tagl:
                    if max_mod.has_key(s_event):
                        if max_mod[s_event]>evalue:
                            max_mod[s_event]=evalue
                            max_pep[s_event]=line[8]
                    else:        
                        max_mod[s_event]=evalue
                        max_pep[s_event]=line[8]
                    
                    if scan_event[indexno]['modified']>evalue:
                        scan_event[indexno]['modified']=evalue
            
    
    score_by_event={}

    for indexno in scan_event:
        evalue=scan_event[indexno]['modified']
        s_event=scan_event[indexno]['ev']
        new_n_eval=scan_event[indexno]['novel']
        
        if new_n_eval!=2.0 and evalue!=2.0:
            score=float(-1.0*float(math.log(new_n_eval/evalue)))
            
            if score_by_event.has_key(s_event):
                if score>score_by_event[s_event]:
                    score_by_event[s_event]=score
            else:
                score_by_event[s_event]=score
        
    
    for s_event in max_mod:
        if not score_by_event.has_key(s_event) and not max_nov.has_key(s_event):
            evalue=max_mod[s_event]
            score_by_event[s_event]=float(-1.0*float(math.log(1.0/evalue)))
    
    for s_event in max_nov:
        if not score_by_event.has_key(s_event):
            evalue=max_nov[s_event]
            score_by_event[s_event]=float(-1.0*float(math.log(evalue)))
    
    with open(outputfile,"w") as output:
        output.write("#Event_ID\tProteogenomics_Peptide\tPEV_sequence\tVariant_Type\tPEV_Score\tIdentification_Class\n")
        for ev in event_list:
            if score_by_event.has_key(ev):
                score=score_by_event[ev]
            else:
                score=999999.9 #eventlist[ev]
            
            if max_pep.has_key(ev):
                newseq=max_pep[ev]
            else:
                newseq='-'
    
            if score>0.0:
                category='True'
            else:
                category='False'
                
            mev=event_list[ev]['type']
            seq=event_list[ev]['seq']
    
            evlist=ev.split(',') 
            for fev in evlist:
                k_ev=int(fev)
                if event_i.has_key(k_ev):
                    if score!=999999.9:
                        output.write("%s\t%s\t%s\t%s\t%f\t%s\n" %(fev,seq,newseq,mev,score,category))
                    else:
                        output.write("%s\t%s\t-\t%s\tNA\tTrue\n" %(fev,seq,mev))
                else:
                    output.write("%s\t%s\t-\t%s\tNA\tTrue\n" %(fev,seq,mev))



if __name__ == "__main__":
    start_time = time.time()
    rescore(sys.argv[1:])
    print("Executes in --- %s seconds ---" %(time.time() - start_time))
