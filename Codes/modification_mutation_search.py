#USAGE: python modification_mutation_search.py outputcptac.txt alternate.fasta80.tab new_mod_alter_output80.txt
#PATH:/Users/anb013/Source_Code_GIT/Proteogenomics_Validation/Data

import sys, getopt
import time

#delta=50.0 #ppm
#keylenth=4
#most_common_mod={15.994915:'M',14.015650:'K',-17.026549:'Q',79.966331:'S,T,Y',43.005814:'*',42.010565:'*',0.984016:'N,Q'}
# Oxidation	+15.994915	M	OPTIONAL
# Lysine Methylation	+14.015650	K	OPTIONAL
# Pyroglutamate Formation	-17.026549	Q	OPTIONAL, N-TERMINAL
#Phosphorylation	+79.966331	STY	OPTIONAL
#N-terminal Carbamylation	+43.005814	*	OPTIONAL, N-TERMINAL
#N-terminal Acetylation	+42.010565	*	OPTIONAL, N-TERMINAL
#Deamidation	+0.984016	NQ	OPTIONAL


def print_help():
    print 'python ./Codes/modification_mutation_search.py -i <input file (From get_keywords.py)> -o <outputfile (PSMs)> -p <Protein Database file> -k <Minimum Key length; default=4> -m <PArent mass tolerance in ppm; defualt is 50>'

#(theoretical mass) / ((1/errorPPM)*1,000,000) = error(Da)event
    

def comput_mass(peptide,mass_dic):
    total_m=19.0
    for d in peptide:
        if mass_dic.has_key(d):
            total_m+=mass_dic[d]['S']
            if d=="C":
                total_m=total_m+57.0            
    return total_m

def comput_delta(total_m,delta):
    val=total_m/((1/delta)*1000000)
    return val
    
def build_key_index(k_index_dic,next_new_node,kw,event):
    node_i=0
    kw=kw+'#'
    for item in kw:
        if item=='#':
            if k_index_dic.has_key(node_i):
                node=k_index_dic[node_i]
                if node.has_key(item):
                    break
                else:
                    node[item]={'k':kw[:-1],'e':event}
                    k_index_dic[node_i]=node
            else:
                node={}
                node[item]={'k':kw[:-1],'e':event}
                k_index_dic[node_i]=node
        else:
            if k_index_dic.has_key(node_i):
                node=k_index_dic[node_i]
                if node.has_key(item):
                    node_i=node[item]
                else:
                    node[item]=next_new_node
                    k_index_dic[node_i]=node
                    node_i=next_new_node
                    next_new_node+=1
            else:
                node={item:next_new_node}
                k_index_dic[node_i]=node
                node_i=next_new_node
                next_new_node+=1
    kw_data = (k_index_dic,next_new_node)        
    return kw_data

def create_pep_table(inputf,keylenth,delta):
    peptide_dic={}
    k_index_dic={}
    pep_dic_by_event={}
    next_new_node=1
    
    events_novel={}
    with open(inputf,"r") as inputfile:   
        for line in inputfile:
            if '#' in line:
                continue
            line=line.strip().split('\t')
            events_novel[line[0]]=line[2]
            
            plist=line[5].split('.')
            kw=max(plist, key =len)
            kw_l=len(kw)
            if pep_dic_by_event.has_key(line[0]):
                if len(pep_dic_by_event[line[0]])<kw_l:
                    pep_dic_by_event[line[0]]=kw
                    peptide_dic[line[0]]={'NP':line[2],'rN':line[1],'scan':line[8],'pos':line[5].index(kw),'lsupport':line[9],'pmass':float(line[3])}
            else:
                pep_dic_by_event[line[0]]=kw
                peptide_dic[line[0]]={'NP':line[2],'rN':line[1],'scan':line[8],'pos':line[5].index(kw),'lsupport':line[9],'pmass':float(line[3])}
    count=0
    total=0
    for ln in pep_dic_by_event:
        kw=pep_dic_by_event[ln]
        if len(kw)<keylenth:
            pep=peptide_dic[ln]['NP']
            start_position=peptide_dic[ln]['pos']
            alter_support=peptide_dic[ln]['lsupport']
            remaining_size=keylenth-len(kw)
            
            left_b=start_position-remaining_size
            
            if left_b<0:
                left_b=0
            right_b=start_position+keylenth
            if right_b>len(pep):
                right_b=len(pep)

            l_size=start_position
            r_size=len(pep)-start_position-len(kw)

            count_l=0
            for i in range(left_b,left_b+keylenth):
                if alter_support[i]=='1':
                    count_l+=1
                
            count_r=0
            for i in range(start_position,right_b):
                if alter_support[i]=='1':
                    count_r+=1

            if count_r==count_l:
                if l_size>r_size:
                    w_start=left_b
                    w_end=w_start+keylenth
                else:
                    w_start=start_position
                    w_end=w_start+keylenth
            elif count_r>count_l:
                w_start=start_position
                w_end=w_start+keylenth
            elif count_l>count_r:
                w_start=left_b
                w_end=w_start+keylenth

            kw=pep[w_start:w_end]

            count+=1
            
        pep_dic_by_event[ln]=kw
        

        total+=1

    for e in pep_dic_by_event:
        kword=pep_dic_by_event[e]#[-keylenth:] #restrict all keywords to length 4
        temp_s=build_key_index(k_index_dic,next_new_node,kword,e)
        k_index_dic=temp_s[0]
        next_new_node=temp_s[1]
    temp_dic_str=(k_index_dic,peptide_dic,events_novel)

    print count,total
    return temp_dic_str

def read_massdic(massfile):
    mass_dic={}
    with open(massfile,"r") as massf:            
        for line in  massf:
            if '#' in line:
                continue
            line=line.strip().split('\t')
            mass_dic[line[2]]={'C':line[1],'S':float(line[3])}
    return mass_dic

def search_key(k_index_dic,sequence):
    found=[]
    length=len(sequence)-4
    for i in range(0,length):
        kw=sequence[i:]
        node_i=0
        for item in kw:
            if k_index_dic.has_key(node_i):
                node=k_index_dic[node_i]
                if node.has_key('#'):
                    found.append(node['#'])
                if node.has_key(item):
                    node_i=node[item]
                else:
                    break
            else:
                break
    return found

def get_left(x,wstart):
    previousl=x[wstart]
    fixed_start=0
    if wstart>1:
        for i,d in enumerate(x[wstart-1::-1]):
            if d=='R' or d=='K':
                if previousl!='P':
                    fixed_start=wstart-i
                    break
            previousl=d
    return fixed_start


def get_left_seq(x,fixed_start,found_s):
    for i in range(8,40):
        new_end=fixed_start+i
        val=x[fixed_start:new_end]
        found_s[val]=1
    return found_s
        

def get_right(x,wend):
    fixed_end=wend
    nexti=wend+1
    xl=len(x)
    if nexti<xl:
        nextd=x[nexti]
        for d in x[wend:]:
            if d=='R' or d=='K':
                if nextd!='P':
                    fixed_end=nexti
                    break
            nexti=nexti+1
            if nexti<xl:
                nextd=x[nexti]
            else:
                fixed_end=nexti
                break
    return fixed_end


def get_right_seq(x,fixed_end,found_s):
    if fixed_end>len(x):
        fixed_end=len(x)
    for i in range(8,40):
        new_start=fixed_end-i
        val=x[new_start:fixed_end]
        found_s[val]=1
    return found_s

    
def get_new_seq(x,wstart,wend):
    found_s={}
    fixed_start=get_left(x,wstart)
    found_s=get_left_seq(x,fixed_start,found_s)
    
    if fixed_start>0:
        fixed_start_1=get_left(x,fixed_start-1)
        found_s=get_left_seq(x,fixed_start_1,found_s)

    fixed_end=get_right(x,wend)
    found_s=get_right_seq(x,fixed_end,found_s)

    if fixed_end<len(x):
        fixed_end_1=get_right(x,fixed_end)
        found_s=get_right_seq(x,fixed_end_1,found_s)

    return found_s
    
    

def find_match(x,y):
    found_s={}
    yl=len(y)
    for i in range(0,len(x)-yl+1):
        wstart=0+i
        wend=0+i+yl
        test=x[wstart:wend]
        if test==y:
            found_s=get_new_seq(x,wstart,wend)
    return found_s


def get_extended_search_result(seq,uid,peptide_dic,mode_detail,modific,mass_dic,alter_pep_list,outf,outf_fasta,k_index_dic,event_with_alternatives,delta):
    found=search_key(k_index_dic,seq)
    for kw_dic in found:
        kw_ev=kw_dic['e']
        tmp_pep_dic=peptide_dic[kw_ev]
        NPeptide=tmp_pep_dic['NP']
        rNPeptide=tmp_pep_dic['rN']
        pmass=tmp_pep_dic['pmass']
        #scanno=tmp_pep_dic['scan']
        #lefti=tmp_pep_dic['pos']
        key_word=kw_dic['k']
        #righti=len(NPeptide)-lefti+len(key_word)
        #newPeptide='.+?'+kw_dic['k']+'.+?'
        if pmass>1.0:
            mut_mass=pmass
        else:
            mut_mass=comput_mass(NPeptide,mass_dic)
            
        v_delta=comput_delta(mut_mass,delta)
        lookforpep = find_match(seq,key_word) 
        #continue
        pep_length=len(NPeptide)
        min_l=pep_length-3
        max_l=pep_length+3
        minl=len(seq)
        maxr=0
        for new_Ref_Ppeptide in lookforpep:
            if key_word in new_Ref_Ppeptide: 
                new_pep_l=len(new_Ref_Ppeptide)
                if new_pep_l>=min_l and new_pep_l<=max_l: 
                    new_Ref_mass=comput_mass(new_Ref_Ppeptide,mass_dic)
                    diff=mut_mass-new_Ref_mass
    
                    idiff=int(diff)
                    idiff1=idiff+1
                    idiff2=idiff-1
                    fl=0
                    if modific.has_key(idiff):
                        score=modific[idiff]
                        fl=1
                    elif modific.has_key(idiff1):
                        score=modific[idiff1]
                        fl=1
                    elif modific.has_key(idiff2):
                        score=modific[idiff2]
                        fl=1
                    if fl==1:
                        if abs(score-diff)<v_delta:
                            mitem=mode_detail[score]
                            for litem in mitem:
                                llist=litem.split(',')
                                w=llist[0]
                                if w=='*' or w in new_Ref_Ppeptide:
                                    outf.write("%s|%s|%s|" %(kw_ev,rNPeptide,key_word))
                                    outf.write("?R+MOD|%s|%f|%s\t%s\n" %(uid,score,litem,new_Ref_Ppeptide))
                                    l_b=seq.index(new_Ref_Ppeptide)
                                    r_e=l_b+new_pep_l
                                    if l_b<minl:
                                        minl=l_b
                                    if r_e>maxr:
                                        maxr=r_e
                                    break
                    elif abs(diff)<v_delta:
                        outf.write("%s|%s|%s|" %(kw_ev,rNPeptide,key_word))
                        outf.write("?R+INST_TOL|%s\t%s\n" %(uid,new_Ref_Ppeptide))
                        l_b=seq.index(new_Ref_Ppeptide)
                        r_e=l_b+new_pep_l
                        if l_b<minl:
                            minl=l_b
                        if r_e>maxr:
                            maxr=r_e
                                    
        if maxr>minl:
            minl-=1
            if minl<0:
                minl=0
            maxr+=1
            if maxr>=len(seq):
                maxr=len(seq)-1
            new_seq=seq[minl:maxr]
            if not alter_pep_list.has_key(new_seq):
                outf_fasta.write(">%s|Modified\n%s\n" %(kw_ev,new_seq))
                alter_pep_list[new_seq]=1
            for ev_it in kw_ev.split(','):
                event_with_alternatives[ev_it]=1

        
    temp_data=(event_with_alternatives, alter_pep_list)
    return temp_data


def generate_alternatives(argv):
    inputf = ''
    outfn = ''
    seqfile = ''
    keylenth = ''
    delta = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:p:k:m:",["ifile=","ofile=","pdfile=","klength=","pmass="])
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
            outfn = arg
        elif opt in ("-p", "--pdfile"):
            seqfile = arg
        elif opt in ("-k", "--klength"):
            keylenth = arg
        elif opt in ("-m", "--pmass"):
            delta = arg
    if inputf == '':
        print("Error: Wrong parameter. -i input file is missing")  
        print_help()
        sys.exit()
    
    
    if outfn == '':
        print("Error: Wrong parameter. -o output file name is missing")  
        print_help()
        sys.exit()
    
    
    if seqfile == '':
        print("Error: Wrong parameter. -p Database file is missing")  
        print_help()
        sys.exit()


    if keylenth == '':
        print("Minimum Key word length set to default 4")  
        keylenth=4
    else:
        keylenth=int(keylenth)
    
    if delta == '':
        print("Parent mass tolerance set to defualt 50ppm") 
        delta=50.0
    else:
        delta=float(delta)

    event_with_alternatives={}

    FASTA=seqfile
    outf=open(outfn,"w")
    outf_fasta=open(outfn+'.fa',"w")
    outf_report=open(outfn+'.report',"w")
    
    print("Output Files:\n1. %s\n2. %s.fa\n3. %s.report\n" %(outfn,outfn,outfn))
    temp_d_s=create_pep_table(inputf,keylenth,delta)


    peptide_dic=temp_d_s[1]
    k_index_dic=temp_d_s[0]
    events_novel=temp_d_s[2]
    
    for event in events_novel:
        outf_fasta.write(">%s|Novel\n%s\n" %(event,events_novel[event]))
    
    
    mass_dic=read_massdic("./Data/rsiduemass")
    modific={}
    mode_detail={}
    
    
    with open("./Data/moddic.txt","r") as modif:
        for line in modif:
            nline=line.strip().split('\t')
            if len(nline)==2:
                mscore=float(nline[0])
                if mode_detail.has_key(mscore):
                    mode_detail[mscore].append(nline[1])
                else:
                    mode_detail[mscore]=[nline[1]]
                ind=int(mscore)
                modific[ind]=mscore
    count=0
    alter_pep_list={}
    
    with open(FASTA,"r") as fastaf:
        for seq_record in fastaf:
            seq_record=seq_record.strip().split('\t')
            Header = seq_record[0]
            UNIPROT_ID = Header
            Protein_sequence = seq_record[1]
            
            (event_with_alternatives, alter_pep_list)=get_extended_search_result(Protein_sequence,UNIPROT_ID,peptide_dic,mode_detail,modific,mass_dic,alter_pep_list,outf,outf_fasta,k_index_dic,event_with_alternatives,delta)
    count=0
    
    for alt_pep in event_with_alternatives:        
        outf_report.write("%s," %(alt_pep))
        count+=1
    outf_report.write("\n\nNumber of events with alternatives=%d\n" %(count))


if __name__ == "__main__":
    start_time = time.time()
    generate_alternatives(sys.argv[1:])
    print("Executes in --- %s seconds ---" %(time.time() - start_time))
