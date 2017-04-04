FASTA='../uniprot_sprot.tab'
UniProt={}
with open(FASTA,"r") as fastaf:
    for seq_record in fastaf:
        seq_record=seq_record.strip().split('\t')
        if len(seq_record)>=2:
            Header = seq_record[0]
            if "_HUMAN" in Header:
                UNIPROT_ID = Header.split('|')[1]
                Gene_Name = ''    
                if Header.find('GN=') >-1:
                    Gene_Name= Header.split('GN=')[1].split(' ')[0]
                Protein_sequence = seq_record[1]+'#'
                UniProt[UNIPROT_ID]={'seq':Protein_sequence,'gene':Gene_Name}

outf=open("alternate.fasta80.tab","w")

for uid in UniProt:
    seq=UniProt[uid]['seq']
    flag_f=0
    first_s=0
    for i,w in enumerate(seq):
        if i!='#':
            if w =='R' or w=='K' and seq[i+1]!='P':
                full_word=seq[first_s:i+1]
                sl=len(full_word)
                if  sl<40:
                    second_s=i
                if sl>80:
                    outf.write("%s\t%s\n" %(uid,full_word))
                    first_s=second_s+1
                    second_s=i
    full_word=seq[first_s:-1]
    if len(full_word)>5:
        outf.write("%s\t%s\n" %(uid,full_word))
            

outf.close()
