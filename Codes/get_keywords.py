import sys, getopt
import heapq
import time



def print_help():
    print 'python ./Codes/get_keywords.py -i <input PSMs file (From exctract_spec.py)> -s <Single spectra file (From exctract_spec.py)> -o <outputfile>'

def getannotdictionary(afile):
    with open(afile,"r") as af:
        d={}
        for line in af:
            if '#' in line:
                continue
            line=line.strip().split("\t")
            d[line[2]]=line[3]
    return d	


def getb(name,val,charge):
	s=len(name)-1
	tlist=[]
	for i in range(0,s):
		l=name[i]
		if l in ms_dict.keys():
			oldval=val
			val=float(ms_dict[l])+oldval
			if l=="C":
				val=val+57
			m=(val+charge)/charge
			tlist.append(m)
			#print("%s\t%f\n" %(l,tlist[i]))
	return tlist


def getpmasstotal(name,ms_dict):
    val=19
    for l in name:
        if l in ms_dict.keys():
            val=float(ms_dict[l])+val
            if l=="C":
                val=val+57
    return val


def nth_largest(n, iter):
    return heapq.nlargest(n, iter)[-1]
    
    
ms_dict=getannotdictionary("./Data/rsiduemass")

def get_position_mass(pep,scan,pm): 
    position={}
    rpep=pep[::-1]
    s=len(pep)
    additional_mass_at_n_terminal=0.0
    count=0
    charge=1
    blist=getb(pep,additional_mass_at_n_terminal,charge)
    ylist=getb(rpep,18.0,charge)
    ylist=ylist[::-1]
    tcount=0
    position_mass_b={}
    position_mass_y={}
    position_mass_t={}
    position_mass_u={}

    for i in range(0,s):
        position_mass_b[i]=0
        position_mass_y[i]=0
        position_mass_t[i]=0
        position_mass_u[i]=0
        if i==0:
            bi=int(blist[i])
            if scan.has_key(bi):
                bival=scan[bi]['x']
                inv=abs(float(bival)-blist[i])
                if inv<pm:
                    position_mass_b[i]=1
                    position_mass_u[i]=1
            elif scan.has_key(bi+1):
                bival=scan[bi+1]['x']
                inv=abs(float(bival)-blist[i])
                if inv<pm:
                    position_mass_b[i]=1
                    position_mass_u[i]=1
        elif i==(s-1):
            yi=int(ylist[i-1])
            if scan.has_key(yi):
                yival=scan[yi]['x']
                inv=abs(float(yival)-ylist[i-1])
                if inv<pm:
                    position_mass_y[i]=1
                    position_mass_u[i]=1                    
            elif scan.has_key(yi+1):
                yival=scan[yi+1]['x']
                inv=abs(float(yival)-ylist[i-1])
                if inv<pm:
                    position_mass_y[i]=1                    
                    position_mass_u[i]=1
                    
        else:
            bi=int(blist[i])
            if scan.has_key(bi):
                bival=scan[bi]['x']
                inv=abs(float(bival)-blist[i])
                if inv<pm:# and position_mass_b[i-1]==1:
                    position_mass_b[i]=1
                    position_mass_u[i]=1
            elif scan.has_key(bi+1):
                bival=scan[bi+1]['x']
                inv=abs(float(bival)-blist[i])
                if inv<pm:# and position_mass_b[i-1]==1:
                    position_mass_b[i]=1
                    position_mass_u[i]=1
            yi=int(ylist[i-1])
            if scan.has_key(yi):
                yival=scan[yi]['x']
                inv=abs(float(yival)-ylist[i-1])
                if inv<pm: # and position_mass_y[i-1]==1:
                    position_mass_y[i]=1
                    position_mass_u[i]=1
            elif scan.has_key(yi+1):
                yival=scan[yi+1]['x']
                inv=abs(float(yival)-ylist[i-1])
                if inv<pm: # and position_mass_y[i-1]==1:
                    position_mass_y[i]=1
                    position_mass_u[i]=1
                 
    for i in range(0,s-1):
        val=position_mass_b[i]
        if val == 1:
            inext=i+1
            if position_mass_b.has_key(inext):
                nval=position_mass_b[inext]
                if nval==1:
                    tcount+=1
                    position_mass_t[inext]=1
                    
    for i in range(0,s-1):
        if position_mass_t[i]==1:
            continue
        val=position_mass_y[i]
        if val == 1:
            inext=i+1
            if position_mass_y.has_key(inext):
                nval=position_mass_y[inext]
                if nval==1:
                    tcount+=1
                    position_mass_t[i]=1
            
    count=tcount
    for k in position_mass_t:
        position[k]=position_mass_t[k]


    pitem=''
    uitem=''
    for k in position:
        v=position[k]
        u=position_mass_u[k]
        pitem+=str(v)
        uitem+=str(u)
        
    
    result=(pitem,count,uitem)
    return result
    


def get_key(argv):
    #ignore equivalent mass residues
    #ignore  phosphorylation residue serine
    replace_d={'GG':2,'GA':2,'GV':2,'GE':2,'AD':2,'SV':2,'SS':2,'W':1,'N':1,'Q':1,'K':1,'R':1}
    inputfile = ''
    scanfile = ''
    outputfile = ''

    try:
        opts, args = getopt.getopt(argv,"hi:o:s:",["ifile=","ofile=","sfile="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print_help()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-s", "--sfile"):
            scanfile = arg

    if inputfile == '':
        print("Error: Wrong parameter. -i input file is missing")  
        print_help()
        sys.exit()
    
    
    if outputfile == '':
        print("Error: Wrong parameter. -o output file name is missing")  
        print_help()
        sys.exit()

    if scanfile == '':
        print("Error: Wrong parameter. -s Spectra files is missing")  
        print_help()
        sys.exit()

    delta=0.5
    mapscan={}
    eventdic={}
    peptide={}
   
    with open(inputfile,"r") as mf:
        for line in mf:
            if '#' in line:
                continue
            if len(line)<1:
                break
            line=line.strip().split("\t")
            mapscan[int(line[1])]=line[0]
            if '.' in line[3]:
                eventdic[line[0]]=line[3][1:-1].translate(None, '1234567890:._-*!@#$?')
            else:
                eventdic[line[0]]=line[3].translate(None, '1234567890:._-*!@#$?')
                
            peptide[eventdic[line[0]]]=line[3]
       

    scan_i=-1
    eventsupport={}

    with open(scanfile,"r") as sf:
        matched=False
        for line in sf:
           if '=' in line:
               continue
           if 'BEGIN IONS' in line:
               scan_i+=1
               scand={}
               slist=[]
               stemp={}
               scount=0
               if mapscan.has_key(scan_i):
                   matched=True
                   event=mapscan[scan_i]
                   pep=eventdic[event]
               else:
                   matched=False
            
           elif 'END IONS' not in line and matched==True:
               line=line.strip().split()
               if len(line)>=2:
                   xval=float(line[0])
                   yval=float(line[1])
                   stemp[scount]={'y':yval,'x':xval}
                   scount+=1
                   slist.append(yval)
           elif 'END IONS' in line and matched==True:
               l50=nth_largest(50, slist)
               for i in stemp:
                   yval=stemp[i]['y']
                   xval=stemp[i]['x']
                   if yval>=l50:
                       li=int(xval)
                       if scand.has_key(li):
                           oyval=scand[li]['y']
                           if oyval<yval:
                               scand[li]={'y':yval,'x':xval}
                       else:
                           scand[li]={'y':yval,'x':xval}
               pm=delta
               pd=get_position_mass(pep,scand,pm)
               
               if eventsupport.has_key(event):
                   eventsupport[event].append({'count':pd[1],'s':pep,'su':pd[0],'sno':scan_i,'ubit':pd[2]})
               else:
                   eventsupport[event]=[{'count':pd[1],'s':pep,'su':pd[0],'sno':scan_i,'ubit':pd[2]}]
                   
    outf=open(outputfile,"w")
    
    outf.write("#Event_ID\tInput_Peptide(I)\tInput_Peptide(II)\tParent_Mass\tSupport_bits\tKey_Words\tSupport_Count\tSupport_Score\tScan_Number\tb_and_y_ions_match\n")
    
    for event in eventsupport:
        for event_dic in eventsupport[event]:
            count=event_dic['count']
            if count==0:
                continue
            pep=event_dic['s']
            support=event_dic['su']
            scan_number=event_dic['sno']
            
            ubit=event_dic['ubit']
            
            printsupport=support
            newpep=''
            support=support+'#'
            total=len(pep)
            '''
            for i,item in enumerate(pep):
                if support[i]=='1':
                    newpep+=item
                else:
                    newpep+='.'
 
            '''
            previous=''
            for i,item in enumerate(pep):
                k_item=previous+item
                
                if not replace_d.has_key(k_item) and not replace_d.has_key(item):
                    if support[i]=='1':
                        newpep+=item
                    else:
                        newpep+='.'
                else: 
                    if support[i]=='1':
                        count-=1
                        if i>0 and replace_d.has_key(k_item):
                            if support[i-1]==1:
                                newpep=newpep[:-1]+'.'
                                count-=1
                    newpep+='.'
                previous=item
            if count<0:
                count=0
            #ends here
                
            fra=100.0*(float(count)/total)
    
            outf.write("%s\t%s\t%s\t%f\t%s\t%s\t%d\t%0.2f\t%d\t%s\n" %(event,peptide[pep],pep,getpmasstotal(pep,ms_dict),printsupport,newpep,count,fra,scan_number,ubit))
    outf.close()

if __name__ == "__main__":
    start_time = time.time()
    get_key(sys.argv[1:])
    print("Executes in --- %s seconds ---" %(time.time() - start_time))
