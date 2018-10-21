f = open('sequence.fasta','r')
line = f.readline() #to ignore Name

chain=''
while line != '':
	line=f.readline()
	chain+=line[:len(line)-1]

#Q1-----------------------------------------------

#Function which compliments a sequence------------
def findCompliment(chain):	
	compliment=''
	for i in range(len(chain)):
		if chain[i] == 'A':
			compliment+='T'
		if chain[i] == 'T':
			compliment+='A'
		if chain[i] == 'G':
			compliment+='C'
		if chain[i] == 'C':
			compliment+='G'
	rcompliment=''.join(reversed(compliment))
	return rcompliment
#-------------------------------------------------

#To format chain to have 70 character in one line--
def format(chain):
	fchain=''
	for i in range(len(chain)):
		fchain+=chain[i]
		if (i+1)%70 == 0:
			fchain+='\n'
	fchain+='\n\n'
	return fchain
#--------------------------------------------------

o = open('1.txt','w')
rf1='> Reading Frame 1\n'
rf1+=format(chain)
o.write(rf1)

rf2='> Reading Frame 2\n'
rf2+=format(chain[1:])
o.write(rf2)

rf3='> Reading Frame 3\n'
rf3+=format(chain[2:])
o.write(rf3)

comp=findCompliment(chain)

rf4='> Reading Frame 4\n'
rf4+=format(comp)
o.write(rf4)

rf5='> Reading Frame 5\n'
rf5+=format(comp[1:])
o.write(rf5)

rf6='> Reading Frame 6\n'
rf6+=format(comp[2:])
o.write(rf6)
#-------------------------------------------------

#To Get Exonic Sequence---------------------------
fint = open('Introns.csv','r')
ex=str(chain)

line = fint.readline() #to ignore heading
line = fint.readline()
while line != '':
	lt = line.split(',')
	for i in range(int(lt[0])-1,int(lt[1])):
		ex=ex[:i]+'X'+ex[i+1:]
	line=fint.readline()

exon=''
for i in range(len(ex)):
	if ex[i] != 'X':
		exon+=ex[i]
#-------------------------------------------------

#To Get Protein from Exons------------------------
def findProtein(exon):
	codontable = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
	'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	}
	codon=[]
	flag=0
	for j in range(0,len(exon),3):
		if j+3<=len(exon):
			codon.append(exon[j:j+3])
	proteinseq=''
	for j in range(len(codon)):
		protein = codontable[codon[j]]
		if protein == '_':
			break
		if protein == 'M':
			flag=1
		if flag == 1:
			proteinseq+=protein
	return (proteinseq)
#-------------------------------------------------
maxm=0
rf=''
pf=''
pf1=findProtein(exon)
if len(pf1)>maxm:
	maxm=len(pf1)
	pf=pf1
	rf=rf1
pf2=findProtein(exon[1:])
if len(pf2)>maxm:
	maxm=len(pf2)
	pf=pf2
	rf=rf2
pf3=findProtein(exon[2:])
if len(pf3)>maxm:
	maxm=len(pf3)
	pf=pf3
	rf=rf3
exoncomp=findCompliment(exon)
pf4=findProtein(exoncomp)
if len(pf4)>maxm:
	maxm=len(pf4)
	pf=pf4
	rf=rf4
pf5=findProtein(exoncomp[1:])
if len(pf5)>maxm:
	maxm=len(pf5)
	pf=pf5
	rf=rf5
pf6=findProtein(exoncomp[2:])
if len(pf6)>maxm:
	maxm=len(pf6)
	pf=pf6
	rf=rf6

o = open('2.txt','w')
o.write(rf)	

o = open('3.txt','w')
o.write("> Protein Sequence\n")
o.write(format(pf))
#----------------------------------------------

#Q4--------------------------------------------
o = open('4.txt','w')
def gcCluster(seq):
	clen=0
	for i in range(0,len(seq)):
		if seq[i]=='G' or seq[i]=='C':
			clen+=1
		else:
			clen=0
		if clen > (0.2*len(seq)):
			return 1
	return 0

def isLongRun(seq):
	a=g=c=t=0
	for i in range(0,len(seq)):
		if seq[i]=='A':
			a+=1
			g=c=t=0
		elif seq[i]=='G':
			g+=1
			a=c=t=0
		elif seq[i]=='C':
			c+=1
			a=g=t=0
		elif seq[i]=='T':
			t+=1
			a=g=c=0
		if a>(0.2*len(seq)) or g>(0.2*len(seq)) or c > (0.2*len(seq)) or t > (0.2*len(seq)):
			return 1
	return 0

def isSelfComp(seq):
	if seq[:8]==findCompliment(seq[-8:]):
		return 1
	else:
		return 0

def primer(seq):
	fprim=[]
	ftm=[]
	for i in range(0,len(seq)-20):
		a=g=c=t=0
		for j in range(i,i+35):

			if j>=len(seq):
				break
			if seq[j] == 'A':
				a+=1
				end=1
			elif seq[j] == 'G':
				g+=1
				end=0
			elif seq[j] == 'C':
				c+=1
				end=1
			elif seq[j] == 'T':
				t+=1
				end=0

			expected=0
			seqlen=j-i+1
			if seqlen >=20 and end==1:
				expected+=2
				tm = 4*(g+c)+2*(a+t)
				gc = ((g+c)*100)/(a+g+c+t)
				if tm>=60 and tm<=70:
					expected+=1
				if gc>=40 and gc<=60:
					expected+=1
				if gcCluster(seq[i:j+1])==0:
					expected+=1
				if isLongRun(seq[i:j+1])==0:
					expected+=1
				if isSelfComp(seq[i:j+1])==0:
					expected+=1
				if expected==7:
					fprim.append(seq[i:j+1])
					ftm.append(tm)

	seqcomp=findCompliment(seq)
	flag=0
	for i in range(0,len(seqcomp)-20):
		a=g=c=t=0
		for j in range(i,i+35):
			if j>=len(seqcomp):
				break
			if seqcomp[j] == 'A':
				a+=1
				end=1
			elif seqcomp[j] == 'G':
				g+=1
				end=0
			elif seqcomp[j] == 'C':
				c+=1
				end=1
			elif seqcomp[j] == 'T':
				t+=1
				end=0

			expected=0
			seqlen=j-i+1
			if seqlen >=20 and end==1:
				expected+=2
				tm = 4*(g+c)+2*(a+t)
				gc = ((g+c)*100)/(a+g+c+t)
				if tm>=60 and tm<=70:
					expected+=1
				if gc>=40 and gc<=60:
					expected+=1
				if gcCluster(seqcomp[i:j+1])==0:
					expected+=1
				if isLongRun(seqcomp[i:j+1])==0:
					expected+=1
				if isSelfComp(seqcomp[i:j+1])==0:
					expected+=1
				if expected==7:
					k=0
					for prim in fprim:
						if abs(ftm[k]-tm)<=5:
							o.write("Forward Primer - "+prim+"\tTm(in deg C) - "+str(ftm[k])+"\n")
							o.write("Reverse Primer - "+seqcomp[i:j+1]+"\tTm(in deg C) - "+str(tm)+"\n\n")
							flag=1
							break
						k+=1
			if flag==1:
				break
		if flag==1:
			break
primer(chain)
primer(chain[1:])
primer(chain[2:])
#Q5--------------------------------------------
o = open('5.txt','w')
o.write('> E.coli cloning sequence \n')
o.write(format(exon))

rf1=rf1[rf1.find('\n')+1:]
o.write('> Pichia pastoris cloning sequence \n')
o.write(rf1)

o.write('> HEK293 cloning sequence \n')
o.write(rf1)
#----------------------------------------------
