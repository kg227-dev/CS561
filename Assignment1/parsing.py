import sys,re
def fasta(file):
    f=open(file,'r')
    lines=f.readlines()
    hre=re.compile('>(\S+)')
    lre=re.compile('^(\S+)$')
    gene={}
    for line in lines:
            outh = hre.search(line)
            if outh:
                    id=outh.group(1)
            else:
                    outl=lre.search(line)
                    if(id in gene.keys()):
                            gene[id] += outl.group(1)
                    else:
                            gene[id]=outl.group(1)
    return gene

x = fasta('close-first.fasta').get('first1')
y = fasta('close-second.fasta').get('second1')
print(x)
print(y)
