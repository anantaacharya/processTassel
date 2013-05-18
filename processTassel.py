#Coded by Ananta Acharya
#contact ananta@uga.edu if any questions
#This can be used as is or you can modify, The results are not guaranted to be correct, use with caution


import os
import commontools



def main():

    print "build hapmap with minimum coverage and genotype correction """
    par=params()

    print """parameters
        HapMap.hmp=%s
        HapMap.hmc=%s
        out=%s
        mindepth_homo=%s
        mindepth_het=%s
        minratio_het=%s
        minratio_homo=%s
        mindepth_minor=%s
    """ %par
    if "Y" in raw_input("continue?").upper():
        rehapmap(*par)
    


        print "Finished processing"
    else:
        pass
        
def params():    

    hap=raw_input("Hapmap File? Hapmap.hmp.txt").strip().strip('"')
    while not commontools.OKtoRead(hap):
        hap=raw_input("Hapmap File? Hapmap.hmp.txt").strip().strip('"')
    hapcount=raw_input("Hapmap Count File, Hapmap.hmc.txt").strip().strip('"')
    while not commontools.OKtoRead(hapcount):
        hapcount=raw_input("Hapmap Count File, Hapmap.hmc.txt").strip().strip('"')
    mindepth_homo=raw_input("minimum depth to call homozygous for each genotype")
    while not commontools.isNumber(mindepth_homo):
        mindepth_homo=raw_input("minimum depth to call homozygous for each genotype")
    mindepth_het=raw_input("minimum depth to call heterozygous for each genotype")
    while not commontools.isNumber(mindepth_het):
        mindepth_het=raw_input("minimum depth to call heterozygous for each genotype")
    minratio_het=raw_input("cutoff ratio of alleles to call heterozygous ie 0.1, if ratio higher than this will be heterozygous")
    while not commontools.checkRatio(minratio_het):
        minratio_het=raw_input("cutoff ratio of alleles to call heterozygous ie 0.1, if ratio higher than this will be heterozygous")
    minratio_homo=raw_input("cutoff ratio of alleles to call homozygous ie 0.1, if ration lower than this will be homozygous")
    while not commontools.checkRatio(minratio_homo):
        minratio_homo=raw_input("cutoff ratio of alleles to call homozygous ie 0.1, if ration lower than this will be homozygous")
    minMinor=raw_input("minimum read depth of minor allele to call heterozygous")
    while not commontools.isNumber(minMinor):
        minMinor=raw_input("minimum read depth of minor allele to call heterozygous")
    hapout=raw_input("Where to write rehapmap file")
    while not commontools.OKtoWrite(hapout):
        hapout=raw_input("Where to write rehapmap file")
    return (hap, hapcount,hapout,  int(mindepth_homo),int(mindepth_het), float(minratio_het), float(minratio_homo),int(minMinor))
            
    
def rehapmap(hapmap,hapcount, outf, mindepth_homo,mindepth_het,minratio_het, minratio_homo,minMinor):
    
    
    pairdict=dict()
    hapfile=open(hapmap,"r")
    #hapcountf=open(hapcount,"r")
    rehap=open(outf,"a")
    lessreadCount=0
    minMinorCount=0

    for line in hapfile:
        linelist=line.strip().split("\t")
        pairdict[linelist[0]]=linelist[1]
    

    with open(hapcount) as hapcountf:
        line1=next(hapcountf)
        writebuffer="\t".join(line1.strip().split("\t")[0:-5])
        writebuffer+="\n"

        for line in hapcountf:
            linebuffer=""
            linewrite=False
            linelist=line.strip().split("\t")
            ID=linelist[0]
            cntl=[each.split("|") for each in linelist[1:-5]]
            sumeach=0
            for each in cntl:
                sumeach+=sum([int(ee) for ee in each])
            if sumeach<>0:
                #print "sumeach not zero", ID
                linebuffer=""
                linewrite=False
                linebuffer+=ID
                for x in range(1,len(linelist)-5):
                   # raw_input()
                    snps=pairdict[ID]

                    snps=snps.split("/")
                    counts=linelist[x].split("|")

                    allele1=int(counts[0])
                    allele2=int(counts[1])
                    if allele1+ allele2>0:
                        #print "sum  not zero", ID, str(x)
                        if allele1==0:
                            
                            if allele2>mindepth_homo:
                                #print "allele1 zero, allele2 not zero"
                                linebuffer+="\t"+snps[1]
                                linewrite=True
                            else:
                             #   print "allele1 zero, allele2 less than min"
                                linebuffer+="\tN"
                        elif allele2==0:
                            if allele1>mindepth_homo:
                                #print "allele1 zero, allele2 not zero"
                                linebuffer+="\t"+snps[0]
                                linewrite=True
                            else:
                              #  #print "allele1 zero, allele2 less than min"
                                linebuffer+="\tN"
                                                            
                        elif (allele1 + allele2)> mindepth_het:
                                                     
                                                    
                            #print "sum  not zero", ID, str(x), str(allele1), str(allele2), ">min het"   
                                                    
                            if min(allele1,allele2)<minMinor:
                                #print "less than minminor"
                                linebuffer+="\tN"
                                
                            elif minratio_het<=(allele1/float(allele1+allele2))<=(1-minratio_het):
                                #print "sum  not zero", ID, str(x), str(allele1), str(allele2), ">min het its het"  
                                linebuffer+="\t"+iupac2(snps)
                                linewrite=True
                            else:
                                minorallele=min(allele1, allele2)
                                if minratio_homo>=(minorallele/float(allele1+allele2)):
                                    #print "possible homo"

                                    
                                    if minorallele==allele1:
                                        minorsnp=snps[0]
                                        majorsnp=snps[1]
                                    else:
                                        minorsnp=snps[1]
                                        majorsnp=snps[0]
                                    linebuffer+="\t"+majorsnp
                                    linewrite=True
                                
                                
                                
                                else:
                                    #print "in between homo and hetero"            
                                        
                                    linebuffer+="\tN"
                                    lessreadCount+=1
                        else:
                            #print "removed not enough depth"
                            linebuffer+="\tN"
                            lessreadCount+=1
                    else:
                        ##print "zero sum"
                        linebuffer+="\tN"
                        lessreadCount+=1
            if linewrite:
                writebuffer+=linebuffer
                writebuffer+="\n"
            
            
    rehap.write(writebuffer)





     
        
def iupac2(listt):

    

    if "A" in listt and "T" in listt:
        return "W"
    elif "A" in listt and "C" in listt:
        return "M"
    elif "A" in listt and "G" in listt:
        return "R"
    elif "C" in listt and "T" in listt:
        return "Y"
    elif "C" in listt and "G" in listt:
        return "S"
    elif "G" in listt and "T" in listt:
        return "K"
    else:
        return "N"
    

main()
        
