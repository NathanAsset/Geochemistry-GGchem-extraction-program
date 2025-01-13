### Library needed for the program

import matplotlib.pyplot as plt         #For the graphes
import numpy as np
import csv                              #For the creation of the files

### The different functions
# function looking at the name of a element in the first line of the GGchem exit file, and extract its place in the line

def rech(code):                        #The entry of the program is the element look for
    key=-1                             #Original exit values that cannot be a position on the line
    for j in range (0,len(keyword)):   #The search part througt the line
        if keyword[j]==code:
            key=j
    return(key)

#Function producing a list of lists composed by molar fraction for specifics element (la première est pour les solides, la deuxième pour les gaz, le troisième est pour toutes les valeurs pour faire une somme)

def creatl(key,info,gas,tot,lH):
    if info==True and gas==False:           #For the solids
        lis=[]
        for i in range (0, len(dat)):
            d=dat[i][key]
            df=((10**d)*lH[i])/(6.02E23)    #Correction to pass from log10(ncond/nHtot) to molar composition of the condensates
            dfi=df/tot[i]                   #Calculated the molar fraction
            if df<(10**-23):                #Limitation of solids molar fraction
                dfi=0
            lis.append(dfi)
    elif info==True and gas==True:          #For the gas
        lis=[]
        for i in range (0, len(dat)):
            d=dat[i][key]
            df=((10**d)*(9.2E13))/(6.02E23) #Correction to pass from log10(natom)x[cm-3] to molar composition of the gas
            dfi=df/tot[i]                   #Calculated the molar fraction
            if dfi<(10**-10):               #Limitation of gas molar fraction
                dfi=0
            lis.append(dfi)
    elif info=='H' and gas==True:
        lis=[]
        for i in dat:
            d=i[key]
            dH=d*(9.2E13)                   #Correction to pass Hydrogen nuclei particle density from to atomic content of total H
            lis.append(dH)
    else:                                   #All elements for the sum and to create the list of temperature
        lis=[]
        for i in dat:
            d=i[key]
            lis.append(d)
    return(lis)

#Function calculating the sum for all elements (condensates and gas)

def somme(som,lH):
    tot=[]
    for i in dat:
        t=0
        j=0
        for so in som:
            key=rech(so)                        #Search of all elements (gas/condensats)
            d=i[key]                            #Values of the elements
            if so[0]=='n':                      #For the condensates
                df=((10**d)*lH[j])/(6.02E23)    #Correction to pass from log10(ncond/nHtot) to molar composition of the condensates
            else:                               #For the gas
                df=((10**d)*(9.2E13))/(6.02E23) #Correction to pass from log10(natom)x[cm-3] to molar composition of the gas
            t=t+df
        j=j+1
        tot.append(t)
    return(tot)

#Function producing the list in molar fraction for the gas bearing specific atoms
def lprod(l,tot2,lH,lgf):
    for gas in l:                               #List of the different gas
        keys=rech(gas)                          #Find the position of the data in the GGchem data for one gas
        lgcr=creatl(keys,True,True,tot2,lH)     #Compalling the data for one gas producing a list
        lgf.append(lgcr)                        #Add the list of one gas to the final list
    return(lgf)

#Function calculating the sum of molar fraction
def sfu(T,lg):                                       of all gas,for each temperature
    lgt=[]
    for i in range (0,len(T)):      #For each temperature
        so=0
        for j in lg:                #For each gas
            so=so+j[i]              #Creation the sum
        lgt.append(so)
    return(lgt)

#Creation and save data file
def sauvfil(name,lng,lg,T):
    fichsauv=open('C:/Users/natha/Desktop/'+name+'.csv','w')    #Opening of a file and definition of its name
    prli="T;"
    for sauvg in range (0,len(lng)):
        prli=prli+str(lng[sauvg])+";"
        prli=prli+"\n"                          #Writing the frist line with the temperature and the name of the gas with oxygen
        fichsauv.write(prli)
    for ligne in range (0,len(T)):          #Writing the lines for each temperature,of the Partial Pressure of the gas with oxygen
        CvT=str(T[ligne])+';'
        val=CvT
        for ligng in range(0,len(lg)-1):
            Cvg=str(lgf[ligng][ligne])+';'
            val=val+Cvg
        Cvgf=str(lgf[len(lg)-1][ligne])+'\n'
        val=val+Cvgf
        fichsauv.write(val)
    fichsauv.close()


###Opening the GGchem file Static-conc.dat

file   = r'C:\Users\natha\Desktop\GGchem\Static_Conc.dat'  #The path to the file with the results of the GGchem code
data   = open(file)                     #Remnant of the GGchem plot code for oppening the file (necessary, but only use here)
dummy  = data.readline()
dimens = data.readline()
dimens = np.array(dimens.split())
NELEM  = int(dimens[0])
NMOLE  = int(dimens[1])
NDUST  = int(dimens[2])
NPOINT = int(dimens[3])
header = data.readline()
data.close()
dat = np.loadtxt(file,skiprows=3)       #Creation of a list of lists with all the data of the condensation simulation
keyword = np.array(header.split())      #Creation of a list with the name of all chemical element (condensates and gas) in the condensation modelisation

###Part that etablish the list of elements from the file Static_Conc.dat and the conversion in molar fraction

lg=[]
lSig=[]
lSite=[]
lOg=[]
lNg=[]
lAlg=[]
lCag=[]
solid=[]
som=[]
keyT=rech('Tg')                  #Creation of the Temperature list
T=creatl(keyT,False,False,[],[])
keyH=rech('nHges')               #Creation of the list of Hydrogen atomic content
lH=creatl(keyH,'H',True,[],[])
for s in keyword:                #Creation of the list with all the solids produce during the condensation
    if s[0]=='n' and s!='nHges':
        solid.append(s)
for so in keyword:               #Creation of the list with all the elements (gas and condensates) present during the condensation
    if so[0]!='S' and so[0]!='e' and so[0]!='d' and so[len(so)-2]!='W' and so!='Tg' and so!='nHges' and so!='pgas': #Excluded code to limit some column of data
        som.append(so)
    if so[0]=='S' and so[1]=='i' or so[0]=='S' and so[1]=='I':  #Correction from the exclusion of the precedent line
        som.append(so)

for g in keyword:                #Creation of the list with all the gas present during the condensation
    if g[0]!='e' and g[0]!='d' and g[0]!='n' and g[len(g)-2]!='W' and g!='Tg' and g!='nHges' and g!='pgas':
        if g[0]!='S':
            lg.append(g)
        if g[0]=='S' and g[1]=='I' or g[0]=='S' and g[1]=='i':
            lg.append(g)


for Sig in keyword:              #Creation of the list with all the gas with silicium, present during the condensation
    for nSi in range(0,len(Sig)):
        if Sig[nSi]=='S' and Sig[0]!='n' and Sig[0]!='e' and Sig[0]!='S':
            if Sig[nSi+1]=="i" or Sig[nSi+1]=="I":
                lSig.append(Sig)
    if Sig[0]=='S':
        if Sig[1]=="i" or Sig[1]=="I":
            lSig.append(Sig)

for Og in keyword:               #Creation of the list with all the gas with oxygen, present during the condensation
    for o in Og:
        if o=='O' and Og[0]!='S' and Og[0]!='n' and Og[0]!='e':
            lOg.append(Og)
        if o=='O' and Og[0]!='n' and Og[0]!='e':
            if Og[0]=='S' and Og[1]=='I':
                lOg.append(Og)

for Ng in keyword:               #Creation of the list with all the gas with nitrogen, present during the condensation
    for n in range(-1, len(Ng)-1):
        if Ng[n]=='N' and Ng[n+1]!="a" and Ng[n+1]!="A" and Ng[0]!='S' and Ng[n]!='W' and Ng[0]!='n' and Ng[0]!='e':
            lNg.append(Ng)

for Alg in keyword:               #Creation of the list with all the gas with aluminium, present during the condensation
    for Al in range(-1, len(Alg)-1):
        if Alg[Al]=='A' and Alg[Al+1]=="L" and Alg[0]!='S' and Alg[Al]!='W' and Alg[0]!='n' and Alg[0]!='e':
            lAlg.append(Alg)

for Cag in keyword:               #Creation of the list with all the gas with calcium, present during the condensation
    for Ca in range(-1, len(Cag)-1):
        if Cag[Ca]=='C'and Cag[0]!='S' and Cag[Ca]!='W' and Cag[0]!='n' and Cag[0]!='e':
            if Cag[Ca+1]=="a" or Cag[Ca+1]=="A":
                lCag.append(Cag)

#Print of the lenght of the list of name, follow by the list of name of gas and condensates

print(len(solid),solid)     #List of condensates
print(len(lSig),lSig)       #list with all the gas with silicium
print(len(lOg),lOg)         #list with all the gas with oxygen
print(len(lNg),lNg)         #list with all the gas with nitrogen
print(len(lAlg),lAlg)       #list with all the gas with aluminnium
print(len(lCag),lCag)       #list with all the gas with calcium

#Lines to produced a list of data (molar fraction) with the data for different system
tot2=somme(som,lH)
lgf=[]                                      #Final list section
lSigf=[]
lOgf=[]
lNgf=[]
lAlgf=[]
lCagf=[]
lgf=lprod(lg,tot2,lH,lgf)                   #List of data for all the gases
lSigf=lprod(lSig,tot2,lH,lSigf)             #List of data with all the gases with silicon
lOgf=lprod(lOg,tot2,lH,lOgf)                #List of data with all the gases with oxygen
lNgf=lprod(lNg,tot2,lH,lNgf)                #List of data with all the gases with nitrogen
lAlgf=lprod(lAlg,tot2,lH,lAlgf)             #List of data with all the gases with aluminium
lCagf=lprod(lCag,tot2,lH,lCagf)             #List of data with all the gases with calcium


lgt=sfu(T,lgf)                              # Calculating the sum of molar fraction of all gas,for each temperature
lSit=sfu(T,lSigf)                           # Calculating the sum of molar fraction of all gas with Si,for each temperature
lOt=sfu(T,lOgf)                             # Calculating the sum of molar fraction of all gas with O,for each temperature


#Lines to produced a list of data (molar fraction) with the data with all tha gas with oxygen, present during the condensation
lsolid=[]
ls=[]
for sol in solid:
    keys=rech(sol)
    lsol=creatl(keys,True,False,tot2,lH)#Creation of the data list for all solid
    if sum(lsol)!=0:
        soli=sol.replace('n',' ')       #Rewrite the name of the condensate by deleted n
        ls.append(soli)
        lsolid.append(lsol)

### Initial temperature, 50% temperature, Final temperature of condensation
#Definition of the initial temperature of condensation, last temperature of condensation, and the 50% condensation temperature
T50=[]
Tini=[]
Tfin=[]
for i in range (0,len(lsolid)):     #Looking for the three temperatures for each condensates
    Tf=6000
    jf=0
    m50=(lsolid[i][len(lsolid[i])-1])/2
    for j50 in range (1,len(lsolid[i])-1):
        if lsolid[i][j50-1]<=m50 and lsolid[i][j50+1]>=m50: #The temperature at 50% condensation definition
            T1=T[j50-1]
            T2=T[j50+1]
            T50val=(T1+T2)/2
            T50.append(ls[i])
            T50.append(T50val)
        if lsolid[i][j50-1]==0 and lsolid[i][j50]!=0:       #The initial temperature of condensation definition
            Tini.append(ls[i])
            Tini.append(T[j50])
            Tini.append(lsolid[i][j50])
        if lsolid[i][j50+1]-lsolid[i][j50]>1E-9 and lsolid[i][j50]!=0 and Tf>lsolid[i][j50]:                                                      #The final temperature of condensation definition
            Tf=T[j50]
            jf=j50
    Tfin.append(ls[i])
    Tfin.append(Tf)
    Tfin.append(lsolid[i][jf])

### Calcul of the fO and fSi
# Creation of the lists of all elements (condensat and gas) with Si or with O
lSigc=[]
lSigcf=[]
lOgc=[]
lOgcf=[]
for gSi in range(0,len(lSig)):      # Add the gas with Si to a list of name and the data for all element with Si
    lSigc.append(lSig[gSi])
    lSigcf.append(lSigf[gSi])
for cSi in range(0,len(lsolid)):    # Add the condensates with Si to a list of name and the data for all element with Si
    for ncSi in range(0,len(ls[cSi])-1):
        if ls[cSi][ncSi]=='S' and ls[cSi][ncSi+1]=='i':
            lSigc.append(ls[cSi])
            lSigcf.append(lsolid[cSi])
for gO in range(0,len(lOg)):        # Add the gas with O to a list of name and the data for all element with O
    lOgc.append(lOg[gO])
    lOgcf.append(lOgf[gO])
for cO in range(0,len(lsolid)):      # Add the condensates with O to a list of name and the data for all element with O
    for ncO in range(0,len(ls[cO])-1):
        if ls[cO][ncO]=='O':
            lOgc.append(ls[cO])
            lOgcf.append(lsolid[cO])

# Calcul of the quantity of Si and O in the different gas and condensates
#Dictionaries with the fraction of Si or O in a given molecule
dfO={'O':1, 'O2':1, 'OH':0.9407303, 'NO':0.53331987, 'ALO':0.3722342, 'SIO':0.3629167, 'CAO':0.2853041, 'ALOH':0.3637046, 'HALO':0.3637046, 'ALO2H':0.5334067, 'ALO2':0.5425229, 'AL2O':0.228678, 'AL2O2':0.3722342, 'CAOH':0.2802663, 'CA(OH)2':0.4318685, 'HNO':0.5158688, 'HONO':0.6806246, 'HNO2':0.6806246, 'HNO3':0.7617157, 'HO2':0.9694601, 'H2O':0.8880933, 'NO2':0.6955376, 'NO3':0.774099, 'N2O':0.3635112, 'N2O3':0.6314531, 'N2O4':0.6955376, 'N2O5':0.7406368, 'SIO2':0.5325589, 'O3':1,' SiO2[l]':0.5325589,' SiO2':0.5325589,' Al2O3[l]':0.4707388,' Al2O3':0.4707388,' CaO[l]':0.2853041,' CaO':0.2853041,' CaSiO3':0.6308547,' Ca2Al2SiO7':0.4084377,' Ca2SiO4':0.3715567,' CaAl2Si2O8':0.4600636,' Ca3Al2Si3O12':0.426219463}
dfSi={'Si':1, 'SI2':1, 'SIH':0.9653531, 'SIN':0.6672376, 'SIO':0.6370833, 'SIH4':0.874461, 'SI2N':0.800411, 'SIO2':0.4674411, 'SI3':1, 'SIH2':0.9330266, 'SIH3':0.902795,' SiO2[l]':0.4674411,' SiO2':0.4674411,' CaSiO3':0.3691453,' Ca2Al2SiO7':0.1024276,' Ca2SiO4':0.3715567,' CaAl2Si2O8':0.2019051,' Ca3Al2Si3O12':0.187055439}

lSie=[]
lOe=[]
neleSi=0
for Sie in lSigcf: #Calcul the Si only composition in all element for all temperature
    lSiein=[]
    for a in range (0, len(Sie)):
        fSi=(Sie[a]*dfSi[lSigc[neleSi]])    #Calcul the molar composition in Si in a molecule with Si
        lSiein.append(fSi)
    lSie.append(lSiein)
    neleSi=neleSi+1

neleO=0
for Oin in lOgcf:   #Calcul the O only composition in all element for all temperature
    lOin=[]
    for a in range (0, len(Oin)):
        fO=(Oin[a]*dfO[lOgc[neleO]])         #Calcul the molar composition in O in a molecule with O
        lOin.append(fO)
    lOe.append(lOin)
    neleO=neleO+1

lSitot=[]
lOtot=[]
for ntot in range(0, len(lSie[0])):
    Otot=0
    Sitot=0
    for Siet in lSie:
        Sitot=Siet[ntot]+Sitot #Calcul the sum of the composition of Si for all element for each temperature
    for Oet in lOe:
        Otot=Otot+Oet[ntot]    #Calcul the sum of the composition of O for all element for each temperature
    lSitot.append(Sitot)
    lOtot.append(Otot)

# Estimate of the maximum content of Si
max=0       #Value of the maximum of Si molar composition
nTSi=0      #Position of the maximum of Si molar composition in the list of the Si sum
for nSit in range(0,len(lSitot)):
    if lSitot[nSit]>=max:
        max=lSitot[nSit]
        nTSi=nSit

# Estimate of the maximum content of O
maxO=0      #Value of the maximum of O molar composition
nTO=0       #Position of the maximum of O molar composition in the list of the O sum
for nOt in range(0,len(lOtot)):
    if lOtot[nOt]>=maxO:
        maxO=lOtot[nOt]
        nTO=nOt

# first step to calculate the fSi anf the fO
lfSie=[]
lfOe=[]
neleSi=0
for Sie in lSigf:   #Calcul the Si only composition in all gas with Si for all temperature
    lSiein=[]
    for a in range (0, len(Sie)):
        fSi=(Sie[a]*dfSi[lSig[neleSi]]) #Calcul the molar composition in Si in a molecule of gas with Si
        lSiein.append(fSi)
    lfSie.append(lSiein)
    neleSi=neleSi+1

neleO=0
for Oin in lOgf:    #Calcul the O only composition in all gas with O for all temperature
    lOin=[]
    for a in range (0, len(Oin)):
        fO=(Oin[a]*dfO[lOg[neleO]]) #Calcul the molar composition in O in a molecule of gas with O
        lOin.append(fO)
    lfOe.append(lOin)
    neleO=neleO+1

lfSiel=[]
lfOel=[]

# Calcul to estimate the fSi for all gas with Si for each temperature
for Sie in lfSie:
    lfSiein=[]
    for nSi in range(0, len(T)):
        if lSitot[nSi]!=0:
            fSi=Sie[nSi]/(lSitot[nSi])
            lfSiein.append(fSi)
        else:
            lfSiein.append(0)
    lfSiel.append(lfSiein)

# Calcul to estimate the fO for all gas with O for each temperature
for Oe in lfOe:
    lfOein=[]
    for nO in range(0, len(T)):
        if lOtot[nO]!=0:
            fO=Oe[nO]/(lOtot[nO])
            lfOein.append(fO)
        else:
            lfOein.append(0)
    lfOel.append(lfOein)

#Sum of the fSi of each gas with Si, for all temperature
lfOtot=[]
lfSitot=[]
for i in range(0,len(T)):
    sfSi=0
    for fSi in lfSiel:
        sfSi=sfSi+fSi[i]
    lfSitot.append(sfSi)

#Sum of the fO of each gas with O, for all temperature
for i in range(0,len(T)):
    sfO=0
    for fO in lfOel:
        sfO=sfO+fO[i]
    lfOtot.append(sfO)

### Calcul of the d30Si
# Definition of the kinetic fractionement D30Si
lfc=[]
# List the alpha (Racine(MSi28O16O16/MSi30O16O16) for all condensate (will become a dictionary)
rm=[0.9660918,0.9826074,0.9672042,0.9770084,0.9780193,0.9701425,0.9860133,0.9837388,0.9883037,0.9682458,0.9692234]
for i in range(0, len(lSigf[0])):
    alfc=0
    for j in range(0, len(lSigf)):
        if lSit[i]!=0:
            alfc=alfc+rm[j]*(lSigf[j][i]/lSit[i]) #Calcul of the alphatot of kinetic isotopic fractionation for all possible Gas condensation alpha Si
        else:
            alfc=alfc+0
    dfc=(alfc-1)*1000 #Calcul of kinetic isotopic fractionation (D30Sigas-condensates)
    if dfc==-1000:
        dfc=-16
    lfc.append(dfc)

# calculating the d30Si
ld30Si=[]
ld30Sig=[]
d30Sigaz=0
for a in range(0, len(lfSitot)):
    d30Sigaz=0+lfc[a]*np.log(lfSitot[a]) #Calcul of d30Si in the gas at each temperature
    ld30Sig.append(d30Sigaz)
    d30Sico=lfc[a]+d30Sigaz              #Calcul of d30Si in the condensates at each temperature
    ld30Si.append(d30Sico)

###Phase of generation of files for save data
num=1                                   #Number of the simulation
sub='SH'+str(num)                       #Name of the simulation + the number of the simulation

sauvfil('Pression Partielle O '+sub,lOg,lOgf,T)         #Save the oxygen partial pressure data
sauvfil('fO '+sub,lOg,lfOel,T)                          #Save the fO data
sauvfil('Composition des solides '+sub,ls,lsolid,T)     #Save the solids composition data
sauvfil('Pression Partielle Si '+sub,lSig,lSigf,T)      #Save the silicon partial pressure data
sauvfil('fSi '+sub,lSig,lfSiel,T)                       #Save the fSi data
sauvfil('Pression Partielle Al '+sub,lAlg,lAlgf,T)      #Save the aluminium partial pressure data
sauvfil('Pression Partielle Ca '+sub,lCag,lCagf,T)      #Save the calcium partial pressure data
sauvfil('Pression Partielle N '+sub,lNg,lNgf,T)         #Save the nitrogen partial pressure data

fichsauv=open('C:/Users/natha/Desktop/Modélisation d30Si '+sub+'.csv','w')
prli="T;fSi;d30Si"+"\n"
fichsauv.write(prli)
for ligne in range (0,len(T)):
    CvT=str(T[ligne])+';'
    CvfSi=str(lfSitot[ligne])+';'
    Cvd30Si=str(ld30Si[ligne])+'\n'
    Val=CvT+CvfSi+Cvd30Si
    fichsauv.write(Val)
fichsauv.close()

fichsauv=open('C:/Users/natha/Desktop/Tini de condensation '+sub+'.csv','w')
prli="Minéral;T;content"+"\n"
fichsauv.write(prli)
ligne=0
while ligne+3<=len(Tini):
    CM=str(Tini[ligne])+';'
    CT=str(Tini[ligne+1])+';'
    CC=str(Tini[ligne+2])+'\n'
    ligne=ligne+3
    Val=CM+CT+CC
    fichsauv.write(Val)
fichsauv.close()

fichsauv=open('C:/Users/natha/Desktop/T50 de condensation '+sub+'.csv','w')
prli="Minéral;T;content"+"\n"
fichsauv.write(prli)
ligne=0
while ligne+3<=len(T50):
    CM=str(T50[ligne])+';'
    CT=str(T50[ligne+1])+';'
    CC=str(T50[ligne+2])+'\n'
    ligne=ligne+3
    Val=CM+CT+CC
    fichsauv.write(Val)
fichsauv.close()

fichsauv=open('C:/Users/natha/Desktop/Tfin de condensation '+sub+'.csv','w')
prli="Minéral;T;content"+"\n"
fichsauv.write(prli)
ligne=0
while ligne+3<=len(Tfin):
    CM=str(Tfin[ligne])+';'
    CT=str(Tfin[ligne+1])+';'
    CC=str(Tfin[ligne+2])+'\n'
    ligne=ligne+3
    Val=CM+CT+CC
    fichsauv.write(Val)
fichsauv.close()


# Phase of Graphics generation
gra=2.3         #Choose the graphics

#List of color for the graphics
lc=['red','orange','yellow','green','cyan','blue','violet','pink','brown','coral','gold','lawngreen','forestgreen','turquoise','blueviolet','deeppink','chocolate','lime','fuschia','skyblue','silver','peachpuff','seagreen','lightsalmon','darkviolet','moccasin','aquamarine','firebrick','royalblue','saddlebrown','beige','rebeccapurple','navy','peru','forestgreen','dimgray','deepskyblue']


if gra==1.1:   #Graphic representing the Partial Pressure of gas with O vs T + Total Pressure of oxides in gas
    plt.plot(T,lOt, label='O Total', color='black')
    for g in range(0,len(lOg)):
        plt.plot(T,lOgf[g], label=lOg[g], color=lc[g])
    plt.xlabel('Temperature (K)')
    plt.ylabel('Partial Pressure (bar)')
    plt.yscale('log')
    plt.legend(loc='upper left', fontsize='large', ncol=2)

if gra==1.2:   #Graphic representing the Partial Pressure of gas with Si vs T + Total Pressure of gas with Si
    plt.plot(T,lSit, label='Si Total', color='black')
    for g in range(0,len(lSig)):
        plt.plot(T,lSigf[g], label=lSig[g], color=lc[g])
    plt.xlabel('Temperature (K)')
    plt.ylabel('fraction moles')
    plt.yscale('log')
    plt.legend(loc='upper left', fontsize='large', ncol=2)


if gra==2:   #Graphic representing the fSitot vs T
    plt.plot(lT,lfSit, label='fSi', color='red')
    plt.xlabel('Temperature (K)')
    plt.ylabel('fSi')
    plt.legend(loc='upper left', fontsize='large', ncol=2)

if gra==2.3: #Graphic representing the fSitot vs T, with the fSiO, fSiO2, fSi, and D30Sigas
    fig, ax1 = plt.subplots()
    plt.rcParams['font.size'] = '18'
    plt.rcParams['figure.figsize'] = 15, 10
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('fSi', color='red')
    ax1.plot(T,lfSitot,'--', label='fSi total', color='red', linewidth=3)
    ax1.plot(T,lfSiel[0],'--', label='fSi', color='orange', linewidth=3)
    ax1.plot(T,lfSiel[4],'--', label='fSiO', color='green', linewidth=3)
    ax1.plot(T,lfSiel[7],'--', label='fSiO2', color='magenta', linewidth=3)
    ax1.tick_params(axis='y', labelcolor='red')
    ax1.legend(loc='upper right',fontsize='x-small', ncol=2)
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel('Δ30Sigas', color='black')  # we already handled the x-label with ax1
    ax2.plot(T, lfc, color='black', linewidth=3)
    ax2.tick_params(axis='y', labelcolor='black')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fname='C:/Users/natha/Desktop/fSi '+sub+'.png'
    plt.savefig(fname)
    plt.close()

gra=2.4

if gra==2.4: #Graphic representing the fOtot vs T, with the fO, fO2, fOH, fSiO, fH2O, and fSiO2
    fig, ax1 = plt.subplots()
    plt.rcParams['font.size'] = '18'
    plt.rcParams['figure.figsize'] = 15, 10
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('fO', color='red')
    ax1.plot(T,lfOtot, label='fO total', color='red', linewidth=3)
    ax1.plot(T,lfOel[0],'--', label='fO', color='orange', linewidth=3)
    ax1.plot(T,lfOel[1],'--', label='fO2', color='green', linewidth=3)
    ax1.plot(T,lfOel[2],'--', label='fOH', color='blue', linewidth=3)
    ax1.plot(T,lfOel[5],'--', label='fSiO', color='magenta', linewidth=3)
    ax1.plot(T,lfOel[21],'--', label='fH2O', color='pink', linewidth=3)
    ax1.plot(T,lfOel[28],'--', label='fSiO2', color='brown', linewidth=3)
    ax1.tick_params(axis='y', labelcolor='red')
    ax1.legend(loc='best',fontsize='x-small')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fname='C:/Users/natha/Desktop/fO '+sub+'.png'
    plt.savefig(fname)
    plt.close()

if gra==2.5: #Graphic representing the fSitot vs d30Si
    plt.plot(lfSit,ld30Si, color='black')
    plt.xlabel('fSi')
    plt.ylabel('d30Si')

if gra==3:   #Graphic representing the fSitot vs d30Si
    plt.plot(lT,ld30Si, label='d30Sic', color='red')
    plt.plot(lT,ld30Sig, label='d30Sig', color='red')
    plt.xlabel('Temperature (K)')
    plt.ylabel('d30Si')
    plt.legend(loc='upper right', fontsize='large', ncol=2)

if gra==4:    #Graphic representing the molar fraction of all the condensates
    for g in range(0,len(ls)):
        plt.plot(T,lsolid[g], label=ls[g], color=lc[g])
    plt.xlabel('Temperature (K)')
    plt.ylabel('fraction moles')
    plt.yscale('log')
    plt.legend(loc='upper right', fontsize='medium', ncol=2)
