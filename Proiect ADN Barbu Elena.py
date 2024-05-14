import numpy as np

ls_baze_azotate_posibile=["A","C","T","G"]
def introducere_corecta(ls, n):
  print("Secventele ADN trebuie introduse folosiind urmatoarele abrevieri:\n A pentru adenina\n C pentru citozina\n T pentru timina\n G pentru guanina ")
  ls = list(input(f"Introduceti {n}. secventa de ADN: ").upper())
  introdus_corect = False
  i = 0
  introduceri_gresite = 0
  while introdus_corect == False:
      introdus_corect = True
      while i < len(ls):
          if ls[i] not in ls_baze_azotate_posibile:
              introdus_corect = False
              print(f"Ati introdus gresit pe pozitia {i + 1}\n")
              print(
                  "Secventele ADN trebuie introduse folosiind doar urmatoarele abrevieri (fara a lasa spatii libere intre ele):\n A pentru adenina\n C pentru citozina\n T pentru timina\n G pentru guanina\n ")
              introduceri_gresite += 1
              assert introduceri_gresite < 3, f"Ati introdus secventa gresit de prea multe ori."
              ls = list(input(f"Introduceti {n}. secventa de ADN: "))
              i = 0
          else:
              i = i + 1
  return ls




#### Frecventa de aparitie
def bazele_azotate_cautate():
  baze_azotate_cautate=list(input("Introduceti bazele azotate pentru care doriti sa calculam frecventa (fara spatii si semne de punctuatie intre ele): "))
  introdus_corect = False
  i = 0
  introduceri_gresite=0
  while introdus_corect == False:
      introdus_corect = True
      while i < len(baze_azotate_cautate):
          if baze_azotate_cautate[i] not in ls_baze_azotate_posibile:
              introdus_corect = False
              introduceri_gresite=introduceri_gresite+1
              if(introduceri_gresite>=3):
                 print("Ati introdus gresit bazele azotate de 3 ori")
                 baze_azotate_cautate=[]
                 return baze_azotate_cautate
              print("Secventele de ADN contin doar urmatoarele baze azotate:\n A pentru adenina\n C pentru citozina\n T pentru timina\n G pentru guanina\n ")
              baze_azotate_cautate = list(input("Introduceti dinnou "))
              i = 0
          else:
              i = i + 1
  return baze_azotate_cautate

def frecv_baze_azotate(ls_secventa,baze_azotate_cautate):
  dct=dict()
  for baza_azotata in baze_azotate_cautate:
    if baza_azotata not in dct.keys():
       dct[baza_azotata]=ls_secventa.count(baza_azotata)
  return(dct)



### Transcriptia ADN-ului in ARN
def transcriptia(adn_ls):
    arn_ls=[]
    for baza_azotata in adn_ls:
        if baza_azotata=="A":
            arn_ls.append("U")
        elif baza_azotata=="C":
            arn_ls.append("G")
        elif baza_azotata=="T":
            arn_ls.append("A")
        elif baza_azotata=="G":
            arn_ls.append("C")
    return arn_ls



# presupunem ca splicing-ul inlatura daca este necesar bazele azotate de pe ultimele pozitii din ARN mesager, astfel incat sa se obtina ARN-ul mesager matur care contine un numar de baze azotate care este un multiplu a lui 3
def splicing(arnM_ls):
    while (len(arnM_ls) % 3 !=0):
        arnM_ls.pop()
    return arnM_ls


# Aminoacizii codificati de codonii din ARN mesager matur
def translatia(arn_ls):
  aminoacizi_ls=[]
  for i in range(0,len(arn_ls),3):
    #ls=[]
    ls=[arn_ls[i],arn_ls[i+1],arn_ls[i+2]]
    aminoacid="".join(map(str,ls))
    if (aminoacid=="UUU" or aminoacid=="UUC") :
      aminoacizi_ls.append("Phe")
    elif (aminoacid=="UUA" or aminoacid=="UUG" or aminoacid=="CUU" or aminoacid=="CUC" or aminoacid=="CUA" or aminoacid=="CUG"):
      aminoacizi_ls.append("Leu")
    elif (aminoacid=="AUU" or aminoacid=="AUC" or aminoacid=="AUA"):
      aminoacizi_ls.append("Ile")
    elif aminoacid=="AUG":
      aminoacizi_ls.append("Met")
    elif (aminoacid=="GUU" or aminoacid=="GUC" or aminoacid=="GUA" or aminoacid=="GUG"):
      aminoacizi_ls.append("Val")
    elif (aminoacid=="UCU" or aminoacid=="UCC" or aminoacid=="UCA" or aminoacid=="UCG" or aminoacid=="AGU" or aminoacid=="AGC"):
      aminoacizi_ls.append("Ser")
    elif (aminoacid=="CCU" or aminoacid=="CCC" or aminoacid=="CCA" or aminoacid=="CCG"):
      aminoacizi_ls.append("Pro")
    elif (aminoacid=="ACU" or aminoacid=="ACC" or aminoacid=="ACA" or aminoacid=="ACG"):
      aminoacizi_ls.append("Thr")
    elif (aminoacid=="GCU" or aminoacid=="GCC" or aminoacid=="GCA" or aminoacid=="GCG"):
      aminoacizi_ls.append("Ala")
    elif (aminoacid=="UAU" or aminoacid=="UAC"):
      aminoacizi_ls.append("Tyr")
    elif (aminoacid=="UAA" or aminoacid=="UAG" or aminoacid=="UGA"):
      aminoacizi_ls.append("STOP")
    elif (aminoacid=="CAU" or aminoacid=="CAC"):
      aminoacizi_ls.append("His")
    elif (aminoacid=="CAA" or aminoacid=="CAG"):
      aminoacizi_ls.append("Gln")
    elif (aminoacid=="AAU" or aminoacid=="AAC"):
      aminoacizi_ls.append("Asn")
    elif (aminoacid=="AAA" or aminoacid=="AAG"):
      aminoacizi_ls.append("Lys")
    elif (aminoacid=="GAU" or aminoacid=="GAC"):
      aminoacizi_ls.append("Asp")
    elif (aminoacid=="GAA" or aminoacid=="GAG"):
      aminoacizi_ls.append("Glu")
    elif (aminoacid=="UGU" or aminoacid=="UGC"):
      aminoacizi_ls.append("Cys")
    elif aminoacid=="UGG":
      aminoacizi_ls.append("Trp")
    elif (aminoacid=="CGU" or aminoacid=="CGC" or aminoacid=="CGA" or aminoacid=="CGG" or aminoacid=="AGA" or aminoacid=="AGG"):
      aminoacizi_ls.append("Arg")
    elif (aminoacid=="GGU" or aminoacid=="GGC" or aminoacid=="GGA" or aminoacid=="GGG"):
      aminoacizi_ls.append("Gly")
  return aminoacizi_ls



# Frecventa de aparitie a tuturor aminoacizilor
aminoacizi_posibili_ls=["Met","Phe","Leu","Ile","Val","Ser","Pro","Thr","Ala","Tyr","His","Gln","Asn","Lys","Asp","Glu","Cys","Trp","Arg","Gly","STOP"]
def frecv_aminoacizi(proteina_ls,aminoacizi_posibili_ls):
    frecvente_aminoacizi_dct=dict()
    for aminoacid in aminoacizi_posibili_ls:
      if aminoacid not in frecvente_aminoacizi_dct.keys():
        frecvente_aminoacizi_dct[aminoacid]=proteina_ls.count(aminoacid)
    return frecvente_aminoacizi_dct


# # calculul distantei Hamming:
def distanta_Hamming(adn1_ls,adn2_ls):
  lst_pozit_diferite=[]
  if (len(adn1_ls)==len(adn2_ls)):
    distanta_hamming=0
    for i in range(0,len(adn1_ls)):
        if adn1_ls[i] != adn2_ls[i]:
            lst_pozit_diferite.append(i+1) # vrem sa afisam pozitia (nu indexul) in care se afla baze azotate diferite
            distanta_hamming+=1
    print(adn1_ls)      
    print(adn2_ls)  
    print("Distanta Hamming: ",distanta_hamming)     
    print("Lista pozitiilor in care difera bazele azotate ale celor doua secvente de ADN ", lst_pozit_diferite)
  else:
    print("Distanta Hamming nu poate fi calculata pentru ca cele doua secvente de ADN nu au lungimi identice.\n")
        


# Calculul celui mai lung subsir comun(cmlsc) folosiind programarea dinamica
def dezv_rel_rec_cmlsc(a,b):
  m=len(a) #numarul de linii
  n=len(b) #numarul de coloane
  D=[]
  # crearea matricei goale
  for i in range(0,m):
    ls=[]
    for j in range(0,n):
      ls.append(0)
    D.append(ls)

  # dezvoltarea relatiei de recurenta
  for i in range(1,m):
    for j in range(1,n):
      if a[i] == b[j]:
        D[i][j]= int(D[i-1][j-1])+1
      elif a[i] != b[j]:
        D[i][j]= max(int((D[i-1][j])), int((D[i][j-1])))
  return D



def dimens_cmlsc(a,b,D):
  # construirea solutiei
  i= len(a)-1
  j= len(b)-1
  solutia_cmlsc=[]
  print("Dimensiunea celui mai lung subsir comun: ", D[i][j])
  while D[i][j]>0:
    if a[i]==b[j]:
      solutia_cmlsc.append(a[i])
      i=i-1
      j=j-1
    elif D[i][j] == D[i-1][j]:
        i=i-1
    else:
        j=j-1
  solutia_cmlsc.reverse()
  print(solutia_cmlsc)




# calculul DISTANTEI DE EDITARE
def dezv_rel_rec_de(a,b):
  m=len(a)
  n=len(b)
  S=[]
  # crearea matricei goale
  for i in range(0,m):
    ls=[]
    for j in range(0,n):
      ls.append(0)
    S.append(ls)
  # completarea primei linii
  for j in range(0,n):
    S[0][j]= 1* j
  # completarea primei coloane
  for i in range(0,m):
    S[i][0]= 1*i
  for i in range(1,m):
    for j in range(1,n):
      if a[i] == b[j]:
        S[i][j]= int(S[i-1][j-1])
      elif a[i] != b[j]:
        S[i][j]= min(int(S[i-1][j-1]), int(S[i][j-1]),int(S[i-1][j]))+1
  return S




def calc_dist_editare(a,b,S):
  i= len(a)-1
  j= len(b)-1
  k= 0
  secventa=[]
  print("Distanta de editare: ", S[i][j])
  substitutii=0
  insertii=0
  stergeri=0
  while S[i][j]>0:
    if a[i]==b[j]:
      secventa.append(b[j])
      i=i-1
      j=j-1
    elif S[i][j]==(int(S[i-1][j-1]))+1:
      secventa.append(b[j])
      i=i-1
      j=j-1
      substitutii+=1
    elif S[i][j]==(int(S[i][j-1]))+1:
      secventa.append("-")
      j=j-1
      insertii+=1
    else:
        i=i-1
        stergeri+=1
  secventa.reverse()
  print(secventa)
  print(f"Au fost efectuate {substitutii} substitutii, {insertii} insertii si {stergeri} stergeri.")



################################### MENIU CU OPTIUNI #########################################

# Introducerea corecta a secventelor de ADN
adn1_ls = list()
adn1_ls = introducere_corecta(adn1_ls, 1)
print("Prima secventa de ADN:\n",adn1_ls)
adn2_ls = list()
adn2_ls = introducere_corecta(adn2_ls, 2)
print("A doua secventa de ADN\n",adn2_ls)




while(True):
  print("----------------------------------------------------------------------------------------------------------------------------------------")
  print("Alegeti una dintre optiunile disponibile (introduceti doar cifra din fata optiunii): \n")
  optiune=(input(" 1. Determinarea frecventelor de aparitie ale bazelor azotate\n 2. Transcriptia secventelor de ADN in secvente de ARN\n 3. Translatia in secvente de aminoacizi\n 4. Determinarea frecventei de aparitie a aminoacizilor\n 5. Calculul distantei Hamming pentru secventele de ADN\n 6. Determinarea celui mai lung subsir comun al secventelor de ADN\n 7. Calculul distantei de editare a secventelor de ADN\n 8. Iesire din program\n"))


  if(optiune=="1" or optiune=="1." or optiune=="1. Determinarea frecventelor de aparitie ale bazelor azotate"):
    # Calcul baze azotate cautate
    baze_azotate_cautate=bazele_azotate_cautate()
    adn1_dct=frecv_baze_azotate(adn1_ls,baze_azotate_cautate)
    print("In prima secventa de ADN frecventa de aparitie a bazelor azotate cautate este: ",adn1_dct)
    adn2_dct=frecv_baze_azotate(adn2_ls,baze_azotate_cautate)
    print("In a doua secventa de ADN frecventa de aparitie a bazelor azotate cautate este: ",adn2_dct)
    print("----------------------------------------------------------------------------------------------------------------------------------------")
  elif(optiune=="2"or optiune=="2." or optiune=="2. Transcriptia secventelor de ADN in secvente de ARN"):
    # obtinerea ARN-ului mesager
    arn1m_ls=transcriptia(adn1_ls)
    print("ARN mesager 1: ",arn1m_ls)
    arn2m_ls=transcriptia(adn2_ls)
    print("ARN mesager 2: ",arn2m_ls)
    print("----------------------------------------------------------------------------------------------------------------------------------------")
  elif(optiune=="3" or optiune=="3." or optiune=="3. Translatia in secvente de aminoacizi"):
    # obtinerea ARN-ului mesager
    arn1m_ls=transcriptia(adn1_ls)
    arn2m_ls=transcriptia(adn2_ls)

    # obtinerea ARN-uluui matur
    print("Presupunem ca splicing-ul inlatura daca este necesar bazele azotate de pe ultimele pozitii din ARN mesager, astfel incat sa se obtina ARN mesager matur care contine un numar de baze azotate care este un multiplu a lui 3")
    arn1M_ls=splicing(arn1m_ls)
    print("ARN mesager matur 1 : ",arn1M_ls)
    arn2M_ls=splicing(arn2m_ls)
    print("ARN mesager matur 2 : ",arn2M_ls)

    # Translatia codonilor in aminoacizi:
    proteina1_ls=translatia(arn1M_ls)
    print("Proteina 1: ",proteina1_ls)
    proteina2_ls=translatia(arn2M_ls)
    print("Proteina 2: ",proteina2_ls)
    print("----------------------------------------------------------------------------------------------------------------------------------------")
  elif(optiune=="4" or optiune=="4."or optiune=="4. Determinarea frecventei de aparitie a aminoacizilor"):
    arn1m_ls=transcriptia(adn1_ls)
    arn2m_ls=transcriptia(adn2_ls)

    # obtinerea ARN mesager matur
    arn1M_ls=splicing(arn1m_ls)
    arn2M_ls=splicing(arn2m_ls)

    # Translatia codonilor in aminoacizi:
    proteina1_ls=translatia(arn1M_ls)
    proteina2_ls=translatia(arn2M_ls)

    frecv_aa_proteina1_dct=frecv_aminoacizi(proteina1_ls,aminoacizi_posibili_ls)
    print("Frecventa de aparitie a aminoacizilor in proteina 1 ",frecv_aa_proteina1_dct)

    frecv_aa_proteina2_dct=frecv_aminoacizi(proteina2_ls,aminoacizi_posibili_ls)
    print("Frecventa de aparitie a aminoacizilor in proteina 2 ",frecv_aa_proteina2_dct)

    # Introducerea frecventelor de aparitie a aminoacizilor in fisier csv

    import csv
    with open("frecvente_aa_priteine.csv","w",newline="") as t1:
        writer=csv.writer(t1)
        header=["Aminoacidul","Frecventa de aparitie in proteina 1","Frecventa de aparitie in proteina 2"]
        writer.writerow(header)
        for aminoacid in frecv_aa_proteina1_dct:
            linie=[aminoacid,str(frecv_aa_proteina1_dct[aminoacid]),str(frecv_aa_proteina2_dct[aminoacid])]
            writer.writerow(linie)


    # crearea graficului care afiseaza frecventa de aparitie a aminoacizilor in cele doua proteine
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib import style

    frecv_aa_proteina1_ls=list(frecv_aa_proteina1_dct.values())
    frecv_aa_proteina2_ls=list(frecv_aa_proteina2_dct.values())

    values=np.arange(len(aminoacizi_posibili_ls)) # pentru a putea calcula pozitia barelor alaturate

    width=0.3
    plt.bar(values,frecv_aa_proteina1_ls,width,label="Proteina 1")
    plt.bar(values+width,frecv_aa_proteina2_ls,width,label="Proteina 2")

    plt.xticks(values,aminoacizi_posibili_ls) # pentru a inlocui acele values cu denumirile aminoacizilor

    plt.legend(loc="upper right")
    plt.xlabel("Aminoacizi")
    plt.ylabel("Frecventa de aparitie")
    plt.title("Frecventa de aparitie a aminoacizilor in cele doua proteine")
    plt.show()
    print("----------------------------------------------------------------------------------------------------------------------------------------")
  elif(optiune=="5" or optiune=="5." or optiune=="5. Calculul distantei Hamming pentru secventele de ADN"):
    distanta_Hamming(adn1_ls,adn2_ls)
    print("----------------------------------------------------------------------------------------------------------------------------------------")
  elif(optiune=="6" or optiune=="6." or optiune=="6. Determinarea celui mai lung subsir comun al secventelor de ADN"):
    # pregatirea sirurilor pentru dezvoltarea relatiei de recurenta si aplicarea ei
    a=["0"]
    b=["0"]
    for element in adn1_ls:
        a.append(element)
    for element in adn2_ls:
        b.append(element)
    # aplicarea relatiei de recurenta:    
    D=dezv_rel_rec_cmlsc(a,b)
    for linie in D:
       print(linie)
    # gasirea solutiei
    dimens_cmlsc(a,b,D)
    print("----------------------------------------------------------------------------------------------------------------------------------------")
  elif(optiune=="7" or optiune=="7." or optiune=="7. Calculul distantei de editare a secventelor de ADN"):
    # pregatirea sirurilor pentru dezvoltarea relatiei de recurenta si aplicarea ei
    a1=["0"]
    b1=["0"]
    for element in adn1_ls:
        a1.append(element)
    for element in adn2_ls:
        b1.append(element)
    S=dezv_rel_rec_de(a1,b1)
    for linie in S:
       print(linie)
    # gasirea solutiei
    calc_dist_editare(a1,b1,S)
    print("----------------------------------------------------------------------------------------------------------------------------------------")
  elif(optiune=="8" or optiune=="8." or optiune=="8. Iesire din program"):
     exit()






















      