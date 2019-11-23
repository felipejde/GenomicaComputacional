library(seqinr)
library(Biostrings)
library(ggplot2)

x2  <- read.fasta("C:/Users/slope/Documents/Genomica/Practica1/Homo_sapiens_nt.fasta")
x3 <- read.fasta("C:/Users/slope/Documents/Genomica/Practica1/Saccharomyces_cerevisiae_nt.fasta")
x4 <- read.fasta("C:/Users/slope/Documents/Genomica/Practica1/Thermus_thermophilus_nt.fasta")
longitud <- lapply(homo,length)
length(x)
longitudes2 <- as.numeric(longitud)

mean(as.numeric(longitudes2))

sd(as.numeric(longitudes2))

#inciso v)
l <- data.frame(longitudesssss = longitudes2)
ggplot(l, aes(x=longitudesssss))+
  geom_bar(width=5, color="tomato3", fill="tomato3")+
  geom_vline(aes(xintercept=mean(longitudesssss)),color="red", linetype="dashed", size=1)


#2
ara <- readDNAStringSet("C:/Users/slope/Documents/Genomica/Practica1/Arabidopsis_thaliana_nt.fasta")
esche <- readDNAStringSet("C:/Users/slope/Documents/Genomica/Practica1/Escherichia_coli_nt.fasta")
homo <- readDNAStringSet("C:/Users/slope/Documents/Genomica/Practica1/Homo_sapiens_nt.fasta")
saccha <- readDNAStringSet("C:/Users/slope/Documents/Genomica/Practica1/Saccharomyces_cerevisiae_nt.fasta")
thermus <-  readDNAStringSet("C:/Users/slope/Documents/Genomica/Practica1/Thermus_thermophilus_nt.fasta")


#Arabidopsis Thaliana
GS_ara<- oligonucleotideFrequency(ara,1)[,3]
CS_ara <- oligonucleotideFrequency(ara,1)[,2]
long_ara <- width(ara)
GC_ara <- (GS_ara+CS_ara)/long_ara
GC_ara <- round(GC_ara,  2) 
nums_ara <- seq(1, 12654)
list2 <- rep("Arabidopsis",length(nums_ara))


l1 <- data.frame(ContenidoGC = GC_ara, Organismo = list2)


#Escherichia_coli
GS_e<- oligonucleotideFrequency(esche,1)[,3]
CS_e <- oligonucleotideFrequency(esche,1)[,2]
long_e <- width(esche)
GC_e <- (GS_e+CS_e)/long_e
GC_e <- round(GC_e,  2) 
nums_e <- seq(1, 4319)
list_e <- rep("Escherichia",length(nums_e))

l2 <- data.frame(ContenidoGC = GC_e, Organismo = list_e)


#Homo Sapiens
GS_h<- oligonucleotideFrequency(homo,1)[,3]
CS_h <- oligonucleotideFrequency(homo,1)[,2]
long_h <- width(homo)
GC_h <- (GS_h+CS_h)/long_h
GC_h <- round(GC_h,  2) 
nums_h <- seq(1, 1297)
list_h <- rep("Homo Sapiens",length(nums_h))

l3 <- data.frame(ContenidoGC = GC_h, Organismo = list_h)



#Saccharomyces cerevisiae
GS_s<- oligonucleotideFrequency(saccha,1)[,3]
CS_s <- oligonucleotideFrequency(saccha,1)[,2]
long_s <- width(saccha)
GC_s <- (GS_s+CS_s)/long_s
GC_s <- round(GC_s,  2) 
nums_s <- seq(1, 767)
list_s <- rep("Saccharomyces",length(nums_s))

l4 <- data.frame(ContenidoGC = GC_s, Organismo = list_s)


#Thermus_thermophilus
GS_t<- oligonucleotideFrequency(thermus,1)[,3]
CS_t <- oligonucleotideFrequency(thermus,1)[,2]

long_t <- width(thermus)
GC_t <- (GS_t+CS_t)/long_t
GC_t <- round(GC_t,  2) 
nums_t <- seq(1, 1908)
list_t <- rep("Thermus",length(nums_t))

l5 <- data.frame(ContenidoGC = GC_t, Organismo = list_t)

final_table <- rbind.data.frame(l1,l2,l3,l4,l5)
final_table_sin_ara  <- rbind.data.frame(l2,l3,l4,l5)


#2 a)
ggplot(final_table, aes(x=ContenidoGC, fill=Organismo)) +
  geom_histogram(binwidth=.01, alpha=.5, position="dodge")

ggplot(final_table, aes(ContenidoGC, fill=Organismo)) +
  geom_histogram(binwidth=.05, position="dodge")

#sin ara porque esta muy grande
ggplot(final_table_sin_ara, aes(ContenidoGC, fill=Organismo)) +
  geom_histogram(binwidth=.05, position="dodge")


#3
ara2<-Biostrings::translate(ara,if.fuzzy.codon="X")
esche2<-Biostrings::translate(esche,if.fuzzy.codon="X")
homo2<-Biostrings::translate(homo,if.fuzzy.codon="X")
saccha2<-Biostrings::translate(saccha,if.fuzzy.codon="X")
thermus2<-Biostrings::translate(thermus,if.fuzzy.codon="X")

af<-alphabetFrequency(ara2,as.prob=TRUE)
arafreq<-apply(af,2,mean)

ef<-alphabetFrequency(esche2,as.prob=TRUE)
efreq<-apply(ef,2,mean)

hf<-alphabetFrequency(homo2,as.prob=TRUE)
hfreq<-apply(hf,2,mean)

sf<-alphabetFrequency(saccha2,as.prob=TRUE)
sfreq<-apply(sf,2,mean)

tf<-alphabetFrequency(thermus2,as.prob=TRUE)
tfreq<-apply(tf,2,mean)

le = data.frame(Frecuencias = efreq, Aminoacidos = names(efreq),Organismo = rep("X",31)) 
la = data.frame(Frecuencias = arafreq, Aminoacidos = names(arafreq),Organismo = rep("X",31))
lh = data.frame(Frecuencias = hfreq, Aminoacidos = names(hfreq),Organismo = rep("X",31)) 
ls = data.frame(Frecuencias = sfreq, Aminoacidos = names(sfreq),Organismo = rep("X",31)) 
lt = data.frame(Frecuencias = tfreq, Aminoacidos = names(tfreq),Organismo = rep("X",31)) 


#Esche
e <- ggplot(le, mapping = aes(x = Aminoacidos, y = Frecuencias, fill = Organismo)) +
  geom_bar(stat = "identity") +
  ggtitle("Esche")

#ara

a <- ggplot(la, mapping = aes(x = Aminoacidos, y = Frecuencias, fill = Organismo)) +
  geom_bar(stat = "identity", color="blue") +
  ggtitle("Aridopsis Thaliana")

s <- ggplot(data = ls, mapping = aes(x = Aminoacidos, y = Frecuencias, fill = Organismo)) +
  geom_bar(stat = "identity", color="pink") +
  ggtitle("Saccha")

h <- ggplot(data = lh, mapping = aes(x = Aminoacidos, y = Frecuencias, fill = Organismo)) +
  geom_bar(stat = "identity", color="green") +
  ggtitle("homo sapiens")

te <-ggplot(data = lt, mapping = aes(x = Aminoacidos, y = Frecuencias, fill = Organismo)) +
      geom_bar(stat = "identity",color="yellow") 
     
#4

esche4 <- DNAString(esche[[1]])
matchProbePair("CTACAGCCGTTGCCGAACGT","AAAAATACTCTGCCTTTGAG",esche4)

