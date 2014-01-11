#What is the genetic network of three co-occurring species within the same plant family from the Southwest?
#Data from Sharon and Erika
#1. Analyze each species separately then compare	
#2. Locus specific features can be looked at with FSAT and Genepop	
#
#Na
#
#Ho
#
#He
#
#Fis, Fst etc…departures from HWE and LD
#3. Null alleles and correction with Microchecker and FreeNA	
#4. Intersite FSt values, AMOVA etc..Arlequin or GenAlex	
#5. Genepop intersite Fst on Geographic distance -	

#Genetic Networks from 
source('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/Genetic_Networks/ind_net.R')
library(igraph)
source('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/Genetic_Networks/igraph_adds.R')
library(RColorBrewer)


#Fremont network
pf.ms=read.csv('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/PVERgeneFLOW/PVER_PF.csv')

#Separate genetic data from env data 
x=pf.ms[,-1:-2]
y=pf.ms[,1:2]
pop=y[,2]

#DEAL WITH MISSING ALLELES
x[is.na(x)==TRUE]=-1
x=allele.patch(x,method='random')

t1=table(pf.ms[,2])
pf.ms=na.omit(pf.ms)
t2=table(pf.ms[,2])
barplot(as.numeric(((t1-t2)/t1))[order(as.numeric(((t1-t2)/t1)),decreasing=TRUE)],names=names(t1)[order(as.numeric(((t1-t2)/t1)),decreasing=TRUE)],las=2,ylab='Percent Missing Alleles')

image(x==-1)


#Take the first three letters of the pop names
pop=as.character(pop)
for (i in 1:length(pop)){
   z=strsplit(pop[i],split='')[[1]]
   pop[i]=paste(z[1],z[2],z[3],z[4],sep='')
	}

#separate the alleles for all the loci into two matrixes
x1x2=allele.split(x)
x1=x1x2[[1]]
x2=x1x2[[2]]

#create a vector for the allele names
alleles=unique(vectorize(cbind(x1,x2)))

#create the codification matrix
C=codify(x1,x2)
colnames(C)=alleles

#Obtain the centroids for populations by averaging the multivariate coding vectors for each population
C.=centroids(C,pop)
pk=table(vectorize(cbind(x1,x2))) #alleleic frequencies
K=length(alleles)

#Genetic distance
#dij^2 = (1/2)*sum((1/K*pk)*(yik-yjk)^2)
dij=gdist(C.,pk,K)

#Calculate the independence network
net.=ind.net(as.matrix(dij),nrow(x),alpha=0.05)

#Convert partial correlations to relative strengths
sij=(-1/2)*log((1-net.^2))
heatmap(sij)

pf.sij=sij

sij=pf.sij

#Fremont Graph
graph.=as.matrix(round(sij,5))
graph.sij=graph.adjacency(graph.,mode='undirected',weighted=TRUE)
#sij.wtc=walktrap.community(graph.sij,steps=4)
sij.sgc=spinglass.community(graph.sij,spins=30,gamma=1)
###sij.mod=sij.wtc$membership+1
#sij.mod=sij.sgc$membership+1

net.mod=0
net.mod[1]=sij.sgc$modularity

#Using the spiinglass modularity output
sij.mod=scan()
1 1 3 4 1 3 1 4 1 4 2 2 2 4 5

pf.mod=sij.mod

#par(bg='black')
g=graph.sij
vertex.label=rownames(graph.)
vertex.label=c('AF','AL','CNWR','FC','GRNCA','GR3W','HP','ML','PVER','RL','SP','SCR','SC','HD','VRBD')
vertex.color=brewer.pal(length(unique(sij.mod)),'Set1')[sij.mod]
my.layout=layout.fruchterman.reingold(g)
vertex.label.color='white'
vertex.frame.color='white'
vertex.label.family='Helvetica'
vertex.size=30
edge.label=graph.[lower.tri(graph.)]
edge.label=edge.label[edge.label>0]
edge.label.color='white'
edge.label.family='Helvetica'
edge.label.cex=0.8
edge.width=(edge.label*10)^1.5
edge.color='lightgrey'

rownames(pf.sij)=colnames(pf.sij)=vertex.label

plot.igraph(x=g,vertex.label=vertex.label,layout=my.layout,vertex.size=vertex.size,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label=round(edge.label,2),edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

#Remove links less than 0.01
graph.=as.matrix(round(sij,5))
graph.[graph.<0.01]=0
graph.sij=graph.adjacency(graph.,mode='undirected',weighted=TRUE)
g=graph.sij
edge.label=graph.[lower.tri(graph.)]
edge.label=edge.label[edge.label>0]
edge.label.cex=0.8
edge.width=(edge.label*10)^1.5

plot.igraph(x=g,vertex.label=vertex.label,vertex.size=vertex.size,layout=my.layout,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)


edge.label=''
vertex.label=''
plot.igraph(x=g,vertex.label=vertex.label,vertex.size=vertex.size,layout=my.layout,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)


#hand positioning
graph.=round(sij,2)
graph.[graph.==0]<-0.01
diag(graph.)=0
graph.sij=graph.adjacency(graph.,mode='undirected',weighted=TRUE)
g=graph.sij
edge.label=graph.[lower.tri(graph.)]
edge.label=edge.label[edge.label>0]
edge.width=5
tkplot(g,layout=my.layout,vertex.label=vertex.label,vertex.color='gray',vertex.label.color='white',edge.width=edge.width,vertex.label.family=vertex.label.family,vertex.size=30,edge.label='',edge.label.family=edge.label.family,edge.label.color=edge.label.color,vertex.frame.color=vertex.frame.color)

#Gooding's Willow network
sg.ms=read.csv('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/PVERgeneFLOW/PVER_SG.csv')

#DEAL WITH MISSING ALLELES
#t1=table(sg.ms[,2])
#sg.ms=na.omit(sg.ms)
#t2=table(sg.ms[,2])
#plot(((t1-t2)/t1))

sg.ms[is.na(sg.ms)]=-1

#Separate genetic data from env data 
x=sg.ms[,-1:-2]
y=sg.ms[,1:2]
pop=y[,2]

par(mfrow=c(1,2))
image(x==-1)
x=allele.patch(x)
image(x==-1)

#Take the first three letters of the pop names
pop=as.character(pop)
for (i in 1:length(pop)){
   z=strsplit(pop[i],split='')[[1]]
   pop[i]=paste(z[1],z[2],z[3],z[4],sep='')
	}

#separate the alleles for all the loci into two matrixes
x1x2=allele.split(x)
x1=x1x2[[1]]
x2=x1x2[[2]]

#create a vector for the allele names
alleles=unique(vectorize(cbind(x1,x2)))

#create the codification matrix
C=codify(x1,x2)
colnames(C)=alleles

#Obtain the centroids for populations by averaging the multivariate coding vectors for each population
C.=centroids(C,pop)
pk=table(vectorize(cbind(x1,x2))) #alleleic frequencies
K=length(alleles)

#Genetic distance
#dij^2 = (1/2)*sum((1/K*pk)*(yik-yjk)^2)
dij=gdist(C.,pk,K)

#Calculate the independence network
net.=ind.net(as.matrix(dij),nrow(x),alpha=0.05)

#Convert partial correlations to relative strengths
sij=(-1/2)*log((1-net.^2))
heatmap(sij)

sg.sij=sij

sij=sg.sij

#Goodingii Graph
graph.=as.matrix(round(sij,5))
graph.sij=graph.adjacency(graph.,mode='undirected',weighted=TRUE)
sij.sgc=spinglass.community(graph.sij,spins=30,gamma=1)
#sij.mod=sij.sgc$membership+1

net.mod[2]=sij.sgc$modularity

sij.mod=scan()
1 3 1 2 1 1 1 3 3 2 1 3 2 2

sg.mod=sij.mod

#par(bg='black')
g=graph.sij
vertex.label=rownames(graph.)
vertex.label=c('AF','AL','CNWR','FC','GRNCA','GR3W','HP','ML','PVER','RL','SP','SCR','SC','HD')
vertex.color=brewer.pal(length(unique(sij.mod)),'Set1')[sij.mod]
my.layout=layout.fruchterman.reingold(g)
vertex.label.color='white'
vertex.frame.color='white'
vertex.label.family='Helvetica'
vertex.size=30
edge.label=graph.[lower.tri(graph.)]
edge.label=edge.label[edge.label>0]
edge.label.color='white'
edge.label.family='Helvetica'
edge.label.cex=0.8
edge.width=(edge.label*10)^1.5
edge.color='lightgrey'

rownames(sg.sij)=colnames(sg.sij)=vertex.label

plot.igraph(x=g,vertex.label=vertex.label,layout=my.layout,vertex.size=vertex.size,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label=round(edge.label,2),edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

#Remove edge labels
edge.label=graph.[lower.tri(graph.)]
edge.label=edge.label[edge.label>0]
edge.label.cex=0.8
edge.width=(edge.label*10)^1.5
edge.label=''

plot.igraph(x=g,vertex.label=vertex.label,vertex.size=vertex.size,layout=my.layout,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

#Remove edges that are less than 0.01
graph.=as.matrix(round(sij,5))
graph.[graph.<0.01]<-0
graph.sij=graph.adjacency(as.matrix(graph.),mode='undirected',weighted=TRUE)
g=graph.sij
edge.label=graph.[lower.tri(graph.)]
edge.label=edge.label[edge.label>0]
edge.label.cex=0.8
edge.width=(edge.label*10)^1.5
edge.label=''

plot.igraph(x=g,vertex.label=vertex.label,vertex.size=vertex.size,layout=my.layout,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

vertex.label=''
plot.igraph(x=g,vertex.label=vertex.label,vertex.size=vertex.size,layout=my.layout,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

#hand positioning

tkplot(g,layout=my.layout,vertex.label=vertex.label,vertex.color=vertex.color,vertex.label.color='white',edge.width=edge.width)

#Coyote Willow network
se.ms=read.csv('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/PVERgeneFLOW/PVER_SE.csv')

#t1=table(se.ms[,2])
#se.ms=na.omit(se.ms)
#t2=table(se.ms[,2])
#plot(((t1-t2)/t1))

se.ms[is.na(se.ms)==TRUE]=-1

#Separate genetic data from env data 
x=se.ms[,-1:-2]
y=se.ms[,1:2]
pop=y[,2]

#DEAL WITH MISSING ALLELES

x=allele.patch(x,method='random')
image(is.na(x))

#Take the first three letters of the pop names
pop=as.character(pop)
z=character()
for (i in 1:length(pop)){
   z=strsplit(pop[i],split='')[[1]]
   pop[i]=paste(z[1],z[2],z[3],z[4],sep='')
	}
pop=factor(pop)

#separate the alleles for all the loci into two matrixes
x1x2=allele.split(x)
x1=x1x2[[1]]
x2=x1x2[[2]]

#create a vector for the allele names
alleles=unique(vectorize(cbind(x1,x2)))

#create the codification matrix
C=codify(x1,x2)
colnames(C)=alleles

#Obtain the centroids for populations by averaging the multivariate coding vectors for each population
C.=centroids(C,pop)
pk=table(vectorize(cbind(x1,x2))) #alleleic frequencies
K=length(alleles)

#Genetic distance
#dij^2 = (1/2)*sum((1/K*pk)*(yik-yjk)^2)
dij=gdist(C.,pk,K)

#Calculate the independence network
net.=ind.net(as.matrix(dij),nrow(x),alpha=0.05)

#Convert partial correlations to relative strengths
sij=(-1/2)*log((1-net.^2))
heatmap(sij)

se.sij=sij
sij=se.sij

#Coyote Graph
graph.=as.matrix(round(sij,5))
graph.sij=graph.adjacency(graph.,mode='undirected',weighted=TRUE)

#sij.wtc=walktrap.community(graph.sij,steps=4)
#sij.mod=sij.wtc$membership+1
sij.sgc=spinglass.community(graph.sij,spins=30,gamma=1)
#sij.mod=sij.sgc$membership+1

net.mod[3]=sij.sgc$modularity

sij.mod=scan()
1 2 3 3 3 2 3 2 2

se.mod=sij.mod

#par(bg='black')
g=graph.sij
vertex.label=rownames(graph.)
vertex.label=c('AL','CNWR','DD','GRNCA','GR3W','ML','PVER','HD','VRBD')
vertex.color=brewer.pal(length(unique(sij.mod)),'Set1')[sij.mod]
my.layout=layout.fruchterman.reingold(g)
vertex.label.color='white'
vertex.frame.color='white'
vertex.label.family='Helvetica'
vertex.size=30
edge.label=graph.[lower.tri(graph.)]
edge.label=edge.label[edge.label>0]
edge.label.color='white'
edge.label.family='Helvetica'
edge.label.cex=0.8
edge.width=(edge.label*10)^1.5
edge.color='lightgrey'

rownames(se.sij)=colnames(se.sij)=vertex.label

plot.igraph(x=g,vertex.label=vertex.label,layout=my.layout,vertex.size=vertex.size,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label=round(edge.label,2),edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

#Remove edge labels
edge.label=''
plot.igraph(x=g,vertex.label=vertex.label,layout=my.layout,vertex.size=vertex.size,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label=edge.label,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

#Remove links less than 0.01
graph.=as.numeric(sij)
graph.[graph.<0.01]<-0
graph.sij=graph.adjacency(graph.,mode='undirected',weighted=TRUE)
g=graph.sij
edge.label=graph.[lower.tri(graph.)]
edge.label=edge.label[edge.label>0]
edge.label.cex=0.8
edge.width=(edge.label*10)^1.5
edge.label=''

plot.igraph(x=g,vertex.label=vertex.label,layout=my.layout,vertex.size=vertex.size,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label=edge.label,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

vertex.label=''
plot.igraph(x=g,vertex.label=vertex.label,layout=my.layout,vertex.size=vertex.size,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label=edge.label,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

#hand positioning

tkplot(g,layout=my.layout,vertex.label=vertex.label,vertex.color=vertex.color,vertex.label.color='white',edge.width=edge.width)


#Format for STRUCTURE
pf=read.csv('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/PVERgeneFLOW/PVER_PF.csv')
se=read.csv('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/PVERgeneFLOW/PVER_SE.csv')
sg=read.csv('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/PVERgeneFLOW/PVER_SG.csv')

pf[is.na(pf)]=-9
se[is.na(se)]=-9
sg[is.na(sg)]=-9

nrow(pf)
nrow(se)
nrow(sg)

(ncol(pf)-2)
(ncol(se)-2)
(ncol(sg)-2)

#pf[,1]=as.numeric(pf[,1])
#pf[,2]=as.numeric(pf[,2])

pf.q=data.frame(pf[,2],as.numeric(pf[,2]))
pf.q[order(pf.q[,2]),]

write.table(pf,file='/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/PVER_STRUCTURE/PVER_pfSTR.txt',sep=' ',row.names=FALSE,col.names=FALSE)


write.table(se,file='/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/PVER_STRUCTURE/PVER_seSTR.txt',sep=' ',row.names=FALSE)
write.table(sg,file='/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/PVER_STRUCTURE/PVER_sgSTR.txt',sep=' ',row.names=FALSE)

#Pf posterior probabilities
pf.lp=scan()
-4656.4 -4456.9 -4313.2 -4275.6 -4152.2 -4121.2  -4143.6 -4332.5 -4160.7

pf.lp=pf.lp-pf.lp[pf.lp==max(pf.lp)]

pf.pp=pf.lp*0

for (i in seq(along=pf.pp)){
	pf.pp[i]=exp(pf.lp[i])/sum(exp(pf.lp))
	}

pf.pp

plot(pf.pp~seq(along=pf.pp),xlab='K',ylab='Posterior Probability',type='l')

pf.clust=read.csv('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/PVER_STRUCTURE/PVER_pfSTRclust.csv')

pf.clust[,1]=levels(pf[,2])

for (i in seq(along=pf.clust[,1])){
	pf.clust[i,8]=(1:6)[pf.clust[i,2:7]==apply(pf.clust[,2:7],1,max)[i]]
	}

#Node statistics
setwd('/Users/Aeolus/Documents/Active_Projects/Cwood2011/Cwood2011Presentation/images')

pdf('pf_degree.pdf',pointsize=15)
par(col.lab='white',col.axis='white',col='white',mai=c(1.3, 1.5, 0.82, 0.42))
pf.deg=degree(pf.sij,rescale=TRUE)
names(pf.deg)=rownames(pf.sij)
barplot(sort(pf.deg,decreasing=TRUE),las=2,col=brewer.pal(length(unique(pf.mod)),'Set1')[pf.mod[order(pf.deg,decreasing=TRUE)]],ylab='Centralization')
axis(2,col='white',labels=FALSE)
dev.off()

pdf('se_degree.pdf',pointsize=15)
par(col.lab='white',col.axis='white',col='white',mai=c(1.3, 1.5, 0.82, 0.42))
se.deg=degree(se.sij,rescale=TRUE)
names(se.deg)=rownames(se.sij)
barplot(sort(se.deg,decreasing=TRUE),las=2,col=brewer.pal(length(unique(se.mod)),'Set1')[se.mod[order(se.deg,decreasing=TRUE)]],ylab='Centralization')
axis(2,col='white',labels=FALSE)
dev.off()

pdf('sg_degree.pdf',pointsize=15)
par(col.lab='white',col.axis='white',col='white',mai=c(1.3, 1.5, 0.82, 0.42))
sg.deg=degree(sg.sij,rescale=TRUE)
names(sg.deg)=rownames(sg.sij)
barplot(sort(sg.deg,decreasing=TRUE),las=2,col=brewer.pal(length(unique(sg.mod)),'Set1')[sg.mod[order(sg.deg,decreasing=TRUE)]],ylab='Centralization')
axis(2,col='white',labels=FALSE)
dev.off()

#Network Centrality

pf.C=centralization(pf.sij,degree,normalize=TRUE)
pf.rgC=centralization(rgraph(nrow(pf.sij),5000),degree,normalize=TRUE)

sg.C=centralization(sg.sij,degree,normalize=TRUE)
sg.rgC=centralization(rgraph(nrow(sg.sij),5000),degree,normalize=TRUE)

se.C=centralization(se.sij,degree,normalize=TRUE)
se.rgC=centralization(rgraph(nrow(se.sij),5000),degree,normalize=TRUE)

net.ses=-1*c((pf.C-mean(pf.rgC))/sd(pf.rgC),(sg.C-mean(sg.rgC))/sd(sg.rgC),(se.C-mean(se.rgC))/sd(se.rgC))

names(net.ses)=c('P.f.','S.g.','S.e.')

pdf('net_centralization.pdf',pointsize=15)
par(col.lab='white',col.axis='white',col='white',mai=c(1.3, 1.5, 0.82, 0.42))
barplot(net.ses,col=brewer.pal(length(net.ses),'Set1')[1:length(net.ses)],ylab='Standardized Centrality',ylim=c(0,3))
axis(2,col='white',labels=FALSE)
dev.off()

#Network modularity
names(net.mod)=c('P.f.','S.g.','S.e.')

pdf('net_modularity.pdf',pointsize=15)
par(col.lab='white',col.axis='white',col='white',mai=c(1.3, 1.5, 0.82, 0.42))
barplot(net.mod,col=brewer.pal(length(net.ses),'Set1')[1:length(net.ses)],ylab='Modularity',ylim=c(0,0.35))
axis(2,col='white',labels=FALSE)
dev.off()

#Edge distributions
pf.hist=hist(pf.sij[lower.tri(pf.sij)],main='',xlab='Edge Weight')
sg.hist=hist(sg.sij[lower.tri(sg.sij)],main='',xlab='Edge Weight')
se.hist=hist(se.sij[lower.tri(se.sij)],main='',xlab='Edge Weight')

pdf('pf_hist.pdf',pointsize=15)
par(col.lab='white',col.axis='white',col='white',mai=c(1.3, 1.5, 0.82, 0.42))
plot(pf.hist,xlab='Edge Weight',main='',axes=FALSE,col=brewer.pal(length(net.ses),'Set1')[1])
axis(1,col='white',labels=TRUE)
axis(2,col='white',labels=TRUE)
dev.off()

pdf('sg_hist.pdf',pointsize=15)
par(col.lab='white',col.axis='white',col='white',mai=c(1.3, 1.5, 0.82, 0.42))
plot(sg.hist,xlab='Edge Weight',main='',axes=FALSE,col=brewer.pal(length(net.ses),'Set1')[2])
axis(1,col='white',labels=TRUE)
axis(2,col='white',labels=TRUE)
dev.off()

pdf('se_hist.pdf',pointsize=15)
par(col.lab='white',col.axis='white',col='white',mai=c(1.3, 1.5, 0.82, 0.42))
plot(se.hist,xlab='Edge Weight',main='',axes=FALSE,col=brewer.pal(length(net.ses),'Set1')[3])
axis(1,col='white',labels=TRUE)
axis(2,col='white',labels=TRUE)
dev.off()

pdf('net_edgeD.pdf',pointsize=15)
par(col.lab='white',col.axis='white',col='white',mai=c(1.3, 1.5, 0.82, 0.42))
plot(log(pf.hist$counts+1)~log(pf.hist$mids),pch=19,xlab='Log Edge Weight',ylab='Log Frequency',axes=FALSE,cex=1.5,xlim=c(-4,0.25),ylim=c(-0.5,5),col=brewer.pal(length(net.ses),'Set1')[1])
abline(lm(log(pf.hist$counts+1)~log(pf.hist$mids)),col=brewer.pal(length(net.ses),'Set1')[1])
axis(1,col='white',labels=TRUE)
axis(2,col='white',labels=TRUE)


points(log(sg.hist$counts+1)~log(sg.hist$mids),pch=18,col=brewer.pal(length(net.ses),'Set1')[2],cex=1.5)
abline(lm(log(sg.hist$counts+1)~log(sg.hist$mids)),col=brewer.pal(length(net.ses),'Set1')[2])


points(log(se.hist$counts+1)~log(se.hist$mids),pch=17,cex=1.5,col=brewer.pal(length(net.ses),'Set1')[3])
abline(lm(log(se.hist$counts+1)~log(se.hist$mids)),lty=2,col=brewer.pal(length(net.ses),'Set1')[3])
dev.off()

summary(lm(log(pf.hist$counts+1)~log(pf.hist$mids)))
summary(lm(log(sg.hist$counts+1)~log(sg.hist$mids)))
summary(lm(log(se.hist$counts+1)~log(se.hist$mids)))


#Edges per node distribution
pf.epn=pf.sij
sg.epn=sg.sij
se.epn=se.sij

pf.epn=apply(pf.epn,2,mean)
sg.epn=apply(sg.epn,2,mean)
se.epn=apply(se.epn,2,mean)

hist(pf.epn)
hist(sg.epn)
hist(se.epn)


