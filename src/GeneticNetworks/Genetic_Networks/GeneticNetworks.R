#Genetic Networks
source('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/Genetic_Networks/ind_net.R')
library(igraph)
source('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/Genetic_Networks/igraph_adds.R')
library(RColorBrewer)


#Relict Data
#With parentals
dir('/Users/aeolus/Documents/Active_Projects/GeneticNetworks/Genetic_Networks')
x=read.csv('/Users/aeolus/Documents/Active_Projects/GeneticNetworks/Genetic_Networks/relicts data revised sw.csv')#[-303:-311,]
table(x[,2])
#removed the parentals
#x=x[x[,2]!='Fremont',]
#x=x[x[,2]!='narrowleaf',]
#x=x[x[,2]!='trichocarpa',]

#Change parental labels
levels(x[,2])[levels(x[,2])=='Fremont']='F'
levels(x[,2])[levels(x[,2])=='narrowleaf']='N'
levels(x[,2])[levels(x[,2])=='trichocarpa']='T'

#extract geographic and genotype removal info
NV.rm=x[,3]
NV.rm[is.na(NV.rm)]=0
x=x[NV.rm==0,]
NV.geo=x[,4:7]
x=x[,-3:-7]

#Take the most abundanct allele from all populations
x.=x[,-1:-2]
x.=allele.patch(x.)

x=data.frame(x[,1:2],x.)

#remove missing alleles (i.e., -1 values)
nr=nrow(x)

x[x==-1]<-NA
x=na.omit(x)
(nr-nrow(x))/nr
colnames(x)

#
y=x[,1:2]
x=x[,-1:-2]


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
pop=y[,2]
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
image(sij)

#Look for communities/modules
#http://cneurocvs.rmki.kfki.hu/igraph/doc/R/plot.common.html

graph.=as.matrix(round(sij,5))
graph.sij=graph.adjacency(graph.,mode='undirected',weighted=TRUE)
sij.wtc=walktrap.community(graph.sij,steps=4)
sij.mod=sij.wtc$membership+1

modularity(graph.sij,sij.mod,unclass(graph.sij)[[9]][4][[1]]$weight)

#par(bg='black')
g=graph.sij
vertex.label=rownames(graph.)
vertex.color=brewer.pal(length(unique(sij.mod)),'Set1')[sij.mod]
my.layout=layout.fruchterman.reingold(g)
vertex.label.color='white'
vertex.frame.color='white'
vertex.label.family='Helvetica'
edge.label=graph.[lower.tri(graph.)]
edge.label=edge.label[edge.label>0]
edge.label.color='violet'
edge.label.family='Helvetica'
edge.label.cex=0.8
edge.width=(edge.label*20)^1.5
edge.color='lightgrey'

plot.igraph(x=g,vertex.label=vertex.label,layout=my.layout,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label=round(edge.label,2),edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

plot.igraph(x=g,vertex.label=vertex.label,layout=my.layout,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

graph.=as.matrix(round(sij,1))
graph.sij=graph.adjacency(graph.,mode='undirected',weighted=TRUE)
g=graph.sij
edge.label=graph.[lower.tri(graph.)]
edge.label=edge.label[edge.label>0]
edge.label.cex=0.8
edge.width=(edge.label*20)^1.5

plot.igraph(x=g,vertex.label=vertex.label,layout=my.layout,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

my.layout=layout.drl(g)
plot.igraph(x=g,vertex.label=vertex.label,layout=my.layout,vertex.color=vertex.color,vertex.label.color=vertex.label.color,vertex.frame.color=vertex.frame.color,vertex.label.family=vertex.label.family,edge.label.color=edge.label.color,edge.label.cex=edge.label.cex,edge.label.family=edge.label.family,edge.width=edge.width,edge.color=edge.color)

#hand positioning
tk.layout=read.csv('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/Genetic_Networks/tklayout.csv')[,-1]
colnames(tk.layout)=c('','')
tkplot(g,layout=tk.layout,vertex.label=vertex.label,vertex.color=vertex.color,vertex.label.color='white',edge.width=edge.width)

#save the tk layout
#tk.layout=tkplot.getcoords(23)
#write.csv(tk.layout,file='tklayout.csv')

#Overlay graph with geographic info
library(maps)
library(gmaps)
library(shape)
attach(NV.geo)
sites=y[,2]
NV.geo=data.frame(sites,Center.Lat,Center.Long)
NV.geoX=data.frame(array(NA,c(length(unique(sites)),3)))
colnames(NV.geoX)=colnames(NV.geo)
NV.geoX[,1]=unique(sites)
for (i in seq(along=NV.geoX[,1])){
	q=NV.geo[NV.geo$sites==NV.geoX$sites[i],2:3]
	NV.geoX[i,2:3]=q[1,]
	}


map('state','nevada')
FNT=locator(3)
NV.geoX[1:3,2:3]=cbind(FNT$y,FNT$x)
text(NV.geoX[,3],NV.geoX[,2],labels=as.character(NV.geoX[,1]),col='black',cex=0.5)
FNT=locator(3)
NV.geoX[1:3,2:3]=cbind(FNT$y,FNT$x)
map('state','nevada')
text(NV.geoX[order(NV.geoX[,1]),3],NV.geoX[order(NV.geoX[,1]),2],labels=as.character(NV.geoX[order(NV.geoX[,1]),1]),col=brewer.pal(length(unique(sij.mod)),'Set1')[sij.mod],cex=0.5)



#Overall Fit Test
#From Dyer and Nasson 2004 using model deviance
#Dm=n_total*log(Sigma/S), where Sigma is the determinant if the MLE estimate of the covariance matrix and S is the determinant of the observed sample covariance matrix
#According to Fortuna et al 2009 and Dyer and Nasson 2004, the fit test isn't exactly necessary because adding additional links to improve fit doesn't necessarily change the patterns of network
#SKIP
#library(ggm)
#Cij=d2cov(dij)
#fitCovGraph(sij,abs(Cij),nrow(x))

#number of excluded edges
#sij.lower=sij[lower.tri(sij)]
#Dm.df=length(sij.lower[sij.lower==0])
#chisq.crit=qchisq(0.05,Dm.df,lower.tail=FALSE)
#plot(0:300,dchisq(0:300,Dm.df),type='l',xlab='Chi^2',ylab='Density')
#points(chisq.crit,dchisq(chisq.crit,Dm.df),pch=19,col='red')
#abline(v=chisq.crit,lty=2)
#text(chisq.crit+25,dchisq(chisq.crit,Dm.df),labels='P > 0.05',cex=0.85)