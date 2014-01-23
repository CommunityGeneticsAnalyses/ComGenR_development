cgSim <-
function(tree.pheno='tree phenotype matrix',insect='community phenotype matrix',reps=10,GG=8,YY=5,VeT=8,Ve=0.1,VeN=15,K=100,cf=0){

if (any(tree.pheno=='tree phenotype matrix')){tree.pheno <- gpmTrees()}
if (any(insect=='community phenotype matrix')){insect <- gpmCom()}

###Experimental design
# reps <- 10 #number of times to run simulation
# GG <- 8 #number of selection scenarios
# YY <- 5 #number of environmental scenarios
###Parameter values
# VeT <- 8 #environmental variance in tree trait influences tree heritability (2 for high H2 and 8 for low H2)
# Ve <- 0.1 #environmental variance for insect trait
# VeN <- 15 #step size for environmental variance in interactions (0 15 30 45 60)
# K <- 100 #insect pop carrying capacity

###Initiate output objects
T <- nrow(tree.pheno) #number of trees
I <- nrow(insect) #number of insects
art_g <- matrix(NA,nrow=T,ncol=I) #insect gene frequency
art_z <- matrix(NA,nrow=T,ncol=I) #insect trait mean
art_Vg <- matrix(NA,nrow=T,ncol=I) #insect gene variability
art_Vz <- matrix(NA,nrow=T,ncol=I) #insect trait variability
gen_load <- matrix(NA,nrow=T,ncol=I) #insect genetic load
dem_load <- matrix(NA,nrow=T,ncol=I) #insect demographic load
Enij <- matrix(NA,nrow=T,ncol=I) #Environmental effect
gamma <- matrix(NA,nrow=T,ncol=I) #selection intensity
art_pop <- matrix(NA,nrow=T,ncol=I) #insect population community matrix
colnames(art_pop) <- paste('S',I(1:I),sep='')
                                        #STORAGE LIST
out <- gg.list <- yy.list <- list()
###Simulation
tic <- Sys.time() #start stopwatch
for (RR in 1:reps){
###generating trees to use
  scores_XX <- matrix(NA,nrow=nrow(tree.pheno),ncol=GG)
  scores_XX[,1] <- tree.pheno[,2] + runif(T,0,1) * VeT - VeT/2 #scores are in second column of trees
###filling phenotypic values for trees (same set of trees for all scenarios
  for (z in 2:GG){
    scores_XX[1:T,z] <- scores_XX[1:T,1]
  }
  trees <- scores_XX #just renaming phenotypes to be trees
  for (y in 1:YY){
                                        #YY variation scenarios of other ecological interactions
    for (z in 1:GG){
                                        #GG selection intensity scenarios
      for (i in 1:T){
                                        #insects on trees
                                        #for each tree i
        for (j in 1:I){
                                        #for each insect j
###Equation 6 from MS - Arthropod (art) Gene frequency (art_g = pij)
          if (trees[i,z] < 2*insect[j,1]){
            art_g[i,j] <- 0
          }else if (trees[i,z] > 2*insect[j,2]){
            art_g[i,j] <- 1
          }else{
            art_g[i,j] <- (trees[i,z] - 2*insect[j,1]) / (2*insect[j,2] - 2*insect[j,1])
          }
###Equation 6 from MS - calculating mean trait Z from art gene frequency
###CHANGE = actually Equation 4
###This is the trait value at equilibrium (zbar^star_ij)
          art_z[i,j] <- 2*insect[j,2]*art_g[i,j]^2 + 
            2*art_g[i,j]*(1-art_g[i,j])*(insect[j,1]+insect[j,2]) + 
              2*insect[j,1]*(1-art_g[i,j])^2 + 
                runif(1)*Ve - Ve/2
##Note: the final term (runif(1)*Ve - Ve/2) is the from the e_z term in Eq. 3.

###Art genetic (Vg = sigma^2_Gij) and trait variance (Vz = sigma^2_zbar_ij)
          art_Vg[i,j] <- 2*art_g[i,j]*(1-art_g[i,j])
          art_Vz[i,j] <- art_Vg[i,j]*(insect[j,2] - insect[j,1])^2 + Ve

###Evolutionary (gen) and demographic (dem) loads from selection
### dem_load[i,j] <- 0.5*(selection intensity)*(arthropod trait variance)
### gen_load = 0.5*(selection intensity)*(arthropod mean trait value - tree trait value)^2
          dem_load[i,j] <- 0.5*(0.00007924*2.511886^(z-1))*(art_Vz[i,j])
          gen_load[i,j] <- 0.5*(0.00007924*2.511886^(z-1))*(art_z[i,j] - trees[i,z])^2
          gamma[i,j] <- 0.5*(0.00007924*2.511886^(z-1))
###Equation 7 from MS - art predicted population size as a function of loads and ecological variance
###art_pop = n^star_ij = population of species j and tree i at equilibrium
###art_pop = (K = carrying capacity) * (community selection) + (E_n_ij = environmental variance)
          En[i,j] <- (runif(1)*VeN*(y-1)-VeN*(y-1)/2)
          art_pop[i,j] <- K * (1 - gen_load[i,j] - dem_load[i,j]) + Enij
###preventing art pops from going negative
###NOTE: this code was changed from the original which made 
###zero or negative species into runif(1)*3)
          if (art_pop[i,j] < 0){art_pop[i,j] <- cf}
        } #end insect loop (j)
      } #end tree loop (i)
###keeping track of which iteration the simulation is on reps*GG*(y-1)+reps*(z-1)+RR
###dlmwrite(int2str(reps*GG*(y-1)+reps*(z-1)+RR),art_pop,'\t');  %%% write arthropod sampling to file %%%
                                        #print(paste(reps*GG*(y-1),reps*(z-1),RR,sep='_'))
      print(paste(RR,y,z,sep=' '))
      gg.list[[z]] <- art_pop
                                        # names(gg.list)[z] <- paste(reps*GG*(y-1),reps*(z-1),RR,sep='')
      names(gg.list)[z] <- paste(RR,y,z,sep='_')
    } #end selection loop
    yy.list[[y]] <- gg.list
  } #end environment loop
  toc <- Sys.time() #stop stopwatch
  out[[RR]] <- yy.list
} #end REP loop and END simulation

return(out)
}
