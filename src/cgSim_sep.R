cgSim <- function(trees='simulated trees',insects='simulated insects',
                  z='selection strength',VeT='tree trait variance',
                  Ve='insect environmental variance',VeN='interaction environmental variance',
                  K='insect carrying capacity',artpop.only=FALSE){

###Default Parameter Settings
  ##environmental variance in tree trait influences tree heritability (2 for high H2 and 8 for low H2)
  if (any(VeT=='tree trait variance')){VeT <- 2}
  ##selection strength = c from gamma^c
  if (any(z=='selection strength')){z <- 0}
  ##environmental variance for insect trait
  if (any(Ve=='insect environmental variance')){Ve <- 0.1}
  ##step size for environmental variance in interactions (0 15 30 45 60)
  if (any(VeN=='interaction environmental variance')){VeN <- 0}
  ##insect pop carrying capacity
  if (any(K=='insect carrying capacity')){K <- 100}
  ##negative abundance correction factor
  if (any(cf=='correction factor')){cf <- 0}
  ##trees
  if (any(trees=='simulated trees')){
    tree.gpm <- gpmTrees()
    trees <- simTrees(tree.gpm,VeT=VeT)
  }
  if (any(insects=='simulated insects')){
    insects <- simSpp()
  }

###Initiate output objects
T <- length(trees) #number of trees
I <- nrow(insects) #number of insects
art_g <- matrix(NA,nrow=T,ncol=I) #insect gene frequency
art_z <- matrix(NA,nrow=T,ncol=I) #insect trait mean
art_Vg <- matrix(NA,nrow=T,ncol=I) #insect gene variability
art_Vz <- matrix(NA,nrow=T,ncol=I) #insect trait variability
gen_load <- matrix(NA,nrow=T,ncol=I) #insect genetic load
dem_load <- matrix(NA,nrow=T,ncol=I) #insect demographic load
En <- matrix(NA,nrow=T,ncol=I) #Environmental effect
#gamma = selection intensity
art_pop <- matrix(NA,nrow=T,ncol=I) #insect population community matrix
                                        #names the species
colnames(art_pop) <- paste('S',I(1:I),sep='')
colnames(art_g) <- colnames(art_z) <- paste('S',I(1:I),sep='')
colnames(art_Vg) <- colnames(art_Vz) <- paste('S',I(1:I),sep='')
colnames(gen_load) <- colnames(dem_load) <- paste('S',I(1:I),sep='')
colnames(En) <- paste('S',I(1:I),sep='')

for (i in 1:T){
                                        #insects on trees
                                        #for each tree i
  for (j in 1:I){
                                        #for each insect j
###Equation 6 from MS - Arthropod (art) Gene frequency (art_g = pij)
    if (trees[i] < 2*insects[j,1]){
      art_g[i,j] <- 0
    }else if (trees[i] > 2*insects[j,2]){
      art_g[i,j] <- 1
    }else{
      art_g[i,j] <- (trees[i] - 2*insects[j,1]) / (2*insects[j,2] - 2*insects[j,1])
    }
###Equation 6 from MS - calculating mean trait Z from art gene frequency
###CHANGE = actually Equation 4
###This is the trait value at equilibrium (zbar^star_ij)
    art_z[i,j] <- 2*insects[j,2]*art_g[i,j]^2 + 
      2*art_g[i,j]*(1-art_g[i,j])*(insects[j,1]+insects[j,2]) + 
        2*insects[j,1]*(1-art_g[i,j])^2 + runif(1)*Ve - Ve/2
    ##Note: the final term (runif(1)*Ve - Ve/2) is the from the e_z term in Eq. 3.

###Art genetic (Vg = sigma^2_Gij) and trait variance (Vz = sigma^2_zbar_ij)
    art_Vg[i,j] <- 2*art_g[i,j]*(1-art_g[i,j])
    art_Vz[i,j] <- art_Vg[i,j]*(insects[j,2] - insects[j,1])^2 + Ve

###Evolutionary (gen) and demographic (dem) loads from selection
### dem_load[i,j] <- 0.5*(selection intensity)*(arthropod trait variance)
### gen_load = 0.5*(selection intensity)*(arthropod mean trait value - tree trait value)^2
    dem_load[i,j] <- 0.5*(0.00007924*2.511886^(z-1))*(art_Vz[i,j])
    gen_load[i,j] <- 0.5*(0.00007924*2.511886^(z-1))*(art_z[i,j] - trees[i])^2
    gamma <- (0.00007924*2.511886^(z-1))
###Equation 7 from MS - art predicted population size as a function of loads and ecological variance
###art_pop = n^star_ij = population of species j and tree i at equilibrium
###art_pop = (K = carrying capacity) * (community selection) + (E_n_ij = environmental variance)
    En[i,j] <- (runif(1)*VeN-VeN/2)
    art_pop[i,j] <- K * (1 - gen_load[i,j] - dem_load[i,j]) + En[i,j]
###preventing art pops from going negative
###NOTE: this code was changed from the original which made 
###zero or negative species into runif(1)*3)
    if (art_pop[i,j] < 0){art_pop[i,j] <- cf}
  } #end insect loop (j)
} #end tree loop (i)

  if (artpop.only){return(art_pop)}else{
    return(list(
                art_pop=art_pop,tree.gpm=tree.gpm,trees=trees,insects=insects,
                art_g=art_g,art_Vg=art_Vg,art_Vz=art_Vz,gamma=gamma,
                dem_load=dem_load,gen_load=gen_load,En=En
                )
           )
  }
}
