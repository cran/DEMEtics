Bootstrapping.D <- function (tab, bt=1000){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           tab <- all.pops.D;
#           allelefrequency,sample.sizes <- allelefreq()
#           HWE <- Hardy.Weinberg();

# Passed:
#           tab2,l -> Hardy.Weinberg();
#           tab3 -> allelefreq();
#           allelefrequency,sample.sizes -> D.calc();

# Output:
#           D.means, D.locis, significance.levels -> Workspace;
#           names(tab2) -> screen (print);
#------------------------------------------------------------------------------------------------------------------------------  
                                  
          # The input for this function has to be a table of the following format:
          
          #     individual population      fragment.length   locus
          # 1        B1.1     Borkum            323          L12
          # 2        B1.1     Borkum            266          L12
          # 3        B1.2     Borkum            325          L12
          # 4        B1.2     Borkum            274          L12
          # 5        B1.3     Borkum            266          L12
          # 6        B1.3     Borkum            323          L12
          # 7        B1.4     Borkum            325          L12 
          
          # The column names must be equal to those in this example. 
          # Certain tables of another format can be converted to this format
          # by the function 'inputformat()' that is included in this package.

          # bt defines the times of bootstrap-resampling.
allelefreq(tab)  
  
Empirical.values <- D.calc(allelefrequency,sample.sizes)

          # D.Chao values are calculated for each locus as well as
          # averaged over all loci from the empirical data table.

D.locus.empirical <- D.values[[1]]
D.means.empirical <- D.values[[2]]  
  

D.locus<-vector("list",length=bt)

          # This vector will be filled with the per locus calculated D values. 
          
D.means<-vector(length=bt)

          # This vector will be filled with the mean D values (mean over all
          # loci).       

tab2 <- split(tab,tab$locus)

          # The table is split according to the several loci.
          
number.loci <- length(tab2)

          # The number of loci that have been examined.           

HWEs <- numeric(0)

          # This vector will be filled with the logical values (TRUE, FALSE) that
          # give the information if all populations are in Hardy Weinberg Equilibrium
          # for a specific locus.

for (l in 1:number.loci) {

                    # For each locus, it has to be found out if all populations are in
                    # Hardy Weinberg Equilibrium.
          
                          Hardy.Weinberg(tab2,l)
                          
                                    # It is tested, if all populations are in Hardy Weinberg equilibrium
                                    # for the actual locus.
                                    # The result is either HWE=TRUE or HWE=FALSE.                 

                           HWEs <- c(HWEs,HWE)  
                           
                                    # The results for the several loci are combined in a single vector.
                                    
                           if (HWE==TRUE){
                                          cat("\n","All of these populations are in Hardy Weinberg Equilibrium with regard to the locus: ",names(tab2)[l],"\n",sep="")
                                          cat("Therefore, alleles are permuted among these populations for this locus.","\n","\n")                                      
                                          }else{
                                          cat("\n","Not all of these populations are in Hardy Weinberg Equilibrium with regard to the locus: ",names(tab2)[l],"\n",sep="")
                                          cat("Therefore, genotypes are permuted among these populations for this locus.","\n","\n")                                      
                                          }

                                                    # User information about the permutation method and its reasons.                                          
                                                                     
                           }        
                           
HWEs <- as.logical(HWEs)

          # zeros are transformed to FALSE, ones to TRUE.  
          
cat("\n","The bootstrapping process takes several minutes...","\n","\n")                                      

          # User information.                                                                                        
   
for (repetition in 1:bt){
              
                            tab2 <- split(tab,tab$locus)

                                      # The table is split according to the several loci.
                                      
                            number.loci <- length(tab2)

                                      # The number of loci that have been examined.           
                                      
                            tab3<-numeric(0)
                       
                            for (l in 1:number.loci){

                                      # The following commands are carried out for each locus separately.
                                      
                                                      if (HWEs[l]==TRUE){
                          
                                                                # The confidence limits of the measure of genetic distance for the
                                                                # several loci and over all loci are determined by a thousandfold
                                                                # resampling of the alleles (for each locus
                                                                # and all populations) if the populations are in Hardy Weinberg equilibrium.
                                                                # Literature: Goudet J, Raymond M, deMeeus T, Rousset F. (1996).
                                                                # Testing differentiation in diploid populations. Genetics 144,1933-1940.
                                                                                                
                                                                            allelepool<-as.numeric(as.vector(tab2[[l]]$fragment.length))

                                                                                      # The alleleles that have been found at a locus in all the populations
                                                                                      # are collected in a common vector named 'allelepool'.
                                                                            
                                                                            resampled.allelepool<-sample(allelepool,length(allelepool),replace=TRUE)

                                                                                      # The alleles from one locus that were found in all populations that
                                                                                      # were sampled, are resampled with replacement.

                                                                            tab2[[l]]$fragment.length<-resampled.allelepool

                                                                                      # The alleles for the actual locus are replaced with alleles from
                                                                                      # the resampled.allelepool.

                                                                            tab3<-rbind(tab3,tab2[[l]])

                                                                                      # The tables with the per locus resampled alleles are bound together
                                                                                      # to a single data frame.
                                                                            }else{
                                          
                                                                                     # The confidence limits of the measure of genetic distance for the
                                                                                     # several loci and over all loci are determined by a thousandfold
                                                                                     # resampling of the genotypes (for each locus and all populations) 
                                                                                     # if the populations are not in Hardy Weinberg equilibrium.
                                                                                     # Literature: Goudet J, Raymond M, deMeeus T, Rousset F. (1996).
                                                                                     # Testing differentiation in diploid populations. Genetics 144,1933-1940.
                                                                                                                                
                                                                                   tab2.genotype <- split(tab2[[l]]$fragment.length,as.vector(tab2[[l]]$individual))

                                                                                              # The genotypes found for the actual locus are filtered out of
                                                                                              # table2.  They are now represented according to the frequency
                                                                                              # with which they occured in the empirical data.

                                                                                   number.genotypes <- length(tab2.genotype)

                                                                                              # The number of genotypes that have to be resampled.

                                                                                   genotypepool<-as.data.frame(as.matrix(tab2.genotype))[1:number.genotypes,]

                                                                                              # The genotypes that have been found for a locus in all the populations
                                                                                              # are collected in a common vector named 'allelepool'.

                                                                                   resampled.genotypepool<-sample(genotypepool,number.genotypes,replace=TRUE)

                                                                                              # The genotypes of the actual locus that were found in all the sampled
                                                                                              # populations, are resampled with replacement.

                                                                                   resampled.genotypepool <- as.data.frame(resampled.genotypepool)

                                                                                              # The list is converted into a data frame format.

                                                                                   resampled <- numeric(0)

                                                                                   for (g in 1:number.genotypes){

                                                                                              # All the resampled genotypes will be combined in a single vector,
                                                                                              # where the two allele lengths of one genotype are placed together one
                                                                                              # under the other.

                                                                                                                  resampled <- c(resampled,as.numeric(as.vector(resampled.genotypepool[1:2,g])))
                                                                                                                  
                                                                                                                  }

                                                                                   tab2[[l]]$fragment.length <- resampled

                                                                                              # The genotypes for the actual locus are replaced with the genotypes from
                                                                                              # the resampled genotypepool.

                                                                                   tab3<-rbind(tab3,tab2[[l]])

                                                                                              # The tables with the per locus resampled alleles are bound together
                                                                                              # to a single data frame.
          
                                                                                  }
                                                                                                                       
                                                      }
    
                            allelefreq(tab3)
                            
                                      # The table that contains the allelefrequencies and the table that
                                      # lists the sample sizes are placed in the R workspace in the 
                                      # object List, but also separately in the object allelefrequency
                                      # and the object sample.sizes by this function.                              
                            
                            D.calc(allelefrequency,sample.sizes)
                            
                                      # The D values for the several loci and over all loci are
                                      # calculated for all the populations that have been examined.
                                      # The results are available from the object 'D.values'. 
                            
                            D.locus[[repetition]]<-D.values[[1]]
                            D.means[repetition]<-D.values[[2]]
                            
                            }

assign("D.means",D.means,pos = ".GlobalEnv")

          # The list "D.locus" and the vector "D.means" are both saved
          # according to their names in the workspace of R in order to be
          # available for further calculations.

critical.values.means <-

  cbind(D.means.empirical-1.96*(sd(as.numeric(as.vector(D.means)))),D.means.empirical+1.96*(sd(as.numeric(as.vector(D.means)))))

colnames(critical.values.means)<- c("0.95.conf.int.lower","0.95.conf.int.upper")

          # The standard bootstrap confidence intervals are calculated.
          
D.loci <- numeric(0)

for (n in 1:bt){
                  D.loci<-rbind(D.loci,D.locus[[n]])
                  }

          # The D values for the several loci are combined in a single
          # data frame.
          
assign("D.loci",D.loci,pos = ".GlobalEnv") 

          # The data.frame "D.loci" is saved with its according name
          # in the workspace of R in order to be available for further calculations.         
          
D.loci2<-split(D.loci,D.loci$locus)     

          # This data frame is splitted so that the data that belong to the same
          # locus are separated from those that belong to a different locus.
          
critical.value.loci<-numeric(0)               

for (l in 1: number.loci){

          # The critical  value is calculated separately for every locus.
          
                          critical.value.loci <-
                             rbind(critical.value.loci,
                                cbind(as.numeric(as.vector(D.locus.empirical$D[l]))-1.96*(sd(as.numeric(as.vector(D.loci2[[l]]$D)))),
                                      as.numeric(as.vector(D.locus.empirical$D[l]))+1.96*(sd(as.numeric(as.vector(D.loci2[[l]]$D))))))


                                    
                                    # The critical values for the several loci are combined together with
                                    # the actual loci names.
                          }

critical.value.loci<-as.data.frame(critical.value.loci)
colnames(critical.value.loci)<-c("0.95.conf.int.lower","0.95.conf.int.upper")


confidence.limits<-list(critical.value.loci,critical.values.means)
names(confidence.limits)<-c("confidence.limits.per.locus","confidence.limit.over.all.loci")
invisible(confidence.limits)
assign("confidence.limits",confidence.limits,pos = ".GlobalEnv")

          # The end result of the bootstrapping is the list called 
          # 'confidence.limits' that contains the significance levels
          # per locus and over all loci.
          # This list is printed and assigned in the workspace  to be available 
          # for further calculations by the name 'confidence.limits'.

}




