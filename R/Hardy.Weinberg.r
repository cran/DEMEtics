Hardy.Weinberg <- function(tab2,l){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#                 tab2 <- Bootstrapping.Chao(), Bootstrapping(), Bootstrapping.Gst();

# Output:
#                 HWE -> Workspace;
#------------------------------------------------------------------------------------------------------------------------------  

          # This function calculates if all populations are in HWE for the
          # actual locus l.

tab2.pop <- split(tab2[[l]],tab2[[l]]$population)

          # The data are splitted so that those belonging to different populations
          # are separated.
          
number.pops <- length(tab2.pop)

          # The number of populations for which data for the actual locus were
          # obtained.
          

Hardy <- numeric(0)

          # This vector will be filled with the p.values for the several
          # populations that give the probabilities that the populations are
          # in HWE for the actual locus.
          
for (p in 1:number.pops){

                          tab2.pop.individual <- as.character(tab2.pop[[p]]$individual)

                                    # The individuals occuring in the actual population.
          
                          tab2.pop.alleles <- as.data.frame.table(table(as.numeric(as.vector(tab2.pop[[p]]$fragment.length))))
                          
                                    # The number of the several alleles of the actual locus occuring in the 
                                    # actual population.   
          
                          number.alleles <- tab2.pop.alleles[,2]      
                          all.alleles <- sum(tab2.pop.alleles[,2])
                          frequency.alleles <- number.alleles/all.alleles
                          
                                    # The frequency of the several alleles is calculated.
          
                          tab2.pop.alleles[,2] <- frequency.alleles
                          
                                    # The numbers of the different alleles are replaced by their 
                                    # frequencies with which they occur in the actual population.
          
                          tab2.pop.number.individuals <- length(tab2.pop.individual)/2

                                    # The number of individuals occuring in the actual population.

                          tab2.pop.genotypes <- split(as.data.frame(as.matrix(tab2.pop[[p]])),as.data.frame(as.matrix(tab2.pop.individual)))

                                    # The genotypes occuring for the actual locus and population.
                          
                                    # For each population it has to be found out seperately, if it is in
                                    # Hardy Weinberg equilibrium for the actual locus.

                          tab2.pop.number.alleles <- length(tab2.pop.alleles[,1])
                          
                                    # The number of alleles for the actual locus and population.
          
                          pairwise.combinations <- sum(seq(1:(tab2.pop.number.alleles)))

                                    # Depending on the number of alleles, the number of possible
                                    # combinations of alleles is calculated.
          
                                    # Construction for pairwise combination:

                          repetition <- seq(tab2.pop.number.alleles,1)
                          
                          allele.one <- numeric(0)

                          for (f in 1:(tab2.pop.number.alleles)){
                                                                  r <- rep(f,repetition[f])
                                                                  allele.one=c(allele.one,r)
                                                                  }

                                    # With these commands, the positions of the first allele in the
                                    # table tab2.pop.alleles that will be combined with another
                                    # allele, are created.
          
                          allele.two <- numeric(0)

                          for (g in 1:tab2.pop.number.alleles){
                                                              se <- seq(g,tab2.pop.number.alleles)
                                                              allele.two=c(allele.two,se)
                                                              }

                                    # The vector of the positions of the alleles in the table tab2.pop.alleles,
                                    # with which the first chosen alleles in the vector 'allele.one'
                                    # will be compared with, is created.
                                    
                          empirical.number.genotypes <- numeric(0)
                          
                                    # This vector will be filled with the number of occurences of the
                                    # several genotypes for the actual locus in the actual population.
                                    
                          estimated.frequency.genotypes <-  numeric(0)
                          
                                    # This vector will be filled with the estimated frequency of genotypes
                                    # that would be expected if the population would be in HWE for the
                                    # actual locus.
                                    
                          for (pc in 1:pairwise.combinations){

                                    # For each combination of alleles (or genotype), the following commands are
                                    # carried out seperately

                                                                first.allele <- as.numeric(as.vector(tab2.pop.alleles[allele.one[pc],1]))

                                                                          # The first allele of the actual genotype that is searched for
                                                                
                                                                second.allele <- as.numeric(as.vector(tab2.pop.alleles[allele.two[pc],1]))
                                                                
                                                                          # The second allele of the actual genotype that is searched for.
          
                                                                estimated.frequency.genotypes <- c(estimated.frequency.genotypes,ifelse(first.allele==second.allele,(tab2.pop.alleles[allele.one[pc],2])^2,2*(tab2.pop.alleles[allele.one[pc],2])*tab2.pop.alleles[allele.two[pc],2]))
          
                                                                          # The estimated frequencies of the different genotypes are calculated.
                                                                          # For homozygotes, p^2 is calculated, for heterozygotes 2*p*q.
                
          
                                                                this.genotype <- numeric(0)
                                                                
                                                                          # This vector will be filled with the number of individuals that possess
                                                                          # the actual genotype
                                                        
                                                                for (ind in 1:tab2.pop.number.individuals){
                                                                
                                                                          # Now, in every individual, the actual allele combination is searched.
                              
                                                                                                             occurrence <- as.numeric(as.vector(tab2.pop.genotypes[[ind]]$fragment.length[1]))==first.allele&&as.numeric(as.vector(tab2.pop.genotypes[[ind]]$fragment.length[2]))==second.allele ||
                                                                                                                           as.numeric(as.vector(tab2.pop.genotypes[[ind]]$fragment.length[2]))==first.allele&&as.numeric(as.vector(tab2.pop.genotypes[[ind]]$fragment.length[1]))==second.allele
                                                                                                             
                                                                                                                      # If the actual allele combination is represented in the actual individual,
                                                                                                                      # a TRUE is returned, otherwise a FALSE.
          
                                                                                                              this.genotype <- c(this.genotype,occurrence)
                                                                                                              
                                                                                                                      # This vector is filled with logical values about the occurence of the
                                                                                                                      # actual genotype.
          
                                                                                                              }
                                                                                                              
                                                                  empirical.number.genotypes <- c(empirical.number.genotypes,length(which(this.genotype==TRUE)))
                                                                  
                                                                            # The number of occurences of the several genotypes are combined
                                                                            # in this vector.
          
                                                                  } 
                                                                  
                          estimated.number.genotypes <- tab2.pop.number.individuals*estimated.frequency.genotypes
                          
                                    # The estimated number of genotypes, if the population were in 
                                    # HWE for the actual locus, is calculated.
          
                          if (length(estimated.number.genotypes)==1 || all(estimated.number.genotypes==estimated.number.genotypes[1]) ||  all(empirical.number.genotypes==empirical.number.genotypes[1])) {
                                                                      p.value <- 1
                                                                      } else {
                                                                      p.value <- chisq.test(estimated.number.genotypes,empirical.number.genotypes,simulate.p.value=TRUE,B=10000)$p.value
                                                                      }                                                                      
                                    # If a locus is made up of just one single allele, a chisquare test
                                    # can't be carried out. In this case, with p.value <- 1, the actual
                                    # population is set to be in Hardy Weinberg equilibrium for this locus.    
          
                          p.value <- ifelse (is.na(p.value), 0,p.value)       
                          
                                    # In the case that a Chisquare test can't be carried out, the randomization
                                    # of alleles is not justified. Therefore, in this case, the p.value is set
                                    # to 0 (not in HW-Equilibrium).                                                                              
          
                          Hardy <- c(Hardy,p.value)
                          
                                    # The p.values (probabilities if HWE=TRUE) are combined. These are
                                    # obtained from a 10.000 fold repeated Monte Carlo simulation.
          
                           }    
                                 
HWE <- ifelse(all(Hardy>0.05),TRUE,FALSE)        

          # If all populations are in HWE for this locus, HWE is set as TRUE,
          # otherwise it is set as FALSE. 
          
assign("HWE",HWE,pos = ".GlobalEnv")          

          # This end result is assigned to the workspace.
          
}                            
                                                                                    
                                                                                                              
                                                                                                                        


          
          
          
