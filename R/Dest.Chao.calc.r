Dest.Chao.calc <- function(tab){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           tab <- pair.pops.Dest.Chao();

# Output:
#           D.Chao.values -> Workspace;
#------------------------------------------------------------------------------------------------------------------------------  

          # Estimator of genetic distance based on Chao et al (2008).
          # Cited in Jost L. (2008). Gst and its relatives do not measure differentiation. 
          # Molecular Ecology 17,4015-4026.

          # tab is a table of the following format:
          
          #     individual population      fragment.length   locus
          # 1        B1.1     Borkum            323          L12
          # 2        B1.1     Borkum            266          L12
          # 3        B1.2     Borkum            325          L12
          # 4        B1.2     Borkum            274          L12
          # 5        B1.3     Borkum            266          L12
          # 6        B1.3     Borkum            323          L12
          # 7        B1.4     Borkum            325          L12 
          
allelefrequency <- tab

          # The function 'allelefreq()' is needed.
          # The result is the table 'allelefrequency' and the table
          # 'sample sizes' that are both ascribed to the workspace.

allelefrequency.loci <- split(allelefrequency,allelefrequency$locus)

          # The data is subdivided in tables that list the values for the
          # several loci.
          
number.loci <- length(allelefrequency.loci)

          # The number of loci that have been examined.
          
Dest.Chao.loci.values <- numeric(0)

          # This vector will be filled with the per locus calculated
          # Dest.Chao values.
          
for (l in 1:number.loci){

          # For each locus, the Dest.Chao value is calculated separately
          
                          allelefrequency.pop <- split(allelefrequency.loci[[l]],allelefrequency.loci[[l]]$population)
                          
                                    # The data is splitted according to the populations for the
                                    # actual locus.
                          
                          number.populations <- length(allelefrequency.pop)
                          
                                    # The number of populations that have been studied for the actual
                                    # locus.
                                    
                          alleles <- as.numeric(as.vector(allelefrequency.pop[[1]]$allele))
                          
                                    # The allele lengths that occur for the actual locus in any of the
                                    # populations.
   
                          a5 <- numeric(0)
                          b2 <- numeric(0)
          
                          for (i in alleles){
                          
                                    # For each allele occuring at least one of the populations, following commands
                                    # are carried out separately.

                                              a1.1 <- vector(length=number.populations)
                                              
                                                        # This vector will encompass a part of the values needed to
                                                        # calculate variable a.
                                                        
                                              a2.1 <- vector(length=number.populations)
                                              
                                                        # This vector will encompass another part of the values that are
                                                        # needed to calculate the variable a.
                                                        
                                              b1.1 <- vector(length=number.populations)
                                              
                                                        # This vector will be filled with the values that are needed for
                                                        # the calculation of the variable b.
          
                                              for (j in 1:number.populations){
                                              
                                                        # Following commands are carried out for the different populations
                                                        # separately.
          
                                                                                Nj <- sum(as.numeric(as.vector(allelefrequency.pop[[j]]$number)))

                                                                                          # The total number of alleles sampled from the actual locus and the
                                                                                          # actual population.
                                                      
                                                                                Nij <- as.numeric(as.vector(subset(allelefrequency.pop[[j]],allelefrequency.pop[[j]]$allele==i)$number))
                                            
                                                                                          # The number of samples of allele i sampled from subpopulation j
                                                                                
                                                                                a1.1[j] <- (Nij/Nj)
                                                                                a2.1[j] <- (Nij/Nj)^2
                                                      
                                                                                b1.1[j] <-(Nij*(Nij-1))/
                                                                                                        (Nj*(Nj-1))
                                                                                }
                                    
                                              a1 <- (sum(a1.1))^2
                                              a2 <- sum(a2.1)
                                              
                                              a3 <- a1-a2
                                        
                                              a4 <- a3/(number.populations-1)
                                              
                                              a5 <- c(a5,a4)
                                              
                                                  # The several values that are calculated for each allele separately
                                                  # are combined to a single vector.
                                                  
                                              b1 <- sum(b1.1)
                                              b2 <- c(b2,b1)
                                             
                                             }
                          a <- sum(a5)
                          b <- sum(b2)
                                    
                          Dest.Chao <- 1-(a/b)
                          Dest.Chao.loci.values <- rbind(Dest.Chao.loci.values,cbind(Dest.Chao,names(allelefrequency.loci)[l]))
                          
                          }
                          
Dest.Chao.loci.values <- as.data.frame(Dest.Chao.loci.values)
colnames(Dest.Chao.loci.values)=c("Dest.Chao","locus")

Dest.Chao.over.loci<-mean(as.numeric(as.vector((Dest.Chao.loci.values$Dest.Chao))))

D.Chao.values<-list(Dest.Chao.loci.values,Dest.Chao.over.loci)
names(D.Chao.values)=c("Dest.Chao.values.for.loci","Mean.Dest.Chao.value")
invisible(D.Chao.values)
assign("D.Chao.values",D.Chao.values,pos = ".GlobalEnv")

}
