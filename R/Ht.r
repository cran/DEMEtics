Ht <- function(Table){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           Table <- Dest.calc(), Gst.est.calc();
#------------------------------------------------------------------------------------------------------------------------------
  

          # function to calculate the total heterozygosity.
    
          # The argument 'Table' has to be of the following format:

          #   allele number population   locus  proportion

          #       4      1     Borkum  LoPi89 0.1250000
          #       5      1     Borkum  LoPi89 0.1250000
          #       6      2     Borkum  LoPi89 0.2500000
          #       7      0     Borkum  LoPi89 0.0000000
          #       8      3     Borkum  LoPi89 0.3750000
          #       9      1     Borkum  LoPi89 0.1250000
          #       4      1   Langeoog  LoPi89 0.1666667
          #       5      1   Langeoog  LoPi89 0.1666667
          #       6      1   Langeoog  LoPi89 0.1666667
          #       7      2   Langeoog  LoPi89 0.3333333
          #       8      1   Langeoog  LoPi89 0.1666667
          #       9      0   Langeoog  LoPi89 0.0000000
          #       4      0 Wangerooge  LoPi89 0.0000000
          #       5      1 Wangerooge  LoPi89 0.1250000
          #       6      0 Wangerooge  LoPi89 0.0000000
          #       7      4 Wangerooge  LoPi89 0.5000000
          #       8      2 Wangerooge  LoPi89 0.2500000
          #       9      1 Wangerooge  LoPi89 0.1250000


          Table2 <- split(Table,Table$locus)

                    # The 'Table' is split according to the several loci that were
                    # examined.

          loci.numbers <- length(Table2)

                    # The number of loci that were examined.

          Ht.values <- numeric(0)

          for (l in 1:loci.numbers){

                    # For each locus, the Ht value will be calculated seperately.

                                      Table3 <- split(Table2[[l]],as.numeric(as.vector(Table2[[l]]$allele)))

                                                # The allelefrequency-values are splitted according to the several
                                                # alleles that were found in this population.

                                      number.alleles <- length(Table3)

                                                # The number of different alleles that were found in the actual locus.

                                      Means  <-  vector(length=number.alleles)

                                                # This vector will be filled with the mean allelefrequencies of one
                                                # allele over all populations that were examined.

                                      for (a in 1:number.alleles){

                                                # For every single allele of the actual locus the following commands are carried out
                                                # separately.

                                                                  Mean  <-  mean(as.numeric(as.vector(Table3[[a]]$proportion)))

                                                                            # The mean allelefrequency over all populations is calculated
                                                                            # for the actual allele.

                                                                  Means[a] <- Mean
                                                                  
                                                                  }
                                                                  
                                      Homozygotes  <-  Means^2

                                      Ht.one.locus <- 1-sum(Homozygotes)

                                             # The Ht value is calculated for the actual locus.
                          
                                      Ht.values <- rbind(Ht.values,(cbind(Ht.one.locus,names(Table2)[l])))
                          
                                             # The Ht values for the several loci are combined with a column that
                                             # contains the names of the loci, the Ht-value belongs to.
                          
                                      }

            Ht.values <- as.data.frame(Ht.values)
            colnames(Ht.values) <- c("Ht.value","locus")
            invisible(Ht.values)
  
}
