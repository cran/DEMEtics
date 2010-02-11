Dest.calc <- function(input,sample.sizes){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           input, sample.sizes <- all.pops.Dest();
#           Hj.one.population <- Hj()
#           Hs.one.locus <- Hs()
#           Dest.one.locus <- Dest()
#           Ht.values <- Ht();

# Passed:
#           allelefrequency3[[pop]]$proportion -> Hj();
#           Hj.values2[[l]]$Hj.value -> Hs();
#           Hs.actual,Ht.actual,sample.size -> Dest();
#           allelefrequency -> Ht();

# Output:
#           D.values -> Workspace;
#------------------------------------------------------------------------------------------------------------------------------
  

          # A function that calculates the mean Dest values over all loci and
          # the Dest values for all loci separately.
          # The arguments of the function:
          # input = table with allelefrequencies as given by the function allelefreq().
          # sample.sizes = table with the sample sizes as given by the function allelefreq().

          allelefrequency <- input

          Hj.all.loci<-numeric(0)

                    # This vector will be filled with the Hj-values calculated for
                    # each population and all loci.

          Locus<-numeric(0)
          
                    # This vector will be filled with the Locus-names.

          allelefrequency2<-split(allelefrequency,allelefrequency$locus)

                    # The data.frame allelefrequency is splitted according to the loci
                    # that have been examined.
                    
          loci.numbers<-length(allelefrequency2)
          
                    # The number of different loci.

          Hj.values<-numeric(0)
          
                    # This vector will be filled with the Hj-values.

          for (l in 1:max(loci.numbers)){
          
                    # The following commands are carried out for each locus separately.

                                          allelefrequency3<-split(allelefrequency2[[l]],allelefrequency2[[l]]$population)
                                          
                                                    # For the actual locus, the table allelefrequency2 is splitted according
                                                    # to the several populations.
                                          
                                          number.populations<-length(allelefrequency3)
                                          
                                                    # The number of populations examined for the actual locus.

                                          for (pop in 1:number.populations){
                              
                                                    # The following commands are carried out for each population separately.
                                          
                                                                  Hj.one.population<-Hj(as.numeric(as.vector(allelefrequency3[[pop]]$proportion)))
                              
                                                                            # The Hj-value for the actual locus and population is calculated.
                                                                  
                                                                  Hj.values<-rbind(Hj.values,cbind(names(allelefrequency2[l]),names(allelefrequency3)[pop],Hj.one.population))
                              
                                                                            # The Hj-values are combined with a column of the names of the actual population
                                                                            # and a column of the names of the actual locus.
                                    
                                                                            }

                                          }

          Hj.values<-as.data.frame(Hj.values)
          colnames(Hj.values)=c("Locus","Population","Hj.value")
          
          
                    # The data are ascribed a data frame that is called Hj.values.

          Hj.values2<-split(Hj.values,Hj.values$Locus)
          
                    # The data.frame Hj.values is split in order to separate the data for
                    # the different loci.

          Hs.values<-numeric(0)
          
                    # This vector will be filled with the per locus calculated Hs values.
          
          for (l in 1:loci.numbers){
          
                    # The calculation of the Hs.value is carried out for every locus
                    # separately.
          
                                    Hs.one.locus<-Hs(as.numeric(as.vector(Hj.values2[[l]]$Hj.value)))
          
                                              # The Hs value is calculated for the actual locus from the Hj values
                                              # of all the populations for the actual locus.
          
                                    Hs.values<-rbind(Hs.values,(cbind(Hs.one.locus,names(Hj.values2[l]))))
                                    
                                              # The Hs values for the different loci are combined and a column
                                              # that contains the actual locus name, is added.
          
                                    }

          Hs.values<-as.data.frame(Hs.values)
          colnames(Hs.values)=c("Hs.value","locus")
          
                    # The Hs values are combined in a data frame and the columns are
                    # named.
          
          Ht.values<-as.data.frame(Ht(allelefrequency))
          
                    # The Ht.values are calculated.
          
          Dest.values<-numeric(0)
          
                    # This vector will be filled with the per locus calculated Dest
                    # values.
          
          sample.sizes2<-split(sample.sizes,sample.sizes$locus)
          
                    # The sample size values are split according to the locus they belong
                    # to.
          
          for (l in 1:max(loci.numbers)){
                    
                    # For every locus, the following commands are carried out seperately.

                                          Hs.actual<-as.numeric(as.vector(Hs.values$Hs.value[l]))
                                          Ht.actual<-as.numeric(as.vector(Ht.values$Ht.value[l]))
                                          sample.size<-as.numeric(as.vector(sample.sizes2[[l]]$sample.size))
                                          
                                          Dest.one.locus<-Dest(Hs.actual,Ht.actual,sample.size)
                                          
                                          Dest.values<-rbind(Dest.values,cbind(Dest.one.locus,names(sample.sizes2)[l]))
                                          
                                                    # The Dest values for each locus are combined and the locus names
                                                    # are added to the Dest value they belong to

                                          }
                                          
          Dest.values<-as.data.frame(Dest.values)
          colnames(Dest.values)=c("Dest","locus")
          
          Dest.over.loci<-mean(as.numeric(as.vector((Dest.values$Dest))))
          
          D.values<-list(Dest.values,Dest.over.loci)
          names(D.values)=c("Dest.values.for.loci","Mean.Dest.value")
          invisible(D.values)
          assign("D.values",D.values,pos = ".GlobalEnv")

}
