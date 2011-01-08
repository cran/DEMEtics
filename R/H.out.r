H.out <- function(tab,data.name=FALSE,filename){

  filename <- filename
          # A function that calculates the Hs Hs_est, Ht and Ht_est values for each locus separately.
          # The arguments of the function:

          # If data.name=TRUE, the name of the data file the result
          # table is saved at, can be exactly determined.  In this
          # case v is thre front part of the name and h the hind
          # part. The mean part will always be: H.values, so that the
          # file name will be as follows: v.H.values.h.txt. The parts
          # v and h have to be quoted, like v="FrontPart",
          # h="HindPart".

tab.pops <- split(tab,tab$population)

          # The table is splitted up into the data that belong to different
          # populations.           

allelefreq(tab)

          # From the table that contains the empirical data, the allelefrequencies
          # as well as the sample.sizes are calculated for each locus.
          # The table that contains the allelefrequencies and the table that
          # lists the sample sizes are assigned to the R workspace in the 
          # object List, but also separately in the object allelefrequency
          # and the object sample.sizes.  


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
          
          
                    # The data are ascribed to a data frame that is called Hj.values.

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

          Hs.est.values<-numeric(0)
          
                    # This vector will be filled with the per locus calculated Hs_est.values


          Ht.est.values<-numeric(0)
          
                    # This vector will be filled with the per locus calculated Ht_est.values

                       
          sample.sizes2<-split(sample.sizes,sample.sizes$locus)
          
                    # The sample size values are split according to the locus they belong
                    # to.

          # To be able to calculate Hs_est and Ht_est values per locus, some functions and sizes first have to be defined:

                    n<-length(sample.sizes)

                    # n is the number of populations that have been sampled.

                    # Function to calculate the harmonic mean of the sample.sizes.

                    harmonic<-function(sample.sizes){
                                      n/
                                      sum(1/sample.sizes)
                                      }
                                      

                    # Function to calculate Hs.est values:

                    Hs.est<-function(N,Hs){
                                      (2*N/(2*N-1))*Hs
                                      }

                    Ht.est<-function(N,n,Ht,Hs.est){
                                            Ht+ (Hs.est/(2*N*n))
                                            }


          for (l in 1:max(loci.numbers)){
                    
                    # For every locus, the following commands are carried out seperately.

                                        # The harmonic mean is calculated from the sample sizes...

                           N<-harmonic(sample.sizes2[[l]]$sample.size)
         
                                                        # ... and ascribed to the variable N.

                           
                                          Hs.actual<-as.numeric(as.vector(Hs.values$Hs.value[l]))
                                          Ht.actual<-as.numeric(as.vector(Ht.values$Ht.value[l]))
                                          sample.size<-as.numeric(as.vector(sample.sizes2[[l]]$sample.size))
                                          
                                          Hs.est.actual <- Hs.est(N,Hs.actual)
                                          Ht.est.actual <- Ht.est(N,n,Ht.actual,Hs.est.actual)
                                          
                                          Hs.est.values<-rbind(Hs.est.values,cbind(Hs.est.actual,names(sample.sizes2)[l]))
                                          Ht.est.values<-rbind(Ht.est.values,cbind(Ht.est.actual,names(sample.sizes2)[l]))

                                        # The Ht.est and Hs.est values
                                                    # for each locus
                                                    # are combined and
                                                    # the locus names
                                                    # are added to the
                                                    # values they
                                                    # belong to

                                          }

          Hs.est.values<-as.data.frame(Hs.est.values)
          colnames(Hs.est.values)=c("Hs.est.value","locus")

          Ht.est.values<-as.data.frame(Ht.est.values)
          colnames(Ht.est.values)=c("Ht.est.value","locus")

          H.output <- cbind(as.character(as.vector(Ht.est.values$locus)),as.numeric(as.vector(Hs.values$Hs.value)),as.numeric(as.vector(Hs.est.values$Hs.est.value)),as.numeric(as.vector(Ht.values$Ht.value)),as.numeric(as.vector(Ht.est.values$Ht.est.value)))

          colnames(H.output)=c("locus","Hs","Hs.est","Ht","Ht.est")

          H.output <- as.data.frame(H.output)
  
cat("\n","++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n",sep="")
cat("+++++++++++++++++++++++++++++++++++++++++ HETEROZYGOSITIES +++++++++++++++++++++++++++++++++++++++++","\n") 
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n","\n")
  
print(H.output)

if (data.name==TRUE){

filename.output <- paste(v,"_H.values_",sep="")
filename.output <- paste(filename.output,h,sep="")
filename.output <- paste(filename.output,".txt",sep="")


}else{

filename.output <- paste("H.values",".txt",sep="")
filename.output <- paste(filename,".",filename.output,sep="")}

          # The filename under which the table that contains the Hs,
          # Hs.est, Ht and Ht.est values for all the loci seperately,
          # is created

write.table(as.data.frame(as.matrix(H.output)),file=filename.output, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

          # The table is saved in the working directory. Two or more
          # analysis that are carried under the same working
          # directory, are saved under the same file name. The data
          # are combined in the order as they have been analysed. No
          # data is lost due to be overwritten.

cat("----------------------------------------------------------------------------------------------------","\n")         
cat("\n","Hs, Hs.est, Ht and Ht.est values for each locus are saved in ","'",filename.output,"'","\n",sep="")


          # User information about the end date of the analysis and the filenames
          # under which the several tables have been saved.

}
