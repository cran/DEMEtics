
allelefreq <- function(tab){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           tab <- Dest.Chao(), all.pops.Dest(), all.pops.Gst();

# Output:
#           sample.sizes, allelefrequency, List -> Workspace;
#------------------------------------------------------------------------------------------------------------------------------
  
          # A function that calculates the allelefrequencies and sample sizes per
          # locus from a table of the following format:
          
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
     
y <- tab

          # the argument, a table,  is assigned to the object 'y'.

y2 <- split(y,y$locus)

          # The data are splitted according to the loci that have been
          # examined.

loci.numbers <- length(y2)

          # The number of loci are assigned to the object 'loci.numbers'.

sample.sizes <- numeric(0)

          # This vector will be filled with the sample sizes per locus and
          # population.

for (l in 1:loci.numbers){

          # For all the loci that have been examined, the following commands
          # are carried out separately.

                            y3 <- split(y2[[l]],y2[[l]]$population)

                                      # For the actual locus, the data that belong to different populations
                                      # are splitted.
                            number.populations <- length(y3)

                                      # The number of populations for which data of the actual locus are
                                      # available.
                            
                            y2.without.zeros <- y2[[l]][y2[[l]]$fragment.length!=0,]

                            x <- (table(y2.without.zeros[,2]))/2

                                      # The sample sizes for the actual locus and each population are 
                                      # assigned to the vector x.

                            tab <- (cbind(as.data.frame.table(x, row.names = NULL, optional = FALSE),names(y2[l])))
                            colnames(tab)=c("population","sample.size","locus")

                                      # The sample sizes for the actual locus and the corresponding populations
                                      # are combined in a data frame.

                            sample.sizes <- rbind(sample.sizes,tab)
                            
                                      # The data for the several loci are combined row by row.                           

                            allele.max <- max(as.numeric(as.vector(y2[[l]]$fragment.length)))
          
                                      # The largest allele for the actual locus is assigned to the object
                                      # 'allele.max'.

                            anzahlk <- numeric(0)

                                      # Besides further data, this object will be filled with the allefrequencies
                                      # for each allele of each locus. An allele that is set to zero (probably
                                      # because of amplification problems), will be ignored.

                            for (pop in 1:number.populations){
                                      
                                      # For all the populations and the actual locus, following commands are
                                      # carried out separately.
                                      
                                                                for (k in 1:allele.max){
          
                                                                          # The following analysis is carried out for each allele of the actual
                                                                          # locus and population separately. The reason why alleles that are set to
                                                                          # zero will be ignored (see above) is, that k starts with 1.

                                                                                        anzahlk <- rbind(anzahlk,cbind(k,length(which((as.numeric(as.vector(y3[[pop]]$fragment.length)))==k)),names(y3)[pop],names(y2)[l]))
          
                                                                                                  # The object 'anzahlk' is filled with the number of alleles of value
                                                                                                  # or length k, that occured in the actual population.

                                                                                        }

                                                                }
                                                                
                            anzahlk=as.data.frame(anzahlk)
                            colnames(anzahlk) <- c("allele","number","population","locus")

                                      # The object 'anzahlk' is set to a data frame format and the columns of this 
                                      # object are named.

                            for (al in 1:allele.max){
                            if (all((subset(anzahlk,allele==al)$number)==0)){
                                                                              anzahlk=anzahlk[-(which(anzahlk$allele==al)),]

                                                                              }
                                                    }
                                                          
                                      # Up to now, all the alleles from 1 to the allele with the maximal 
                                      # value for the actual locus have been taken into account. Those alleles, 
                                      # which values don't occur in any of the populations that have been 
                                      # sampled, will be removed from the object 'anzahlk'.

                            if (l<2) {
                                      anzahlkohnenull <- anzahlk
                                      } else
                                            anzahlkohnenull <- rbind(anzahlkohnenull,anzahlk)


                                      # If the object anzahlkohnenull is created for the first locus (l<2), then it is
                                      # defined as 'anzahlk'.
                                      # If however, the object 'anzahlkohnenull' has already been created 
                                      # for  the first locus, the object 'anzahlkohnenull' for the new locus 
                                      # will be attached to the object 'anzahlkohnenull' that already exists.

                                      # If another locus has been examined, the same commands are now carried
                                      # out for this one until the data for all the loci are present in the
                                      # object 'anzahlkohnenull' 

                            }

anzahlkohnenull=as.data.frame(anzahlkohnenull)

          # The object 'anzahlkohnenull' is set to the format of a data frame.

anzahlkohnenull2 <- split(anzahlkohnenull,anzahlkohnenull$locus)

          # The data for the different loci in the data frame 'anzahlkohnenull' 
          # are splitted.
          
proportion=numeric(0)

          # This vector will be filled with the frequencies with which the several 
          # alleles of a locus are present in each population.

for (l in 1:loci.numbers){

          # For each locus, the following commands are carried out separately.

                          anzahlkohnenull3 <- split(anzahlkohnenull2[[l]],anzahlkohnenull2[[l]]$population)

                                    # For the actual locus, the values that belong to different populations
                                    # are splitted.

                          number.populations <- length(anzahlkohnenull3)
      
                                    # The number of populations for which data for the 
                                    # actual locus are available.        

                          for (pop in 1:number.populations){

                                    # For the actual locus, the following commands are carried out seperately for
                                    # every population.
          
                                                                    allele.numbers=as.numeric(as.vector(anzahlkohnenull3[[pop]]$number))
                                                                    hundred.percent=sum(allele.numbers)
                                                                    percent=allele.numbers/hundred.percent
                                                                    proportion=c(proportion,percent)
                                                                    
                                                                    }
                          }

allelefrequency <- cbind(anzahlkohnenull,proportion)
colnames(allelefrequency)=names(allelefrequency)
rownames(allelefrequency)=seq(1:length(allelefrequency[,1]))   
         
List<-list(allelefrequency,sample.sizes)
names(List)=c("Table","Sample sizes")
invisible(List)

          # The matrix 'allelefrequency' is combined with the matrix 'sample.sizes'
          # in a list.
          
assign("sample.sizes",sample.sizes,pos = ".GlobalEnv")
assign("allelefrequency",allelefrequency,pos = ".GlobalEnv")
assign("List",List,pos = ".GlobalEnv")

          # The objects 'sample.sizes' and 'allelefrequency'are ascribed to the
          # workspace in order to be available for further calculations.

}
