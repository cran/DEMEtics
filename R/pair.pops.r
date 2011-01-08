pair.pops <- function (tab,statistics,bt,x,filename){

  filename <- filename

  # x defines whether D, Dest, Gst or Gst.est is calculated
# tab is the table
# statistics either has the value p, CI or both
# bt defines the number of bootstraps to be carried out

   
          # Function that calculates the Gst or D value for each locus separately and
          # the mean Gst.est value over all loci, for each possible pairwise population 
          # comparison.
allelefreq(tab)

          # The function 'allelefreq' must be present in the r workspace.
          # In this function, the object 'allelefrequency' and the object
          # 'sample sizes' is calculated.
          
y <- tab

          # The table is saved in the object y and needed to calculate the
          # Bootstrap-value.
          
y.pop=split(y,y$population)

          # The table is spltted according to the several populations
          # that have been examined.

allelefrequency.pop=split(allelefrequency,allelefrequency$population)

          # The allelefrequency is splitted in the several populations
          # that have been examined.

sample.sizes.pop=split(sample.sizes,sample.sizes$population)

          # The sample sizes are splitted for the several populations.

number.of.populations=length(allelefrequency.pop)

          # Number of populations that have been studied.

pairwise.comparisons=sum(seq(1:(number.of.populations-1)))

          # Depending on the number of populations, the number of possible
          # comparisons between populations is calculated.

          # Construction for pairwise comparison:

repetition=seq(number.of.populations-1,1)

population.one=numeric(0)

for (f in 1:(number.of.populations-1)){
                                     r=rep(f,repetition[f])
                                     population.one=c(population.one,r)
                                     }

          # With these commands, the positions of the first population in the
          # table allelefrequency.pop that will be compared with another
          # population, are created.

population.two=numeric(0)

for (g in 2:number.of.populations){
                                     se=seq(g,number.of.populations)
                                     population.two=c(population.two,se)
                                     }

          # The vector of the positions of the populations in the table allelefrequency.pop,
          # with which the first chosen population in the vector 'population.one'
          # will be compared with, is created.
          

##################### Calculations before bootstrapping

          
v.locis=numeric(0)

v.mean=numeric(0)

for (i in 1:pairwise.comparisons){

          # For every pairwise comparison, the following commands are carried out
          # seperately.
          
                                    v.loci=numeric(0)
                                    
                                              # This vector will be filled with the values for the several
                                              # comparisons.
                                    
                                    allelefrequency.pair=as.data.frame(as.matrix(rbind(allelefrequency.pop[[(population.one[i])]],allelefrequency.pop[[(population.two[i])]])))
                                    
                                              # The allelefrequencies of the populations that shall be compared
                                              # pairwise with one another, are selected from the table
                                              # 'allelefrequency.pop' and then combined to a table that is named
                                              # allelefrequency.pair.
                                              # In order to get rid of the population levels that are not included
                                              # in the actual comparison, the transformation to a matrix  and then
                                              # to a data frame is carried out.
                                              
                                    sample.sizes.pair=as.data.frame(as.matrix(rbind(sample.sizes.pop[[(population.one[i])]],sample.sizes.pop[[(population.two[i])]])))
                                                
                                              # The sample.sizes for the actual pair of populations that are
                                              # compared to one another, are combined.
                                              # In order to get rid of the population levels that are not included
                                              # in the actual comparison, the transformation to a matrix  and then
                                              # to a data frame is carried out.
                                              
                                    y.pair=as.data.frame(as.matrix(rbind(y.pop[[(population.one[i])]],y.pop[[(population.two[i])]])))
                                    
                                              # The raw data for the actual two populations are combined to a new
                                              # table.
                                              # In order to get rid of the population levels that are not included
                                              # in the actual comparison, the transformation to a matrix  and then
                                              # to a data frame is carried out.
                                    
                                    allelefrequency.pair2=split(allelefrequency.pair,allelefrequency.pair$locus)
                                    
                                              # The table 'allelefrequency.pair' is splitted to get separated
                                              # the values for the different loci. 
                                              
                                    allelefrequency.pair.pops=split(allelefrequency.pair,as.character(allelefrequency.pair$population))
                                    
                                              # The table 'allelefrequency.pair' is split according to the populations
                                              # that are actually compared.          
                                    
                                    number.loci=length(allelefrequency.pair2)
                                    
                                              # The number of loci that have been examined for the actual
                                              # populationpair
                                    
                                              
                                    calc(allelefrequency.pair,sample.sizes.pair,x)
                                    
                                              # The Gst or D values for every locus and the mean Gst or D value over all
                                              # loci is calculated.
                                              # The result is saved in a list called 'values'.          
                                    
                                    for (l in 1:number.loci){
                                    
                                                              allelefrequency.pair3=split(allelefrequency.pair2[[l]],as.character(allelefrequency.pair2[[l]]$population))
                                                    
                                                                        # The table 'allelefrequency.pair2 is split according to the populations
                                                                        # that are actually compared
                                                                        # It is defined again for the case that for one locus, data for one
                                                                        # of the populations would be lacking.
                                                              
                                                              v.loci=rbind(v.loci,cbind(as.numeric(as.vector(values[[1]][l,1])),
                                                              names(allelefrequency.pair2)[l],names(allelefrequency.pair3)[1],
                                                              names(allelefrequency.pair3)[2]))
                                                    
                                                                        # In this vector,the Gst or D value for
                                                                        # the actual locus, the actual locus, population one and populations
                                                                        # two that are actually compared with one another are combined.
                                                                                                            
                                                              }
                                    
                                     Result.actual.comparison.locis <- v.loci
                                      Result.actual.comparison.locis <- as.data.frame(Result.actual.comparison.locis)
                                      colnames(Result.actual.comparison.locis) <- c(paste(x,".for.locus",sep=""),"locus","population1","population2")            
                                    
                                    
                                    
                                    Result.actual.comparison.mean <-cbind(values[[2]],names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2])
                                      Result.actual.comparison.mean <-as.data.frame(Result.actual.comparison.mean)
                                      colnames(Result.actual.comparison.mean)<-c(paste("Mean.",x,sep=""),"population1","population2")
                                           
                                           
                                    interm.result <- list(Result.actual.comparison.locis,Result.actual.comparison.mean)
                                    names(interm.result) <- c(paste(x,".loci.pairwise.comparison",sep=""),paste(x,".mean.pairwise.comparison",sep=""))
                                    
                                    v.locis=rbind(v.locis,v.loci)
                                        
                                              # The data frames for
                                              # the pairwise
                                              # comparisons are
                                              # combined.
                                    
                                    v.mean=rbind(v.mean,cbind(values[[2]],names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2]))
                                    
                                              # In this vector, the Gst or D value over all loci, population one and
                                              # population two that are actually compared with one another and thes
                                              # two named as populationpair are combined.
                                                            
                                              
                                                            }

v.locis=as.data.frame(v.locis)
colnames(v.locis)=c(paste(x,".locus",sep=""),"Locus","Population1","Population2")            

v.mean=as.data.frame(v.mean)
colnames(v.mean)=c(paste(x,".mean",sep=""),"Population1","Population2")
       
v.pairwise=list(v.locis,v.mean)
names(v.pairwise)=c(paste(x,".loci.pairwise.comparison",sep=""),paste(x,".mean.pairwise.comparison",sep=""))
       

pairwise.adjusted <- v.pairwise



names(pairwise.adjusted)=c(paste(x,".loci.pairwise.comparison",sep=""),paste(x,".mean.pairwise.comparison",sep=""))

          
colnames(pairwise.adjusted[[1]]) <- c(paste(x,".for.locus",sep=""),"locus","population1","population2")  
       
colnames(pairwise.adjusted[[2]]) <- c(paste("Mean.",x,sep=""),"population1","population2")
       
rownames(pairwise.adjusted[[1]]) <- seq(1,length(pairwise.adjusted[[1]][,1]))            
rownames(pairwise.adjusted[[2]]) <- seq(1,length(pairwise.adjusted[[2]][,1]))  

cat("\n","++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n",sep="")
cat("+++++++++++++++++++++++++++++++++++  RESULTS WITHOUT STATISTICS ++++++++++++++++++++++++++++++++++++","\n")
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n","\n")         
 
print(pairwise.adjusted)
     
assign(paste(x,"pairwise.adjusted",sep=""),pairwise.adjusted,pos = ".GlobalEnv")

          # The list 'pairwise.adjusted' is assigned to the workspace
          # and therefore available for further calculations.
          
actual.date <- as.Date(Sys.time())
filename.for.loci <- paste("pairwise.",x,".loci.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

filename.for.mean <- paste("pairwise.",x,".mean.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")

                 
write.table(as.data.frame(as.matrix(pairwise.adjusted[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(pairwise.adjusted[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
cat("----------------------------------------------------------------------------------------------------","\n")         
cat("The ",x,".mean values were calculated as the arithmetic means of the ",x,".loci values","\n","for the respective pairwise comparison","\n",sep="")


cat("\n",x, " values for pairwise comparisons between populations for every locus are saved in ","\n","'",filename.for.loci,"'","\n",sep="")
cat("\n",x, " values averaged over all loci for pairwise comparisons between populations saved in ","\n","'",filename.for.mean,"'","\n",sep="")


############### Calculations with Bootstrapping
if (statistics=="p"||statistics=="all"){

v.locis=numeric(0)

v.mean=numeric(0)

for (i in 1:pairwise.comparisons){

          # For every pairwise comparison, the following commands are carried out
          # seperately.
          
                                    v.loci=numeric(0)
                                    
                                              # This vector will be filled with the values for the several
                                              # comparisons.
                                    
                                    allelefrequency.pair=as.data.frame(as.matrix(rbind(allelefrequency.pop[[(population.one[i])]],allelefrequency.pop[[(population.two[i])]])))
                                    
                                              # The allelefrequencies of the populations that shall be compared
                                              # pairwise with one another, are selected from the table
                                              # 'allelefrequency.pop' and then combined to a table that is named
                                              # allelefrequency.pair.
                                              # In order to get rid of the population levels that are not included
                                              # in the actual comparison, the transformation to a matrix  and then
                                              # to a data frame is carried out.
                                              
                                    sample.sizes.pair=as.data.frame(as.matrix(rbind(sample.sizes.pop[[(population.one[i])]],sample.sizes.pop[[(population.two[i])]])))
                                                
                                              # The sample.sizes for the actual pair of populations that are
                                              # compared to one another, are combined.
                                              # In order to get rid of the population levels that are not included
                                              # in the actual comparison, the transformation to a matrix  and then
                                              # to a data frame is carried out.
                                              
                                    y.pair=as.data.frame(as.matrix(rbind(y.pop[[(population.one[i])]],y.pop[[(population.two[i])]])))
                                    
                                              # The raw data for the actual two populations are combined to a new
                                              # table.
                                              # In order to get rid of the population levels that are not included
                                              # in the actual comparison, the transformation to a matrix  and then
                                              # to a data frame is carried out.
                                    
                                    allelefrequency.pair2=split(allelefrequency.pair,allelefrequency.pair$locus)
                                    
                                              # The table 'allelefrequency.pair' is splitted to get separated
                                              # the values for the different loci. 
                                              
                                    allelefrequency.pair.pops=split(allelefrequency.pair,as.character(allelefrequency.pair$population))
                                    
                                              # The table 'allelefrequency.pair' is split according to the populations
                                              # that are actually compared.          
                                    
                                    number.loci=length(allelefrequency.pair2)
                                    
                                              # The number of loci that have been examined for the actual
                                              # populationpair
                                    
                                    time1 <- Sys.time()
                                    
                                              # Time before bootstrapping.


cat("\n","====================================================================================================","\n",sep="")
cat("=============================== Bootstrapping for P-value Calculation ==============================","\n")
cat("====================================================================================================","\n","\n")
cat("\n","Pairwise comparison:",names(allelefrequency.pair.pops),"\n",sep=" ")
print(paste("Start of analysis: ",Sys.time(),sep=""))                                    

cat("\n","\n","WARNING: Depending on the size of your input data, the performance of your computer and the number of ","\n","bootstrap resamplings you have chosen, bootstrapping can take very long (hours to days).","\n",sep="")
                                    

                                    
                                    Bootstrapping.p(y.pair,bt,x)
                                    
                                              # Confidence levels for the Gst or D values for the actual pair of
                                              # populations is returned in the object 'confidence.limits'.
                                              # The bootstrap values are available from the objects 'loci'
                                              # and 'means'.
                                              
                                    calc(allelefrequency.pair,sample.sizes.pair,x)
                                    
                                              # The Gst or D values for every locus and the mean Gst.est value over all
                                              # loci is calculated.
                                              # The result is saved in a list called 'values'.          
                                    
                                    for (l in 1:number.loci){
                                    
                                                              allelefrequency.pair3=split(allelefrequency.pair2[[l]],as.character(allelefrequency.pair2[[l]]$population))
                                                    
                                                                        # The table 'allelefrequency.pair2 is split according to the populations
                                                                        # that are actually compared
                                                                        # It is defined again for the case that for one locus, data for one
                                                                        # of the populations would be lacking.
                                                              
                                                              v.loci=rbind(v.loci,cbind(as.numeric(as.vector(values[[1]][l,1])),
                                                              names(allelefrequency.pair2)[l],names(allelefrequency.pair3)[1],
                                                              names(allelefrequency.pair3)[2]))
                                                    
                                                                        # In this vector,the Gst or D  value for
                                                                        # the actual locus, the actual locus, population one and populations
                                                                        # two that are actually compared with one another are combined.
                                    
                                                              }
                                    
                                              # now, the according p-values are calculated and added to the tables.
                                    
                                    loci2 <- split(loci,loci$locus)
                                    
                                              # The bootstrap values are separated in different tables according
                                              # to the loci they belong to.         
                                    
                                    number.loci=length(loci2)
                                    
                                              # The number of loci that have been studied.
                                    
                                    p.values.loci=numeric(0)
                                    
                                              # This vector will be filled with as many p.values as loci have
                                              # been examined.
                                              
                                    time2 <- Sys.time() 
                                    
                                              # Time after bootstrapping
                                              
                                    bootstrap.time <- difftime(time2,time1,units="secs")   
                                    bootstrap.time.output <- round((as.numeric(bootstrap.time/60)),2)     


cat("\n","... Bootstrapping for p-value calculation (comparison: ",names(allelefrequency.pair.pops),") is terminated.","\n",sep="")
print(paste("Duration of the bootstrapping analysis for this comparison (min):",bootstrap.time.output,sep=" "))
cat("\n","\n","Out of ",pairwise.comparisons," pairwise comparisons, ",i," have already been analysed.","\n",sep="")
print(paste("Estimated end of the whole analysis for p-value calculation:",((pairwise.comparisons-i)*(as.numeric(bootstrap.time))+Sys.time()),sep=" ")) 
                                                                      
                                              # The estimated end of the whole analysis  is as many times the bootstrap.time
                                              # as comparisons still have to be carried out.       
                                              
                                    for (l in 1:number.loci){
                                    
                                              # The following commands are carried out separately for the different
                                              # loci.
                                              
                                                              bootstrapped.values.loci=as.numeric(as.vector(loci2[[l]][,1]))
                                                              empirical.value.loci=as.numeric(as.vector(values[[1]][l,1]))
                                                              
                                                              if (is.nan(empirical.value.loci)){
                                                                                                    p.value <- NA
                                                                                                  } else 
                                                                                                        p.val(empirical.value.loci,bootstrapped.values.loci)
                                                              
                                                                                                                # This function gives the p-value for the actual empirical value
                                                                                                                # in the object 'p.value'.
                                                                        
                                                                        # In the case that no bootstrapped values had been obtained
                                                                        # when populations that are compared contain all
                                                                        # the same allele at one locus (allelefrequency 100%), a Gst or D value
                                                                        # can't be calculated; (0/0) is calculated.
                                                                        # The empirical value is also NA in this case and then the p.value
                                                                        # is also given as NA. Otherwise a p.value is calculated from the
                                                                        # bootstrapped values.
                                                              
                                                              p.values.loci=rbind(p.values.loci,cbind(round(p.value,4),names(loci2)[l]))
                                                             
                                                              }
                                    
                                    bootstrapped.values=means
                                    empirical.value=values[[2]]
                                    
                                    p.value.over.all=p.val(empirical.value,bootstrapped.values)
                                    p.value.over.all=round(p.value,4)
                                    
                                    
                                    
                                    Result.actual.comparison.locis <- cbind(v.loci,as.numeric(as.vector(p.values.loci[,1])))
                                      Result.actual.comparison.locis <- as.data.frame(Result.actual.comparison.locis)
                                      colnames(Result.actual.comparison.locis) <- c(paste(x,".locus",sep=""),"Locus","Population1","Population2","P.value")            
                                    
                                    
                                    
                                    Result.actual.comparison.mean <- cbind(values[[2]],names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],p.value.over.all)
                                      Result.actual.comparison.mean <- as.data.frame(Result.actual.comparison.mean)
                                      colnames(Result.actual.comparison.mean) <- c(paste(x,".mean",sep=""),"Population1","Population2","P.value")  
                                           
                                           
                                    actual.date <- as.Date(Sys.time())
                                    filename.for.loci <- paste("pairwise.",x,".loci.p.",actual.date,sep="")
                                    filename.for.loci <- paste(filename.for.loci,".txt",sep="")
                                    filename.for.loci <- paste(filename,".",filename.for.loci,sep="")
                                    
                                    filename.for.mean <- paste("pairwise.",x,".mean.p.",actual.date,sep="")
                                    filename.for.mean <- paste(filename.for.mean,".txt",sep="")
                                    filename.for.mean <- paste(filename,".",filename.for.mean,sep="")
                                    
                                    write.table(Result.actual.comparison.locis,file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)         
                                    write.table(Result.actual.comparison.mean,file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)         
                                    
                                    interm.result <- list(Result.actual.comparison.locis,Result.actual.comparison.mean)
                                    names(interm.result) <- c(paste(x,".loci.pairwise.comparison",sep=""),paste(x,".mean.pairwise.comparison",sep=""))
                                    
cat("\n","++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n",sep="")
cat("++++++++++++++++++++++++++++++++ INTERMEDIATE RESULTS WITH P-VALUES ++++++++++++++++++++++++++++++++","\n")
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n","\n")    
                                    
                                    print(interm.result)
cat("----------------------------------------------------------------------------------------------------","\n")
print(paste("This analysis was finished on",Sys.time()))
cat("The ",x,".mean value was calculated as the arithmetic mean of the ",x,".loci values","\n",sep="")              

                                    cat("\n",x," values for the comparison between ",names(allelefrequency.pair.pops)[[1]]," and ",names(allelefrequency.pair.pops)[[2]],"\n"," for every locus are saved in ","\n","'",filename.for.loci,"'","\n",sep="")
cat("\n",x," values averaged over loci for the comparison between ",names(allelefrequency.pair.pops)[[1]]," and ",names(allelefrequency.pair.pops)[[2]],"\n"," are saved in ","\n","'",filename.for.mean,"'","\n",sep="")
                                    
                                    v.locis=rbind(v.locis,cbind(v.loci,as.numeric(as.vector(p.values.loci[,1]))))
                                        
                                              # The according p-values are added to the data.frame v.loci and the
                                              # data frames for the pairwise comparisons are combined.
                                    
                                    v.mean=rbind(v.mean,cbind(values[[2]],names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],p.value.over.all))
                                    
                                              # In this vector, the
                                              # Gst or D value over
                                              # all loci, population
                                              # one and population two
                                              # that are actually
                                              # compared with one
                                              # another are combined.
                                              # The p-values are
                                              # added.
                                              
                                  }

v.locis=as.data.frame(v.locis)
colnames(v.locis)=c(paste(x,".for.locus",sep=""),"locus","population1","population2","p.values")            

v.mean=as.data.frame(v.mean)
colnames(v.mean)=c(paste("Mean.",x,sep=""),"population1","population2","p.values")
       
v.pairwise=list(v.locis,v.mean)
names(v.pairwise)=c(paste(x,".loci.pairwise.comparison",sep=""),paste(x,".mean.pairwise.comparison",sep=""))
       
p.value.correcture(v.pairwise)

          # A function that adjusts the p-values after Bonferroni, Holm, Hommel
          # and Benjamini and Hochberg.
          # The result is a list called 'Dv.pairwise.adjusted'.

pairwise.adjusted <- Dv.pairwise.adjusted



names(pairwise.adjusted)=c(paste(x,".loci.pairwise.comparison",sep=""),paste(x,".mean.pairwise.comparison",sep=""))

          
colnames(pairwise.adjusted[[1]]) <- c(paste(x,".for.locus",sep=""),"locus","population1","population2","p.values","p.bonferroni","p.holm","p.hommel","pBH")  
       
colnames(pairwise.adjusted[[2]]) <- c(paste("Mean.",x,sep=""),"population1","population2","p.values","p.bonferroni","p.holm","p.hommel","pBH")   
       
rownames(pairwise.adjusted[[1]]) <- seq(1,length(pairwise.adjusted[[1]][,1]))            
rownames(pairwise.adjusted[[2]]) <- seq(1,length(pairwise.adjusted[[2]][,1]))  


cat("\n","++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n",sep="")
cat("+++++++++++++++++++++++++++++++++++++ END RESULTS WITH P-VALUES ++++++++++++++++++++++++++++++++++++","\n")
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n","\n")    
       
print(pairwise.adjusted)
     
assign(paste(x,".pairwise.adjusted",sep=""),pairwise.adjusted,pos = ".GlobalEnv")

          # The list 'pairwise.adjusted' is assigned to the workspace
          # and therefore available for further calculations.
          
actual.date <- as.Date(Sys.time())
filename.for.loci <- paste("pairwise.",x,".loci.p.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

filename.for.mean <- paste("pairwise.",x,".mean.p.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")

                 
write.table(as.data.frame(as.matrix(pairwise.adjusted[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(pairwise.adjusted[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

cat("\n","----------------------------------------------------------------------------------------------------","\n","\n")
print(paste("This analysis was finished on",Sys.time()))
cat("The ",x,".mean values were calculated as the arithmetic means of the ",x,".loci values","\n","for the respective pairwise comparison","\n",sep="")
cat("\n",x, " values for pairwise comparisons between populations for every locus are saved in ","\n","'",filename.for.loci,"'","\n",sep="")
cat("\n",x, " values averaged over all loci for pairwise comparisons between populations are saved in ","\n","'",filename.for.mean,"'","\n",sep="")



}
if (statistics=="CI"||statistics=="all"){
v.locis=numeric(0)

v.mean=numeric(0)

for (i in 1:pairwise.comparisons){

          # For every pairwise comparison, the following commands are carried out
          # seperately.
          
                                    v.loci=numeric(0)
                                    
                                              # This vector will be filled with the values for the several
                                              # comparisons.
                                    
                                    allelefrequency.pair=as.data.frame(as.matrix(rbind(allelefrequency.pop[[(population.one[i])]],allelefrequency.pop[[(population.two[i])]])))
                                    
                                              # The allelefrequencies of the populations that shall be compared
                                              # pairwise with one another, are selected from the table
                                              # 'allelefrequency.pop' and then combined to a table that is named
                                              # allelefrequency.pair.
                                              # In order to get rid of the population levels that are not included
                                              # in the actual comparison, the transformation to a matrix  and then
                                              # to a data frame is carried out.
                                              
                                    sample.sizes.pair=as.data.frame(as.matrix(rbind(sample.sizes.pop[[(population.one[i])]],sample.sizes.pop[[(population.two[i])]])))
                                                
                                              # The sample.sizes for the actual pair of populations that are
                                              # compared to one another, are combined.
                                              # In order to get rid of the population levels that are not included
                                              # in the actual comparison, the transformation to a matrix  and then
                                              # to a data frame is carried out.
                                              
                                    y.pair=as.data.frame(as.matrix(rbind(y.pop[[(population.one[i])]],y.pop[[(population.two[i])]])))
                                    
                                              # The raw data for the actual two populations are combined to a new
                                              # table.
                                              # In order to get rid of the population levels that are not included
                                              # in the actual comparison, the transformation to a matrix  and then
                                              # to a data frame is carried out.
                                    
                                    allelefrequency.pair2=split(allelefrequency.pair,allelefrequency.pair$locus)
                                    
                                              # The table 'allelefrequency.pair' is splitted to get separated
                                              # the values for the different loci. 
                                              
                                    allelefrequency.pair.pops=split(allelefrequency.pair,as.character(allelefrequency.pair$population))
                                    
                                              # The table 'allelefrequency.pair' is split according to the populations
                                              # that are actually compared.          
                                    
                                    number.loci=length(allelefrequency.pair2)
                                    
                                              # The number of loci that have been examined for the actual
                                              # populationpair
                                    
                                    time1 <- Sys.time()
                                    
                                              # Time before bootstrapping.


cat("\n","====================================================================================================","\n",sep="")
cat("========================= Bootstrapping for Confidence Interval Estimation =========================","\n")
cat("====================================================================================================","\n")    
cat("\n","Pairwise comparison:",names(allelefrequency.pair.pops),"\n",sep=" ")
print(paste("Start of analysis: ",Sys.time(),sep=""))                                    

cat("\n","\n","WARNING: Depending on the size of your input data, the performance of your computer and the number of ","\n","bootstrap resamplings you have chosen, bootstrapping can take very long (hours to days).","\n",sep="")                                    
                                    
                                    Bootstrapping.CI(y.pair,bt,x)
                                    
                                              # Confidence levels for the Gst or D values for the actual pair of
                                              # populations is returned in the object 'confidence.limits'.
                                              # The bootstrap values are available from the objects 'loci'
                                              # and 'means'.
                                              
                                    calc(allelefrequency.pair,sample.sizes.pair,x)
                                    
                                              # The Gst or D values for every locus and the mean Gst or D value over all
                                              # loci is calculated.
                                              # The result is saved in a list called 'values'.          
                                    
                                    for (l in 1:number.loci){
                                    
                                                              allelefrequency.pair3=split(allelefrequency.pair2[[l]],as.character(allelefrequency.pair2[[l]]$population))
                                                    
                                                                        # The table 'allelefrequency.pair2 is split according to the populations
                                                                        # that are actually compared
                                                                        # It is defined again for the case that for one locus, data for one
                                                                        # of the populations would be lacking.
                                                              
                                                              v.loci=rbind(v.loci,cbind(as.numeric(as.vector(values[[1]][l,1])),
                                                              names(allelefrequency.pair2)[l],names(allelefrequency.pair3)[1],
                                                              names(allelefrequency.pair3)[2],confidence.limits[[1]][l,c(1:2)]))
                                                    
                                                                        # In this vector,the Gst or D value for
                                                                        # the actual locus, the actual locus, population one and populations
                                                                        # two that are actually compared with one another and these two
                                                                        # named as a populationpair are combined.
                                                                        # The confidence.levels for the actual locus is added.
                                    
                                                              }
                                    
                                    
                                    loci2 <- split(loci,loci$locus)
                                    
                                              # The bootstrap values are separated in different tables according
                                              # to the loci they belong to.         
                                    
                                    number.loci=length(loci2)
                                    
                                              # The number of loci that have been studied.
                                    
                                    time2 <- Sys.time() 
                                    
                                              # Time after bootstrapping
                                              
                                    bootstrap.time <- difftime(time2,time1,units="secs")   
                                    bootstrap.time.output <- round((as.numeric(bootstrap.time/60)),2)     

cat("\n","... Bootstrapping for confidence interval estimation (comparison: ",names(allelefrequency.pair.pops),") is terminated.","\n",sep="")
print(paste("Duration of the bootstrapping analysis for this comparison (min):",bootstrap.time.output,sep=" "))
cat("\n","\n","Out of ",pairwise.comparisons," pairwise comparisons, ",i," have already been analysed.","\n",sep="")
print(paste("Estimated end of the whole analysis for confidence interval estimation:",((pairwise.comparisons-i)*(as.numeric(bootstrap.time))+Sys.time()),sep=" "))                                     
                                                                      
                                              # The estimated end of the whole analysis  is as many times the bootstrap.time
                                              # as comparisons still have to be carried out.       
                                              
     
                                                            
                                    Result.actual.comparison.locis <- v.loci
                                      Result.actual.comparison.locis <- as.data.frame(Result.actual.comparison.locis)
                                      colnames(Result.actual.comparison.locis) <- c(paste(x,".for.locus",sep=""),"locus","population1","population2","0.95.conf.int.lower","0.95.conf.int.upper")            
                                    
                                    
                                    
                                    Result.actual.comparison.mean <- cbind(values[[2]],names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],confidence.limits[[2]])
                                      Result.actual.comparison.mean <- as.data.frame(Result.actual.comparison.mean)
                                      colnames(Result.actual.comparison.mean) <- c(paste("Mean.",x,sep=""),"population1","population2","0.95.conf.int.lower","0.95.confint.upper")  
                                           
                                           
                                    actual.date <- as.Date(Sys.time())
                                    filename.for.loci <- paste("pairwise.",x,".loci.ci.",actual.date,sep="")
                                    filename.for.loci <- paste(filename.for.loci,".txt",sep="")
                                    filename.for.loci <- paste(filename,".",filename.for.loci,sep="")
                                    
                                    filename.for.mean <- paste("pairwise.",x,".mean.ci.",actual.date,sep="")
                                    filename.for.mean <- paste(filename.for.mean,".txt",sep="")
                                    filename.for.mean <- paste(filename,".",filename.for.mean,sep="")
                                    
                                    write.table(Result.actual.comparison.locis,file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)         
                                    write.table(Result.actual.comparison.mean,file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)         
                                    
                                    interm.result <- list(Result.actual.comparison.locis,Result.actual.comparison.mean)
                                    names(interm.result) <- c(paste(x,".loci.pairwise.comparison",sep=""),paste(x,".mean.pairwise.comparison",sep=""))

cat("\n","++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n",sep="")
cat("++++++++++++++++++++++++++ INTERMEDIATE RESULTS CONFIDENCE INTERVAL LIMITS +++++++++++++++++++++++++","\n")
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n","\n")    
                                    
                                    print(interm.result)
cat("----------------------------------------------------------------------------------------------------","\n")
print(paste("This analysis was finished on",Sys.time()))
cat("The ",x,".mean value was calculated as the arithmetic mean of the ",x,".loci values","\n",sep="")
                                    
cat("\n",x," values for the comparison between ",names(allelefrequency.pair.pops)[[1]]," and ",names(allelefrequency.pair.pops)[[2]],"\n"," for every locus are saved in ","\n","'",filename.for.loci,"'","\n",sep="")
cat("\n",x," values averaged over loci for the comparison between ",names(allelefrequency.pair.pops)[[1]]," and ",names(allelefrequency.pair.pops)[[2]],"\n"," are saved in ","\n","'",filename.for.mean,"'","\n",sep="")
                                    
                                    v.locis=rbind(v.locis,cbind(v.loci))
                                        
                                              # The data frames for
                                              # the pairwise
                                              # comparisons are
                                              # combined.
                                    
                                    v.mean=rbind(v.mean,cbind(values[[2]],names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],confidence.limits[[2]]))
                                    
                                              # In this vector, the Gst or D value over all loci, population one and
                                              # population two that are actually compared with one another and thes
                                              # tw named as populationpair are combined.
                                              # The confidence.levels for the actual locus and the according
                                              # p-values are added.             
                                              
                                  }

                                  
                                  
v.locis=as.data.frame(v.locis)
colnames(v.locis)=c(paste(x,".locus",sep=""),"Locus","Population1","Population2","Lower.0.95.CI","Upper.0.95.CI")            

v.mean=as.data.frame(v.mean)
colnames(v.mean)=c(paste(x,".mean",sep=""),"Population1","Population2","Lower.0.95.CI","Upper.0.95.CI")
       
v.pairwise=list(v.locis,v.mean)
names(v.pairwise)=c(paste(x,".loci.pairwise.comparison",sep=""),paste(x,".mean.pairwise.comparison",sep=""))


rownames(v.pairwise[[1]]) <- seq(1,length(v.pairwise[[1]][,1]))            
rownames(v.pairwise[[2]]) <- seq(1,length(v.pairwise[[2]][,1]))  


cat("\n","++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n",sep="")
cat("++++++++++++++++++++++++++++ END RESULTS WITH CONFIDENCE INTERVAL LIMITS +++++++++++++++++++++++++++","\n")
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n","\n")    
       
print(v.pairwise)
     
assign("v.pairwise",v.pairwise,pos = ".GlobalEnv")

          # The list 'v.pairwise' is assigned to the workspace
          # and therefore available for further calculations.
          
actual.date <- as.Date(Sys.time())
filename.for.loci <- paste("pairwise.",x,".loci.ci.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

filename.for.mean <- paste("pairwise.",x,".mean.ci.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")

                 
write.table(as.data.frame(as.matrix(v.pairwise[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(v.pairwise[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

cat("\n","----------------------------------------------------------------------------------------------------","\n","\n")
print(paste("This analysis was finished on",Sys.time()))
cat("The ",x,".mean values were calculated as the arithmetic means of the ",x,".loci values","\n","for the respective pairwise comparison","\n",sep="")
cat("\n",x, " values for pairwise comparisons between populations for every locus are saved in ","\n","'",filename.for.loci,"'","\n",sep="")
cat("\n",x, " values averaged over all loci for pairwise comparisons between populations are saved in ","\n","'",filename.for.mean,"'","\n",sep="")

     
     }

cat("\n","====================================================================================================","\n",sep="")
cat("====================================================================================================","\n")            
cat("========================================= End of Analysis ==========================================","\n")   
cat("====================================================================================================","\n")
cat("====================================================================================================","\n")    


}
