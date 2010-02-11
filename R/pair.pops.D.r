pair.pops.D <- function (filename,object=FALSE,format.table=TRUE, p.Val=TRUE, bt=1000){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#              filename = object, file (table);

# Passed:
#              filename -> inputformat();
#              tab -> allelefreq();
#              y.pair -> Bootstrapping.Gst();
#              allelefrequency.pair, sample.sizes.pair -> D.calc();
#              empirical.value.loci,bootstrapped.values.loci,empirical.value,bootstrapped.values -> p.value();
#              Gv.pairwise -> p.value.correcture();


# Output:
#              bootstrap.time.output, pairwise.comparisons, bootstrap.time -> screen (print);
#              names(allelefrequency.pair.pops) -> screen (print);
#              interm.result -> screen (print);
#              pairwise.D.loci.[date].txt -> working directory (file);
#              pairwise.D.mean.[date].txt -> working directory (file);
#              Dv.pairwise.adjusted -> screen (print);
#              Dv.pairwise.adjusted -> Workspace;
#              pairwise.D.loci.[date].txt -> working directory (file);
#              pairwise.D.mean.[date].txt -> working directory (file);
#------------------------------------------------------------------------------------------------------------------------------  

          # Function that calculates the D value for each locus separately and
          # the mean D value over all loci, for each possible pairwise population 
          # comparison.

          # bt defines the times that the alleles are resampled by
          # bootstrapping.
  


if (p.Val==TRUE){



if (object==TRUE){

                  tab <- get(filename)
                  
          # The data table can be either an object in the R Workspace (a data table
          # that is already loaded in the Workspace)...
                           
                  }else{
                        tab <- read.table(filename, header=TRUE, sep="")}

          # Or a data table that is saved in a .txt-file and that has to be 
          # assigned to the workspace.

if (format.table==TRUE){
                        inputformat(filename,object) 
                        
                        tab <- read.table("Output-Inputformat.txt", header=TRUE, sep="")
                        
                        }
                       
          # If the argument 'format.table is set as true, the format of the
          # table is changed in that format, that is needed for further 
          # calculations.

          # The table 'tab' has to look like this:

          #     individual population      fragment.length   locus
          # 1        B1.1     Borkum            323          L12
          # 2        B1.1     Borkum            266          L12
          # 3        B1.2     Borkum            325          L12
          # 4        B1.2     Borkum            274          L12
          # 5        B1.3     Borkum            266          L12
          # 6        B1.3     Borkum            323          L12
          # 7        B1.4     Borkum            325          L12

allelefreq(tab)

          # The function 'allelefreq' must be present in the r workspace.
          # In this function, the object 'allelefrequency' and the object
          # 'sample sizes' is calculated.
          
y <- tab

          # The table is saved in the object y and needed to calculate the
          # Bootstrap-value.
          
y.pop=split(y,y$population)

          # The table is splitted according to the several populations
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
          
Dv.locis=numeric(0)

Dv.mean=numeric(0)

for (i in 1:pairwise.comparisons){

          # For every pairwise comparison, the following commands are carried out
          # seperately.
          
                                  Dv.loci=numeric(0)
                                  
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
                                            # populationpair.
                                  
                                  time1 <- Sys.time()
                                  
                                            # Time before bootstrapping.
                                      
                                  cat("\n","====================================================================================================","\n") 
                                  cat("\n","\n","\n","Bootstrapping is carried out for the populations:",names(allelefrequency.pair.pops),"...","\n",sep=" ")
                                  print(Sys.time())
                                  
                                  Bootstrapping.D(y.pair, bt)
                                  
                                            # Confidence levels for the D values for the actual pair of
                                            # populations is returned in the object 'confidence.limits'.
                                            # The bootstrap values are available from the objects 'D.loci'
                                            # and 'D.means'.
                                            
                                  D.calc(allelefrequency.pair,sample.sizes.pair)
                                  
                                            # The D values for every locus and the mean D value over all
                                            # loci is calculated.
                                            # The result is saved in a list called 'D.values'.
                                  
                                  for (l in 1:number.loci){
                                  
                                                            allelefrequency.pair3=split(allelefrequency.pair2[[l]],as.character(allelefrequency.pair2[[l]]$population))
                                                  
                                                                      # The table 'allelefrequency.pair2 is split according to the populations
                                                                      # that are actually compared
                                                                      # It is defined again for the case that for one locus, data for one
                                                                      # of the populations would be lacking.
                                                  
                                                            Dv.loci=rbind(Dv.loci,cbind(as.numeric(as.vector(D.values[[1]]$D[l])),
                                                            names(allelefrequency.pair2)[l],names(allelefrequency.pair3)[1],
                                                            names(allelefrequency.pair3)[2],paste(names(allelefrequency.pair3)[1],
                                                            names(allelefrequency.pair3)[2],sep=""),confidence.limits[[1]][l,c(1:2)]))
                                                  
                                                                      # In this vector,the D value for
                                                                      # the actual locus, the actual locus, population one and populations
                                                                      # two that are actually compared with one another and these two
                                                                      # named as a populationpair are combined.
                                                                      # The confidence.levels for the actual locus is added.
                                  
                                                          }
                                  
                                            # now, the according p-values are calculated and added to the tables.
                                  
                                  D.loci2=split(D.loci,D.loci$locus)
                                  
                                            # The bootstrap values are separated in different tables according
                                            # to the loci they belong to.          
                                  
                                  number.loci=length(D.loci2)
                                  
                                            # The number of loci that have been studied.
                                  
                                  p.values.loci=numeric(0)
                                  
                                            # This vector will be filled with as many p.values as loci have
                                            # been examined.
                                            
                                  time2 <- Sys.time() 
                                  
                                            # Time after bootstrapping
                                  
                                  bootstrap.time <- difftime(time2,time1,units="secs")   
                                  bootstrap.time.output <- round((as.numeric(bootstrap.time/60)),2) 
                                   
                                  cat("\n","----------------------------------------------------------------------------------------------------","\n")  
                                  
                                  cat("\n","Bootstrapping for the comparison between the populations:",names(allelefrequency.pair.pops),"is terminated.","\n","\n",sep=" ")
                                  print(paste("Duration of the bootstrapping analysis for this comparison (min):",bootstrap.time.output,sep=" "))
                                  cat("\n","Out of",pairwise.comparisons,"comparisons",i,"have already been analysed.","\n","\n",sep=" ")
                                  
                                  print(paste("Estimated end of the whole analysis:",((pairwise.comparisons-i)*(as.numeric(bootstrap.time))+Sys.time()),sep=" ")) 
                                                                    
                                            # The estimated end of the whole analysis  is as many times the bootstrap.time
                                            # as comparisons still have to be carried out.                               
                                                
                                  for (l in 1:number.loci){
                                  
                                            # The following commands are carried out separately for the different
                                            # loci.
                                            
                                                              bootstrapped.values.loci=as.numeric(as.vector(D.loci2[[l]]$D))
                                                              empirical.value.loci=as.numeric(as.vector(D.values[[1]]$D[l]))
                                                              
                                                              p.val(empirical.value.loci,bootstrapped.values.loci)
                                                              
                                                                        # This function gives the p-value for the actual empirical value
                                                                        # in the object 'p.value'.
                                                              
                                                              p.values.loci=rbind(p.values.loci,cbind(round(p.value,4),names(D.loci2)[l]))
                                                            
                                                            }
                                  
                                  bootstrapped.values=D.means
                                  empirical.value=D.values[[2]]
                                  
                                  p.value.over.all=p.val(empirical.value,bootstrapped.values)
                                  p.value.over.all=round(p.value,4)
                                  
                                  
                                  Result.actual.comparison.locis <- cbind(Dv.loci,as.numeric(as.vector(p.values.loci[,1])))
                                    Result.actual.comparison.locis <- as.data.frame(Result.actual.comparison.locis)
                                    colnames(Result.actual.comparison.locis) <- c("D.for.locus","locus","population1","population2",
                                         "populationpair","0.95.conf.int.lower","0.95.conf.int.upper","p.values")            
                                  
                                  
                                  
                                  Result.actual.comparison.mean <- cbind(D.values[[2]],names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],
                                            paste(names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],sep=""),
                                            confidence.limits[[2]],p.value.over.all)
                                    Result.actual.comparison.mean <- as.data.frame(Result.actual.comparison.mean)
                                    colnames(Result.actual.comparison.mean) <- c("Mean.D","population1","population2",
                                         "populationpair","0.95.conf.int.lower","0.95.conf.int.upper","p.values")  
                                         
                                         
                                  actual.date <- as.Date(Sys.time())
                                  filename.for.loci <- paste("pairwise.D.loci.",actual.date,sep="")
                                  filename.for.loci <- paste(filename.for.loci,".txt",sep="")
                                  filename.for.loci <- paste(filename,".",filename.for.loci,sep="")
                                  
                                  filename.for.mean <- paste("pairwise.D.mean.",actual.date,sep="")
                                  filename.for.mean <- paste(filename.for.mean,".txt",sep="")
                                  filename.for.mean <- paste(filename,".",filename.for.mean,sep="")
                                  
                                  write.table(Result.actual.comparison.locis,file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)         
                                  write.table(Result.actual.comparison.mean,file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)         
                                  
                                  interm.result <- list(Result.actual.comparison.locis,Result.actual.comparison.mean)
                                  names(interm.result) <- c("D.loci.pairwise.comparison","D.mean.pairwise.comparison")
                                  
                                  cat("\n","----------------------------------------------------------------------------------------------------","\n")  
                                  
                                  cat("\n","\n","++++++++++++++++++++++++++++++++++++++  INTERMEDIATE RESULT  +++++++++++++++++++++++++++++++++++++++","\n","\n") 
                                  
                                  print(interm.result)
                                  
                                  cat("\n","----------------------------------------------------------------------------------------------------","\n")  
                                  
                                  cat("\n","D values for the comparison between",names(allelefrequency.pair.pops)[[1]],
                                      "and",names(allelefrequency.pair.pops)[[2]],"\n","for every locus saved in",
                                      "'",filename.for.loci,"'","\n",sep=" ")
                                  cat("\n","Mean D values for the comparison between",names(allelefrequency.pair.pops)[[1]],
                                      "and",names(allelefrequency.pair.pops)[[2]],"\n","saved in","'",filename.for.mean,
                                      "'","\n",sep=" ")    
                                  
                                  Dv.locis=rbind(Dv.locis,cbind(Dv.loci,as.numeric(as.vector(p.values.loci[,1]))))
                                      
                                            # The according p-values are added to the data.frame Dv.loci and the
                                            # data frames for the pairwise comparisons are combined.
                                  
                                  Dv.mean=rbind(Dv.mean,cbind(D.values[[2]],names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],
                                            paste(names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],sep=""),
                                            confidence.limits[[2]],p.value.over.all))
                                  
                                            # In this vector, the D value over all loci, population one and
                                            # population two that are actually compared with one another and thes
                                            # tw named as populationpair are combined.
                                            # The confidence.levels for the actual locus and the according
                                            # p-values are added.
               
                                }

Dv.locis=as.data.frame(Dv.locis)
colnames(Dv.locis)=c("D.for.locus","locus","population1","population2",
"populationpair","0.95.conf.int.lower","0.95.conf.int.upper","p.values")            
       
Dv.mean=as.data.frame(Dv.mean)
colnames(Dv.mean)=c("Mean.D","population1","population2",
"populationpair","0.95.conf.int.lower","0.95.conf.int.upper","p.values")
       
Dv.pairwise=list(Dv.locis,Dv.mean)
names(Dv.pairwise)=c("D.loci.pairwise.comparison","D.mean.pairwise.comparison")
       
p.value.correcture(Dv.pairwise)

          # A function that adjusts the p-values after Bonferroni, Holm, Hommel
          # and Benjamini and Hochberg.
          # The result is a list called 'Dv.pairwise.adjusted'.

D.pairwise.adjusted <- Dv.pairwise.adjusted


names(D.pairwise.adjusted)=c("D.loci.pairwise.comparison","D.mean.pairwise.comparison")

          
colnames(D.pairwise.adjusted[[1]]) <- c("D.for.locus","locus","population1","population2",
       "populationpair","0.95.conf.int.lower","0.95.conf.int.upper","p.values","p.bonferroni","p.holm","p.hommel","pBH")  
       
colnames(D.pairwise.adjusted[[2]]) <- c("Mean.D","population1","population2",
       "populationpair","0.95.conf.int.lower","0.95.conf.int.upper","p.values","p.bonferroni","p.holm","p.hommel","pBH")   
       
rownames(D.pairwise.adjusted[[1]]) <- seq(1,length(D.pairwise.adjusted[[1]][,1]))            
rownames(D.pairwise.adjusted[[2]]) <- seq(1,length(D.pairwise.adjusted[[2]][,1]))  
cat("\n","----------------------------------------------------------------------------------------------------","\n")         

cat("\n","\n","+++++++++++++++++++++++++++++++++++++++++++  END RESULT  +++++++++++++++++++++++++++++++++++++++++++","\n","\n")  

        
print(D.pairwise.adjusted)
       
assign("D.pairwise.adjusted",D.pairwise.adjusted,pos = ".GlobalEnv")

          # The list 'D.pairwise.adjusted' is assigned to the workspace
          # and therefore available for further calculations.
          
actual.date <- as.Date(Sys.time())
filename.for.loci <- paste("pairwise.D.loci.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

filename.for.mean <- paste("pairwise.D.mean.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")
 
write.table(as.data.frame(as.matrix(D.pairwise.adjusted[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(D.pairwise.adjusted[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
cat("\n","----------------------------------------------------------------------------------------------------","\n","\n")         
 
print(paste("This analysis was carried out on",Sys.time()))
cat("\n","D values for pairwise comparisons between populations for every locus are saved in ","'",filename.for.loci,"'","\n",sep="")
cat("\n","Mean D values over all loci for pairwise comparisons between populations saved in ","'",filename.for.mean,"'","\n","\n","\n",sep="")
cat("\n","====================================================================================================","\n")  


} else {

  # The following calculations are made when the p.values are not needed.
if (object==TRUE){

                  tab <- get(filename)
                  
          # The data table can be either an object in the R Workspace (a data table
          # that is already loaded in the Workspace)...
                           
                  }else{
                        tab <- read.table(filename, header=TRUE, sep="")}

          # Or a data table that is saved in a .txt-file and that has to be 
          # assigned to the workspace.

if (format.table==TRUE){
                        inputformat(filename,object) 
                        
                        tab <- read.table("Output-Inputformat.txt", header=TRUE, sep="")
                        
                        }
                       
          # If the argument 'format.table is set as true, the format of the
          # table is changed in that format, that is needed for further 
          # calculations.

          # The table 'tab' has to look like this:

          #     individual population      fragment.length   locus
          # 1        B1.1     Borkum            323          L12
          # 2        B1.1     Borkum            266          L12
          # 3        B1.2     Borkum            325          L12
          # 4        B1.2     Borkum            274          L12
          # 5        B1.3     Borkum            266          L12
          # 6        B1.3     Borkum            323          L12
          # 7        B1.4     Borkum            325          L12

allelefreq(tab)

          # The function 'allelefreq' must be present in the r workspace.
          # In this function, the object 'allelefrequency' and the object
          # 'sample sizes' is calculated.
          
y <- tab

          # The table is saved in the object y and needed to calculate the
          # Bootstrap-value.
          
y.pop=split(y,y$population)

          # The table is splitted according to the several populations
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
          
Dv.locis=numeric(0)

Dv.mean=numeric(0)

for (i in 1:pairwise.comparisons){

          # For every pairwise comparison, the following commands are carried out
          # seperately.
          
                                  Dv.loci=numeric(0)
                                  
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
                                            # populationpair.
                                  
                                            
                                  D.calc(allelefrequency.pair,sample.sizes.pair)
                                  
                                            # The D values for every locus and the mean D value over all
                                            # loci is calculated.
                                            # The result is saved in a list called 'D.values'.
                                  
                                  for (l in 1:number.loci){
                                  
                                                            allelefrequency.pair3=split(allelefrequency.pair2[[l]],as.character(allelefrequency.pair2[[l]]$population))
                                                  

                                        # The table
                                        # 'allelefrequency.pair2 is
                                        # split according to the
                                        # populations that are
                                        # actually compared. It is
                                        # defined again for the case
                                        # that for one locus, data for
                                        # one of the populations would
                                        # be lacking.
                                                  
                                                            Dv.loci=rbind(Dv.loci,cbind(as.numeric(as.vector(D.values[[1]]$D[l])),
                                                            names(allelefrequency.pair2)[l],names(allelefrequency.pair3)[1],
                                                            names(allelefrequency.pair3)[2],paste(names(allelefrequency.pair3)[1],
                                                            names(allelefrequency.pair3)[2],sep="")))
                                                  
                                         # In this vector,the D value
                                         # for the actual locus, the
                                         # actual locus, population
                                         # one and population two that
                                         # are actually compared with
                                         # one another and these two
                                         # named as a populationpair
                                         # are combined.
                                  
                                                          }
                                                

                                  Result.actual.comparison.locis <- Dv.loci
                                    Result.actual.comparison.locis <- as.data.frame(Result.actual.comparison.locis)
                                    colnames(Result.actual.comparison.locis) <- c("D.for.locus","locus","population1","population2","populationpair")            
                                  
                                  
                                  
                                  Result.actual.comparison.mean <- cbind(D.values[[2]],names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],paste(names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],sep=""))
                                    Result.actual.comparison.mean <- as.data.frame(Result.actual.comparison.mean)
                                    colnames(Result.actual.comparison.mean) <- c("Mean.D","population1","population2",
                                         "populationpair")
                                         
                                         
                                  
                                  interm.result <- list(Result.actual.comparison.locis,Result.actual.comparison.mean)
                                  names(interm.result) <- c("D.loci.pairwise.comparison","D.mean.pairwise.comparison")
                                  
                                  
                                  
                                  Dv.locis=rbind(Dv.locis,Dv.loci)
                                      
                                            # The data frames for the
                                            # pairwise comparisons are
                                            # combined.
                                  
                                  Dv.mean=rbind(Dv.mean,cbind(D.values[[2]],names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],paste(names(allelefrequency.pair.pops)[1],names(allelefrequency.pair.pops)[2],sep="")))
                                  
                                            # In this vector, the D value over all loci, population one and
                                            # population two that are actually compared with one another and these
                                            # two named as populationpair are combined.
                                                          
                                }

Dv.locis=as.data.frame(Dv.locis)
colnames(Dv.locis)=c("D.for.locus","locus","population1","population2",
"populationpair")            
       
Dv.mean=as.data.frame(Dv.mean)
colnames(Dv.mean)=c("Mean.D","population1","population2",
"populationpair")
       
Dv.pairwise=list(Dv.locis,Dv.mean)
names(Dv.pairwise)=c("D.loci.pairwise.comparison","D.mean.pairwise.comparison")
       

D.pairwise.adjusted <- Dv.pairwise


names(D.pairwise.adjusted)=c("D.loci.pairwise.comparison","D.mean.pairwise.comparison")

          
colnames(D.pairwise.adjusted[[1]]) <- c("D.for.locus","locus","population1","population2",
       "populationpair")  
       
colnames(D.pairwise.adjusted[[2]]) <- c("Mean.D","population1","population2",
       "populationpair")   
       
rownames(D.pairwise.adjusted[[1]]) <- seq(1,length(D.pairwise.adjusted[[1]][,1]))            
rownames(D.pairwise.adjusted[[2]]) <- seq(1,length(D.pairwise.adjusted[[2]][,1]))  
cat("\n","----------------------------------------------------------------------------------------------------","\n")         

cat("\n","\n","+++++++++++++++++++++++++++++++++++++++++++  END RESULT  +++++++++++++++++++++++++++++++++++++++++++","\n","\n")  

        
print(D.pairwise.adjusted)
       
assign("D.pairwise.adjusted",D.pairwise.adjusted,pos = ".GlobalEnv")

          # The list 'D.pairwise.adjusted' is assigned to the workspace
          # and therefore available for further calculations.
          
actual.date <- as.Date(Sys.time())
filename.for.loci <- paste("pairwise.D.loci.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

filename.for.mean <- paste("pairwise.D.mean.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")
 
write.table(as.data.frame(as.matrix(D.pairwise.adjusted[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(D.pairwise.adjusted[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
cat("\n","----------------------------------------------------------------------------------------------------","\n","\n")         
 
print(paste("This analysis was carried out on",Sys.time()))
cat("\n","D values for pairwise comparisons between populations for every locus are saved in ","'",filename.for.loci,"'","\n",sep="")
cat("\n","Mean D values over all loci for pairwise comparisons between populations saved in ","'",filename.for.mean,"'","\n","\n","\n",sep="")
cat("\n","====================================================================================================","\n")  





                                
}}
