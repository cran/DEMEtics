all.pops.Dest <- function (filename,object=FALSE,format.table=TRUE, p.Val=TRUE, bt=1000){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           filename = object, file (table);

# Received:
#           tab <- Output-Inputformat.txt (Working directory);
#           allelefrequency,sample.sizes <- allelefreq();

# Passed:
#           filename -> inputformat();
#           tab -> Bootstrapping(), allelefreq();
#           allelefrequency, sample.sizes -> Dest.calc();
#           empirical.value.loci,bootstrapped.values.loci,empirical.value,bootstrapped.values -> p.value();

# Output:
#           Dest.all.pops -> Workspace;
#           names(tab.pops), bootstrap.time.output, Dest.all.pops -> screen (print);
#           all.pops.Dest.loci.[date].txt -> working directory (file);
#           all.pops.Dest.mean.[date].txt -> working directory (file);
#           Dest.sample.sizes.[date].txt -> working directory (file);
#           Dest.allelefrequencies.[date].txt -> working directory (file);
#------------------------------------------------------------------------------------------------------------------------------  

if (p.Val==TRUE){
  
  # bt defines the times that the alleles are resampled by bootstrapping.

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

          # This function will calculate the Dest values for every locus and
          # the mean Dest value over all loci for all the populations that 
          # have been sampled.
          
          # The table 'tab' has to be of the following format:
          
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
          
tab.pops <- split(tab,tab$population)

          # The table is splitted up into the data that belong to different
          # populations.           

cat("\n","====================================================================================================","\n")          
cat("\n","\n","\n","Bootstrapping is carried out for the populations:",names(tab.pops),"...","\n",sep=" ")

          # User information.
          
print(Sys.time())

          # The actual time is printed.

time1 <- Sys.time()

          # The time before the bootstrapping.
          
Bootstrapping.Dest(tab, bt)

          # The confidence limits of the Dest values for the several loci and 
          # over all loci are determined by a thousandfold resampling of the
          # alleles (for each locus and all populations) and the calculation of 
          # the Dest value from those resampled data tables.
          # The results are available from the object 'confidence.limits'
          # The bootstrap values are available from the objects 'Dest.loci'
          # and 'Dest.means'.    
          
time2 <- Sys.time()

          # The time after the bootstrapping.

bootstrap.time <- difftime(time2,time1,units="secs") 

          # Duration of the Bootstrapping analysis (in seconds).
            
bootstrap.time.output <- round((as.numeric(bootstrap.time/60)),2) 

          # Duration of the Bootstrapping analysis (in minutes), rounded (2 digits).
 
cat("\n","----------------------------------------------------------------------------------------------------","\n")                   
cat("\n","Bootstrapping for the comparison between the populations:",names(tab.pops),"is terminated.","\n",sep=" ")
print(paste("Duration of the bootstrapping analysis for this comparison (min):",bootstrap.time.output,sep=" "))

          # User information
          
allelefreq(tab)

          # From the table that contains the empirical data, the allelefrequencies
          # as well as the sample.sizes are calculated for each locus.
          # The table that contains the allelefrequencies and the table that
          # lists the sample sizes are assigned to the R workspace in the 
          # object List, but also separately in the object allelefrequency
          # and the object sample.sizes.         
          

Dest.calc(allelefrequency,sample.sizes)

          # The Dest values for the several loci and the mean Dest value of 
          # all loci are calculated.
          # The results are available from the object 'D.values'. 
                 
          # From the bootstrapping data, the p-values are assigned to the
          # empirical Dest values.
          
Dest.loci2 <- split(Dest.loci,Dest.loci$locus)

          # The bootstrap values are separated in different tables according
          # to the loci they belong to.
          
number.loci <- length(Dest.loci2)

          # The number of loci that have been studied.

p.values.loci <- numeric(0)

          # This vector will be filled with as many p-values as loci have
          # been examined.
         
for (l in 1:number.loci){

          # The following commands are carried out separately for the different
          # loci.
          
                          bootstrapped.values.loci=as.numeric(as.vector(Dest.loci2[[l]]$Dest))

                                    # The bootstrapping values (Dest values) that have been obtained for 
                                    # the several loci are assigned to a vector.
          
                          empirical.value.loci=as.numeric(as.vector(D.values[[1]]$Dest[l]))
                          
                                    # The empirical Dest values for the several loci are assigned to a
                                    # vector.

                          p.val(empirical.value.loci,bootstrapped.values.loci)

                                    # This function returns the p-value for the actual empirical Dest value
                                    # in the object 'p.value'.                          
          
                          p.values.loci=rbind(p.values.loci,cbind(round(p.value,4),names(Dest.loci2)[l]))
                          
                                    # The p-values (rounded up to 4 decimal places) for the several loci
                                    # are combined.
                                    # The p-value can't be more exact, when the bootstrapping is based on
                                    # thousand resamplings.                          
                          
                          }
                          
bootstrapped.values=Dest.means

          # The mean Dest values over all loci that have been obtained from
          # the bootstrapping.

empirical.value=D.values[[2]]

          # The empirical mean Dest value over all loci. 

p.value.over.all=p.val(empirical.value,bootstrapped.values)

          # This function returns the p-value for the actual empirical Dest value
          # in the object 'p.value'.
          
p.value.over.all=round(p.value,4)

          # This p.value is rounded up to 4 decimal places.

D.loci=cbind(D.values[[1]],confidence.limits[[1]][,c(1:2)],as.numeric(as.vector(p.values.loci[,1])))
D.loci=as.data.frame(D.loci)
colnames(D.loci)=c("Dest","locus","0.95.conf.int.lower","0.95.conf.int.upper","p.value")

          # A table is created that contains the Dest values for the several
          # loci with the according confidence levels and the p-value obtained
          # from the bootstrapping analysis.

D.mean=cbind(D.values[[2]],confidence.limits[[2]],p.value.over.all)
D.mean=as.data.frame(D.mean)
colnames(D.mean)=c("Dest.mean","0.95.conf.int.lower","0.95.conf.int.upper","p.value")

          # A table is created that contains the mean Dest value over all
          # loci with the according confidence levels and the p-value obtained
          # from the bootstrapping analysis.

Dest.all.pops=list(D.loci,D.mean)
names(Dest.all.pops)=c("Dest.loci.over.all.populations","Dest.mean.over.all.populations")

          # These two tables are combined in a single list.

assign("Dest.all.pops",Dest.all.pops,pos = ".GlobalEnv")

          # The object 'Dest.all.pops', a list, is assigned to the R workspace.

cat("\n","----------------------------------------------------------------------------------------------------","\n")         
cat("\n","\n","+++++++++++++++++++++++++++++++++++++++++++  END RESULT  +++++++++++++++++++++++++++++++++++++++++++","\n","\n")  
print(Dest.all.pops)

          # The result tables are printed.

actual.date <- as.Date(Sys.time())

          # The actual date is determined.

filename.for.loci <- paste("all.pops.Dest.loci.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

          # The filename under which the table that contains the Dest values
          # for all the loci seperately, is created. It will encompass the actual
          # date at the end.

filename.for.mean <- paste("all.pops.Dest.mean.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")

          # The filename under which the table that contains the mean Dest value
          # over all the loci, is created. It will encompass the actual
          # date at the end.

filename.for.sample.sizes <- paste("Dest.sample.sizes.",actual.date,sep="")
filename.for.sample.sizes <- paste(filename.for.sample.sizes,".txt",sep="")
filename.for.sample.sizes <- paste(filename,".",filename.for.sample.sizes,sep="") 

          # The filename under which the table that contains the sample sizes
          # for the several loci, is created. It will encompass the actual
          # date at the end.

filename.for.allelefrequencies <- paste("Dest.allelefrequencies.",actual.date,sep="")
filename.for.allelefrequencies <- paste(filename.for.allelefrequencies,".txt",sep="")
filename.for.allelefrequencies <- paste(filename,".",filename.for.allelefrequencies,sep="")

          # The filename under which the table that contains the allelefrequencies
          # for the several loci, is created. It will encompass the actual
          # date at the end.

write.table(as.data.frame(as.matrix(Dest.all.pops[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(Dest.all.pops[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(List[[1]])),file=filename.for.allelefrequencies, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(List[[2]])),file=filename.for.sample.sizes, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

          # The tables are saved in the working directory. Two or more analysis
          # that are carried under the same working directory, are saved in the
          # under the same file name. The data are combined in the order as they
          # have been analysed. No data is lost due to be overwritten.

cat("\n","----------------------------------------------------------------------------------------------------","\n","\n")         
print(paste("This analysis was carried out on",Sys.time()))

cat("\n","\n","Dest values for pairwise comparisons between populations for every locus are saved in ","'",filename.for.loci,"'","\n",sep="")
cat("\n","Mean Dest values over all loci for pairwise comparisons between populations saved in ","'",filename.for.mean,"'","\n",sep="")
cat("\n","Table with allelefrequencies is saved in ","'",filename.for.allelefrequencies,"'","\n",sep="")
cat("\n","Per locus calculated sample sizes are saved in ","'",filename.for.sample.sizes,"'","\n","\n","\n",sep="")
cat("\n","====================================================================================================","\n")   

          # User information about the end date of the analysis and the filenames
          # under which the several tables have been saved.

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

          # This function will calculate the Dest values for every locus and
          # the mean Dest value over all loci for all the populations that 
          # have been sampled.
          
          # The table 'tab' has to be of the following format:
          
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

Dest.calc(allelefrequency,sample.sizes)

          # The Dest values for the several loci and the mean Dest value of 
          # all loci are calculated.
          # The results are available from the object 'D.values'. 
                 
          # From the bootstrapping data, the p-values are assigned to the
          # empirical Dest values.
          
Dest.loci2 <- split(tab,tab$locus)

          # The bootstrap values are separated in different tables according
          # to the loci they belong to.
          
number.loci <- length(Dest.loci2)

          # The number of loci that have been studied.

         
for (l in 1:number.loci){

          # The following commands are carried out separately for the different
          # loci.
         
                          empirical.value.loci=as.numeric(as.vector(D.values[[1]]$Dest[l]))
                          
                                    # The empirical Dest values for the several loci are assigned to a
                                    # vector.                          
                          
                          }
                          


empirical.value=D.values[[2]]

          # The empirical mean Dest value over all loci. 


D.loci=D.values[[1]]
D.loci=as.data.frame(D.loci)
colnames(D.loci)=c("Dest","locus")

          # A table is created that contains the Dest values for the several
          # loci 

D.mean=D.values[[2]]
D.mean=as.data.frame(D.mean)
colnames(D.mean)=c("Dest.mean")

          # A table is created that contains the mean Dest value over all
          # loci with the according confidence levels and the p-value obtained
          # from the bootstrapping analysis.

Dest.all.pops=list(D.loci,D.mean)
names(Dest.all.pops)=c("Dest.loci.over.all.populations","Dest.mean.over.all.populations")

          # These two tables are combined in a single list.

assign("Dest.all.pops",Dest.all.pops,pos = ".GlobalEnv")

          # The object 'Dest.all.pops', a list, is assigned to the R workspace.

cat("\n","----------------------------------------------------------------------------------------------------","\n")         
cat("\n","\n","+++++++++++++++++++++++++++++++++++++++++++  END RESULT  +++++++++++++++++++++++++++++++++++++++++++","\n","\n")  
print(Dest.all.pops)

          # The result tables are printed.

actual.date <- as.Date(Sys.time())

          # The actual date is determined.

filename.for.loci <- paste("all.pops.Dest.loci.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

          # The filename under which the table that contains the Dest values
          # for all the loci seperately, is created. It will encompass the actual
          # date at the end.

filename.for.mean <- paste("all.pops.Dest.mean.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")

          # The filename under which the table that contains the mean Dest value
          # over all the loci, is created. It will encompass the actual
          # date at the end.

filename.for.sample.sizes <- paste("Dest.sample.sizes.",actual.date,sep="")
filename.for.sample.sizes <- paste(filename.for.sample.sizes,".txt",sep="")
filename.for.sample.sizes <- paste(filename,".",filename.for.sample.sizes,sep="") 

          # The filename under which the table that contains the sample sizes
          # for the several loci, is created. It will encompass the actual
          # date at the end.

filename.for.allelefrequencies <- paste("Dest.allelefrequencies.",actual.date,sep="")
filename.for.allelefrequencies <- paste(filename.for.allelefrequencies,".txt",sep="")
filename.for.allelefrequencies <- paste(filename,".",filename.for.allelefrequencies,sep="")

          # The filename under which the table that contains the allelefrequencies
          # for the several loci, is created. It will encompass the actual
          # date at the end.

write.table(as.data.frame(as.matrix(Dest.all.pops[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(Dest.all.pops[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(List[[1]])),file=filename.for.allelefrequencies, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(List[[2]])),file=filename.for.sample.sizes, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

          # The tables are saved in the working directory. Two or more analysis
          # that are carried under the same working directory, are saved in the
          # under the same file name. The data are combined in the order as they
          # have been analysed. No data is lost due to be overwritten.

cat("\n","----------------------------------------------------------------------------------------------------","\n","\n")         
print(paste("This analysis was carried out on",Sys.time()))

cat("\n","\n","Dest values for pairwise comparisons between populations for every locus are saved in ","'",filename.for.loci,"'","\n",sep="")
cat("\n","Mean Dest values over all loci for pairwise comparisons between populations saved in ","'",filename.for.mean,"'","\n",sep="")
cat("\n","Table with allelefrequencies is saved in ","'",filename.for.allelefrequencies,"'","\n",sep="")
cat("\n","Per locus calculated sample sizes are saved in ","'",filename.for.sample.sizes,"'","\n","\n","\n",sep="")
cat("\n","====================================================================================================","\n")   

          # User information about the end date of the analysis and the filenames
          # under which the several tables have been saved.



}
}
  
    

                 




