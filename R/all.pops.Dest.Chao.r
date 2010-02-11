all.pops.Dest.Chao <- function (filename,object=FALSE,format.table=TRUE, p.Val=TRUE, bt=1000){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#              filename = object, file (table);

# Passed:
#              tab -> Bootstrapping.Chao(), Dest.Chao();
#              filename -> inputformat()
#              empirical.value.loci,bootstrapped.values.loci,empirical.value,bootstrapped.values -> p.value();

# Output:
#              Dest.Chao.all.pops -> Workspace;
#              bootstrap.time.output -> screen (print);
#              names(tab.pops) -> screen (print);
#              Dest.Chao.all.pops -> screen (print);
#              all.pops.Dest.Chao.loci.[date].txt -> working directory (file);
#              all.pops.Dest.Chao.mean.[date].txt -> working directory (file);
#              Dest.Chao.sample.sizes.[date].txt -> working directory (file);
#              Dest.Chao.allelefrequencies.[date].txt -> working directory (file);
#------------------------------------------------------------------------------------------------------------------------------  

  # bt defines the times that the alleles are resampled by bootstrapping.  

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

          # This function will calculate the Dest.Chao values for every locus and
          # the mean Dest.Chao value over all loci for all the populations that 
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
   
Bootstrapping.Chao(tab, bt)

          # The confidence.limits of the Dest.Chao values for the several loci and 
          # over all loci are determined by a thousandfold resampling of the
          # alleles (for each locus and all populations) and the calculation of 
          # the Dest.Chao value from those resampled data tables.
          # The results are available from the object 'confidence.limits'
          # The bootstrap values are available from the objects 'Dest.Chao.loci'
          # and 'Dest.Chao.means'. 
          
time2 <- Sys.time() 

          # The time after the bootstrapping.

bootstrap.time <- difftime(time2,time1,units="secs") 

          # Duration of the Bootstrapping analysis (in seconds).
  
bootstrap.time.output <- round((as.numeric(bootstrap.time/60)),2) 

          # Duration of the Bootstrapping analysis (in minutes), rounded (2 digits).
  
cat("\n","----------------------------------------------------------------------------------------------------","\n")                   
cat("\n","\n","\n","Bootstrapping for the comparison between the populations:",names(tab.pops),"is terminated.","\n",sep=" ")
print(paste("Duration of the bootstrapping analysis for this comparison (min):",bootstrap.time.output,sep=" "))

          # User information
                             
Dest.Chao(tab)        
                    
          # The Dest.Chao values for the several loci and the mean Dest.Chao value of 
          # all loci are calculated.
          # The results are available from the object 'D.Chao.values'. 
          
          # The table that contains the allelefrequencies and the table that
          # lists the sample sizes are present in the R workspace in the 
          # object List, but also separately in the object allelefrequency
          # and the object sample.sizes 
                   
          # From the bootstrapping data, the p-values are assigned to the
          # empirical Dest.Chao values.
          
Dest.Chao.loci2=split(Dest.Chao.loci,Dest.Chao.loci$locus)

          # The bootstrap values are separated in different tables according
          # to the loci they belong to.
          
number.loci=length(Dest.Chao.loci2)

          # The number of loci that have been studied.

p.values.loci=numeric(0)

          # This vector will be filled with as many p.values as loci have
          # been examined.

for (l in 1:number.loci){

          # The following commands are carried out separately for the different
          # loci.
          
                          bootstrapped.values.loci=as.numeric(as.vector(Dest.Chao.loci2[[l]]$Dest.Chao))
                          
                                    # The bootstrapping values (Dest.Chao values) that have been obtained for 
                                    # the several loci are assigned to a vector.
                          
                          empirical.value.loci=as.numeric(as.vector(D.Chao.values[[1]]$Dest.Chao[l]))
                          
                                    # The empirical Dest.Chao values for the several loci are assigned to a
                                    # vector.                          

                          p.val(empirical.value.loci,bootstrapped.values.loci)

                                    # This function returns the p-value for the actual empirical Dest.Chao value
                                    # in the object 'p.value'.
          
                          p.values.loci=rbind(p.values.loci,cbind(round(p.value,4),names(Dest.Chao.loci2)[l]))
                          
                                    # The p-values (rounded up to 4 decimal places) for the several loci
                                    # are combined.
                                    # The p-value can't be more exact, when the bootstrapping is based on
                                    # thousand resamplings.                           
                          }

bootstrapped.values=Dest.Chao.means

          # The mean Dest.Chao values over all loci that have been obtained from
          # the bootstrapping.

empirical.value=D.Chao.values[[2]]

          # The empirical mean Dest.Chao value over all loci. 

p.value.over.all=p.val(empirical.value,bootstrapped.values)

          # This function returns the p-value for the actual empirical Dest.Chao value
          # in the object 'p.value'.

p.value.over.all=round(p.value,4)

          # This p.value is rounded up to 4 decimal places.

D.loci=cbind(D.Chao.values[[1]],confidence.limits[[1]][,c(1:2)],as.numeric(as.vector(p.values.loci[,1])))
D.loci=as.data.frame(D.loci)
colnames(D.loci)=c("Dest.Chao","locus","0.95.conf.int.lower","0.95.conf.int.upper","p.value")

          # A table is created that contains the Dest.Chao values for the several
          # loci with the according confidence levels and the p-value obtained
          # from the bootstrapping analysis.

D.mean=cbind(D.Chao.values[[2]],confidence.limits[[2]],p.value.over.all)
D.mean=as.data.frame(D.mean)
colnames(D.mean)=c("Dest.Chao.mean","0.90.conf.int.lower","0.95.conf.int.upper","p.value")

          # A table is created that contains the mean Dest value over all
          # loci with the according confidence levels and the p-value obtained
          # from the bootstrapping analysis.

Dest.Chao.all.pops=list(D.loci,D.mean)
names(Dest.Chao.all.pops)=c("Dest.Chao.loci.over.all.populations","Dest.Chao.mean.over.all.populations")

          # These two tables are combined in a single list.

assign("Dest.Chao.all.pops",Dest.Chao.all.pops,pos = ".GlobalEnv")

          # The object 'Dest.Chao.all.pops', a list, is assigned to the R workspace.

cat("\n","----------------------------------------------------------------------------------------------------","\n")         
cat("\n","\n","+++++++++++++++++++++++++++++++++++++++++++  END RESULT  +++++++++++++++++++++++++++++++++++++++++++","\n","\n")  
print(Dest.Chao.all.pops)

          # The result tables are printed.
          
actual.date <- as.Date(Sys.time())

          # The actual date is determined.

filename.for.loci <- paste("all.pops.Dest.Chao.loci.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

          # The filename under which the table that contains the Dest.Chao values
          # for all the loci seperately, is created. It will encompass the actual
          # date at the end.

filename.for.mean <- paste("all.pops.Dest.Chao.mean.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")

          # The filename under which the table that contains the mean Dest.Chao value
          # over all the loci, is created. It will encompass the actual
          # date at the end.

filename.for.sample.sizes <- paste("Dest.Chao.sample.sizes.",actual.date,sep="")
filename.for.sample.sizes <- paste(filename.for.sample.sizes,".txt",sep="")
filename.for.sample.sizes <- paste(filename,".",filename.for.sample.sizes,sep="") 

          # The filename under which the table that contains the sample sizes
          # for the several loci, is created. It will encompass the actual
          # date at the end.

filename.for.allelefrequencies <- paste("Dest.Chao.allelefrequencies.",actual.date,sep="")
filename.for.allelefrequencies <- paste(filename.for.allelefrequencies,".txt",sep="")
filename.for.allelefrequencies <- paste(filename,".",filename.for.allelefrequencies,sep="")

          # The filename under which the table that contains the allelefrequencies
          # for the several loci, is created. It will encompass the actual
          # date at the end.

write.table(as.data.frame(as.matrix(Dest.Chao.all.pops[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(Dest.Chao.all.pops[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(List[[1]])),file=filename.for.allelefrequencies, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(List[[2]])),file=filename.for.sample.sizes, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

          # The tables are saved in the working directory. Two or more analysis
          # that are carried under the same working directory, are saved in the
          # under the same file name. The data are combined in the order as they
          # have been analysed. No data is lost due to be overwritten.
          
cat("\n","----------------------------------------------------------------------------------------------------","\n","\n")         
print(paste("This analysis was carried out on",Sys.time()))

cat("\n","\n","Dest.Chao values for pairwise comparisons between populations for every locus are saved in ","'",filename.for.loci,"'","\n",sep="")
cat("\n","Mean Dest.Chao values over all loci for pairwise comparisons between populations saved in ","'",filename.for.mean,"'","\n",sep="")
cat("\n","Table with allelefrequencies is saved in ","'",filename.for.allelefrequencies,"'","\n",sep="")
cat("\n","Per locus calculated sample sizes are saved in ","'",filename.for.sample.sizes,"'","\n","\n","\n",sep="")  
cat("\n","====================================================================================================","\n")        

          # User information about the end date of the analysis and the filenames
          # under which the several tables have been saved.

}else{

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

          # This function will calculate the Dest.Chao values for every locus and
          # the mean Dest.Chao value over all loci for all the populations that 
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
         

Dest.Chao(tab)        
                    
          # The Dest.Chao values for the several loci and the mean Dest.Chao value of 
          # all loci are calculated.
          # The results are available from the object 'D.Chao.values'. 
          
          # The table that contains the allelefrequencies and the table that
          # lists the sample sizes are present in the R workspace in the 
          # object List, but also separately in the object allelefrequency
          # and the object sample.sizes 
                   
      
          
Dest.Chao.loci2=split(tab,tab$locus)

          # The bootstrap values are separated in different tables according
          # to the loci they belong to.
          
number.loci=length(Dest.Chao.loci2)

          # The number of loci that have been studied.

for (l in 1:number.loci){

          # The following commands are carried out separately for the different
          # loci.
          
                           empirical.value.loci=as.numeric(as.vector(D.Chao.values[[1]]$Dest.Chao[l]))
                          
                                    # The empirical Dest.Chao values for the several loci are assigned to a
                                    # vector.                          

                          }

empirical.value=D.Chao.values[[2]]

          # The empirical mean Dest.Chao value over all loci. 


D.loci=D.Chao.values[[1]]
D.loci=as.data.frame(D.loci)
colnames(D.loci)=c("Dest.Chao","locus")

          # A table is created that contains the Dest.Chao values for the several
          # loci with the according confidence levels and the p-value obtained
          # from the bootstrapping analysis.

D.mean=D.Chao.values[[2]]
D.mean=as.data.frame(D.mean)
colnames(D.mean)=c("Dest.Chao.mean")

          # A table is created that contains the mean Dest value over all
          # loci with the according confidence levels and the p-value obtained
          # from the bootstrapping analysis.

Dest.Chao.all.pops=list(D.loci,D.mean)
names(Dest.Chao.all.pops)=c("Dest.Chao.loci.over.all.populations","Dest.Chao.mean.over.all.populations")

          # These two tables are combined in a single list.

assign("Dest.Chao.all.pops",Dest.Chao.all.pops,pos = ".GlobalEnv")

          # The object 'Dest.Chao.all.pops', a list, is assigned to the R workspace.

cat("\n","----------------------------------------------------------------------------------------------------","\n")         
cat("\n","\n","+++++++++++++++++++++++++++++++++++++++++++  END RESULT  +++++++++++++++++++++++++++++++++++++++++++","\n","\n")  
print(Dest.Chao.all.pops)

          # The result tables are printed.
          
actual.date <- as.Date(Sys.time())

          # The actual date is determined.

filename.for.loci <- paste("all.pops.Dest.Chao.loci.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

          # The filename under which the table that contains the Dest.Chao values
          # for all the loci seperately, is created. It will encompass the actual
          # date at the end.

filename.for.mean <- paste("all.pops.Dest.Chao.mean.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")

          # The filename under which the table that contains the mean Dest.Chao value
          # over all the loci, is created. It will encompass the actual
          # date at the end.

filename.for.sample.sizes <- paste("Dest.Chao.sample.sizes.",actual.date,sep="")
filename.for.sample.sizes <- paste(filename.for.sample.sizes,".txt",sep="")
filename.for.sample.sizes <- paste(filename,".",filename.for.sample.sizes,sep="") 

          # The filename under which the table that contains the sample sizes
          # for the several loci, is created. It will encompass the actual
          # date at the end.

filename.for.allelefrequencies <- paste("Dest.Chao.allelefrequencies.",actual.date,sep="")
filename.for.allelefrequencies <- paste(filename.for.allelefrequencies,".txt",sep="")
filename.for.allelefrequencies <- paste(filename,".",filename.for.allelefrequencies,sep="")

          # The filename under which the table that contains the allelefrequencies
          # for the several loci, is created. It will encompass the actual
          # date at the end.

write.table(as.data.frame(as.matrix(Dest.Chao.all.pops[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(Dest.Chao.all.pops[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(List[[1]])),file=filename.for.allelefrequencies, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(List[[2]])),file=filename.for.sample.sizes, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

          # The tables are saved in the working directory. Two or more analysis
          # that are carried under the same working directory, are saved in the
          # under the same file name. The data are combined in the order as they
          # have been analysed. No data is lost due to be overwritten.
          
cat("\n","----------------------------------------------------------------------------------------------------","\n","\n")         
print(paste("This analysis was carried out on",Sys.time()))

cat("\n","\n","Dest.Chao values for pairwise comparisons between populations for every locus are saved in ","'",filename.for.loci,"'","\n",sep="")
cat("\n","Mean Dest.Chao values over all loci for pairwise comparisons between populations saved in ","'",filename.for.mean,"'","\n",sep="")
cat("\n","Table with allelefrequencies is saved in ","'",filename.for.allelefrequencies,"'","\n",sep="")
cat("\n","Per locus calculated sample sizes are saved in ","'",filename.for.sample.sizes,"'","\n","\n","\n",sep="")  
cat("\n","====================================================================================================","\n")        

          # User information about the end date of the analysis and the filenames
          # under which the several tables have been saved.
  
  
}  
}
  
    




