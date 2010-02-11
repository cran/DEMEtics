D<-function(Hs,Ht,values){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           Hs,Ht,values <- D.calc();
#------------------------------------------------------------------------------------------------------------------------------  

          # A function that calculates the D value.
          # The arguments are the Hs-value, the Ht-value and the vector 'values'
          # (a vector that contains the number of individuals that have been 
          # sampled from each population - the sample sizes).
          
          # The Hs- and Ht-values can be calculated by the functions 'Hs()' and 'Ht()'.
           
          # See: Jost L. (2008). Gst and its relatives do not measure differentiation. 
          # Molecular Ecology 17,4015-4026.


          n<-length(values)

                    # n is the number of populations that have been sampled.

          harmonic<-function(values){
                                      n/
                                      sum(1/values)
                                      }
                                      
                    # The harmonic mean is calculated from the sample sizes...
                                                           
          N<-harmonic(values)
         
                    # ... and ascribed to the variable N.


          Hs.est<-function(values,Hs){
                                      (2*N/(2*N-1))*Hs
                                      }
                                                       
          Hs.est<-Hs.est(values,Hs)
          
                    # The Hs.est values are calculated and ascribed to the variable 'Hs.est'.
                    # See: Jost L. (2008). Gst and its relatives do not measure differentiation. 
                    # Molecular Ecology 17,4015-4026.  

          Ht.est<-function(values,n,Ht,Hs){
                                            Ht+ (Hs.est/(2*N*n))
                                            }
          Ht.est<-Ht.est(values,n,Ht,Hs)
          
                    # The Ht.est values are calculated and ascribed to the variable 'Ht.est'.
                    # See: Jost L. (2008). Gst and its relatives do not measure differentiation. 
                    # Molecular Ecology 17,4015-4026. 

          ((Ht-Hs)/
                          (1-Hs))*
                                      (n/(n-1))
                                      
                    # Out of Ht.est, Hs.est and n, the D value is calculated.                                      
                                                  
}
