suppressMessages(library(ff))
#---------------------------------------------
#Script: compute effective population size from Ldbycorrelation.R output
#Derived form ld_average.R of J.P.S. and Y.B. Summer 2011 by F.F.Cardoso 2013
#Use as you see fit, please keep this header in every file using part or all of this code.
#----------------------------------------------

#--------------------------------------------------
## outout files -
# A) Table with Ne at different generation intervals given by the "gen" vector. 1 column per population.
#       >Ne_by_generation.txt, rows:
#       interval of t=1/2c values to which Ne was calculated using average r2 and c in the bin 
# B) Plot of Ne by generation Table with average LD between adjacent markers in each chromosome. 1 row per chromosome, 1 column per population.
#      > LD_chr.txt
# E) Table with average LD at different intervals (rows). 1 column per population.
#       >LD_target_distances.txt, rows:
#       average LD at distances specified in vector "target" above
#--------------------------------------------------------

x<-getwd();
x <- paste(x,"input_information.R",sep="/")
x
source(x)
n_chr<-length(chr)
n_breeds<-length(breed_fold)
source("./ld_functions.R")

# Function to estimate average Ne for consecutive generations bins given by a vector
Calc_Ne_average<-function(ld_matrix, gen_vector) {
    Ne_out<-c()
   # Calculate the generation (t=1/2*c) corresponding to each marker pair distance (c) 
    tgen <- 100000000/(2*ld_matrix[,4])
   ## histogram of marker pairs per generation bin
   pdf(paste("./", breed_fold[k], "/", "hist_markerpairs_Ne_bin.pdf", sep=''))
   hist_gen<-hist(tgen[tgen<=max(gen_vector)], breaks=gen_vector, xlab="Generation")
   dev.off()
   hist_gen<-cbind(as.data.frame(hist_gen$breaks)[2:nrow(as.data.frame(hist_gen$breaks)),],
         as.data.frame(hist_gen$counts))
   write.table(hist_gen, paste("./", breed_fold[k], "/", "hist_markerpairs_Ne_bin.txt", sep=''), col.names=T, row.names=T, quote=F, append=F, sep='\t')

    for (m in 1:(length(gen_vector)-1)){
        if(gen_vector[m+1]>gen_vector[m])
        {
   # Select just r2 and c values in the desired bin range
        idx<-(tgen>=gen_vector[m])*(tgen<gen_vector[m+1])
        r2 <- ld_matrix[idx==1,3]^2
        c <- ld_matrix[idx==1,4]/100000000 
   # Limit r2 values since Ne = Inf for r2=0 and Ne = 0 for r2=1 
        idx2 <- r2 > 0.0 & r2 < 1.0 
   # Calculation of Ne averaging r2 and c  
        Ne_out<-c(Ne_out, (1-mean(r2))/(4*mean(r2)*mean(c)))
        } }
    return(Ne_out)
    }


for (p in 1:length(breed_fold)) {
    if(file.exists(paste("./", breed_fold[p], "/", file_names, ".ffData", sep=''))==FALSE) {
    print(paste("In breed ", breed_name[p], " the file LD.ff does not exist" , sep=''))
    print(paste("Please run LDbycorrelation.R in the subdirectory of population ", breed_name[p], sep=""))
    stop()
    }
    }

#FC target bins of generation in the past to calculate Ne
gen <- c(seq(0.5,10.5,1),seq(12.5,102.5,5),seq(125,1025,50))

# matrix of Ne by generation
ne_gen<-matrix(99, ncol=n_breeds, nrow=(length(gen)-1))
colnames(ne_gen)<-breed_name
rownames(ne_gen)<-paste(gen[1:(length(gen)-1)], "-", gen[2:length(gen)])

## check of breed number, names and folders
if(n_breeds!=length(breed_fold))
   {print("Breed number and the folder names do not match")
   stop(.call=T)
   }
if(length(breed_fold)!=length(breed_name))
   {print("Number of breed names and folders do not match")
   stop(.call=T)
   }

# read the map file 
# map<-read.table(map_path, header=F)

for (k in 1:n_breeds) {
     # the ffload doesn't work and I am not sure why
#     k <- 1
     suppressMessages(file_nam<-paste("./", breed_fold[k], "/", file_names, sep=''))
     suppressWarnings(suppressMessages(ffload(file_nam))) #takes from a few secs to a couple minutes, much less that read.ffdf

     #only consider markers with pairwise distance < 100Mb (100Mb ~ 1Morgan)
    suppressMessages(ld_long<-as.ffdf(ld_long[ld_long[,4]<=100000000,] ))
    print(paste("In population ", breed_name[k], " the matrix length is: ", nrow(ld_long), sep=''))  # again checking that the input used by this program makes sense in your case

## FC - calculate Ne by generation - fill ne_gen
    ne_gen[,k]<-Calc_Ne_average(ld_long, gen)

    close(ld_long)
    rm(ld_long)   # we have to do this, otherwise the next ffload will not overwrite the existing ld_long objects
    }


## save the table constructed
write.table(ne_gen, "Ne_by_generation.txt", col.names=T, row.names=T, quote=F, append=F, sep='\t')

## plotting the Ne by generation up "maxgenplot"
plot_ne_by_gen <- c("Effective_Pop_Size.pdf")
pdf(plot_ne_by_gen)

# Calculates the generation (t=1/2*c) corresponding to each bin rounded to nearest integer     
nbin <- length(gen)-1
tgen <- round((gen[1:nbin]+gen[2:(nbin+1)])/2, digits = 0)
maxgenplot<-50
idx<-tgen <= maxgenplot

plot(0, xlim=c(0,max(tgen[idx])), ylim=c(0,max(ne_gen[idx,])*1.2), col="white", xlab="Generation in the past, years", ylab="Ne")
legend(x="topleft", bty="n",legend=breed_name, col=cols, pch=18)

for (k in 1:n_breeds) {
points(tgen[idx],ne_gen[idx,k], col=cols[k], pch=18, type="b")
}

dev.off()


