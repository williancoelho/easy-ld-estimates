suppressMessages(library(ff))
options(error = expression(NULL))
#---------------------------------------------
#Script: compute LD descriptive statistics from Ldbycorrelation.R output
#by J.P.S. and Y.B. Summer 2011.
#Use as you see fit, please keep this header in every file using part or all of this code.
#----------------------------------------------

#--------------------------------------------------
## outout files -
# A) average LD between adjacent markers for four breeds  (pdf)
#>plot_by_chr
# average LD for markers at various distances within a 50kb window around a target distance   (pdf)
#>plot_by_dist
# B) Table with LD statistics with 1 column per population and one row per statistic.
#      >LD_table1.txt, statistics:
#      1) average LD between adjacent markers
#      2) % adjacent marker pairs with r2>0.2
#      3) % adjacent marker pairs with r2>0.3
#      4) number of markers included in calculations
#      5) average distance between adjacent markers
# C) Table with average LD at different intervals (rows). 1 column per population.
#       >LD_average.txt, rows:
#       average LD at distances specified in vector "distances" above
# D) Table with average LD between adjacent markers in each chromosome. 1 row per chromosome, 1 column per population.
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

for (p in 1:n_breeds) {
    if(file.exists(paste("./", breed_fold[p], "/", file_names, ".ffData", sep=''))==FALSE) {
    print(paste("In breed ", breed_name[p], " the file LD.ff does not exist" , sep=''))
    print(paste("Please run LDbycorrelation.R in the subdirectory of population ", breed_name[p], sep=""))
    stop()
    }
    }

if(length(target_dist)!=0) {
   distances<-c(target_dist, seq(0,10000000,20000))  } else {
   distances<-distances}

n_int<-distances[2:length(distances)]>distances[1:(length(distances)-1)]
# matrix of LD by distance for all intervals from 0-10Mb
# FC sd included
ld_average<-matrix(99, ncol=n_breeds*2, nrow=(sum(n_int)))
colnames(ld_average)<-paste(rep(breed_name,each=2),c("mean", "sd"))
rownames(ld_average)<-paste(distances[1:(length(distances)-1)]/1000000, "-", distances[2:length(distances)]/1000000, sep='')[n_int]
# Table of LD for adjacent markers, % r2> 0.2, % r2>0.3, number of markers used in a breed, average distance between adjacent markers
# will also include average LD between adjacent markers per chr
ld_table1<-matrix(99, ncol=n_breeds, nrow=(6+n_chr*2))
colnames(ld_table1)<-breed_name
## FC last field in row names should be 1:n_chr
rownames(ld_table1)<-c("Adj. marker LD - Mean","Adj. marker LD - SD", "% >0.2", "% >0.3", "# markers", "Average marker distance", 
                       paste(1:n_chr,"Mean"),paste(1:n_chr,"SD"))
# Table for LD between adjacent markers in sparse sets - if computation was wanted
 if(run_sparse==T) {
    ld_sparse<-matrix(99, ncol=(n_breeds+1), nrow=(6*length(sparsing)))
    colnames(ld_sparse)<-c(breed_name, "sparsing")
    rownames(ld_sparse)<-rep(c("Adj. marker LD - Mean","Adj. marker LD - SD", "% >0.2", "% >0.3", "# markers", "Average marker distance"), length(sparsing))
    }

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
map<-read.table(map_path, header=F)

for (k in 1:n_breeds) {
     # the ffload doesn't work and I am not sure why
     suppressWarnings(suppressMessages(file_nam<-paste("./", breed_fold[k], "/", file_names, sep='')))
     suppressWarnings(suppressMessages(ffload(file_nam))) #takes from a few secs to a couple minutes, much less that read.ffdf

 #    str(ld_long)     # this will output the structure of your ff-object and can be used in case troubleshooting is necessary

    #only consider markers with pairwise distance <10Mb
    suppressMessages(ld_long<-as.ffdf(ld_long[ld_long[,4]<=10000000,] ))
# FC 
## histogram of marker distances
   pdf(paste("./", breed_fold[k], "/", "hist_distances.pdf", sep=''))
   hist_distances<-hist(ld_long[,4], breaks=10000, xlab="Marker distances",xlim=c(0,100000),)
   dev.off()
   hist_distances<-cbind(as.data.frame(hist_distances$breaks)[2:nrow(as.data.frame(hist_distances$breaks)),],
         as.data.frame(hist_distances$counts))
   write.table(hist_distances, paste("./", breed_fold[k], "/", "hist_distances.txt", sep=''), col.names=T, row.names=T, quote=F, append=F, sep='\t')
# FC
	

    print(paste("In population ", breed_name[k], " the matrix length is: ", nrow(ld_long), sep=''))  # again checking that the input used by this program makes sense in your case

    ## calculate average LD for markers by distance - fill ld_average
    # make a vector that will save the results of the average distance for now
   
    ld_average[,(k*2-1):(k*2)]<-Calc_LD_average(ld_long, distances)

    snp<-read.table(paste("./", breed_fold[k], "/",snp_file, sep=''), header=F)
    snps<-paste(snp[1:(nrow(snp)-1),1], snp[2:nrow(snp),1], sep='')
    
    #aqui
    #erro
    #Error in ld_table1[, k] <- Calc_LD_table1(ld_long, snps, snp, n_chr, 1) : número de itens para para substituir não é um múltiplo do comprimento do substituto
    
    ld_table1[,k]<-Calc_LD_table1(ld_long, snps, snp, n_chr,  1)

    ## estimation of LD for sparse marker sets
    if(run_sparse==T) {
        for (g in 1:length(sparsing)) {
            snps<-Render_snps_sparse(map, sparsing[g], snp, n_chr)
            snps<-paste(snps[,1], snps[,2], sep='')
            ld_sparse[seq(1,nrow(ld_sparse),  6)[g]:(seq(1,nrow(ld_sparse),  6)[g]+5),k]<-Calc_LD_table1(ld_long, snps, snp, n_chr, 0)
            ld_sparse[seq(1,nrow(ld_sparse),  6)[g]:(seq(1,nrow(ld_sparse),  6)[g]+5),k+1]<-sparsing[g]
            }
        }
    close(ld_long)
    rm(ld_long, snp, snps)   # we have to do this, otherwise the next ffload will not overwrite the existing ld_long objects
    }

    ## at this point we should have ld_chr, ld_target, ld_average and ld_table1 so that we can plot from the matrixes 

## save the tables constructed
write.table(ld_average, "LD_average.txt", col.names=T, row.names=T, quote=F, append=F, sep='\t')
write.table(ld_table1, "LD_table1.txt", col.names=T, row.names=T, quote=F, append=F, sep='\t')
if(run_sparse==T) {
write.table(ld_sparse, "LD_sparse.txt", col.names=T, row.names=T, quote=F, append=F, sep='\t')
} #Willian Coelho - this if condition was included to avoid erros when run_sparse is FALSE
## plotting the average by distance and chr
pdf(plot_by_dist)

ylim<-round(max(ld_average[(length(target_dist)):nrow(ld_average),seq(1,k*2-1,2)])+0.1,1)

plot(0, xlim=c(0,10000), ylim=c(0,ylim), col="white", xlab="Marker distance in kb", ylab=expression(Average*" "*r^2)  )
legend(x="topright", bty="n",legend=breed_name, col=cols, pch=18)

for (k in 1:n_breeds) {
dist<-distances[n_int]
dist<-(dist[(length(target_dist)+1):length(dist)])/1000-50
points(ld_average[(length(target_dist)):nrow(ld_average),(k*2-1)]~dist, col=cols[k], pch=18, type="b")
}

dev.off()

## plotting the average by distance and chr
pdf(plot_by_chr)

ylim<-round(max(ld_table1[7:(6+n_chr),])+0.2,1)

plot(0, xlim=c(0,n_chr), ylim=c(0,ylim), col="white", xlab="Chromosome", ylab="Average r2")
legend(x="topright", bty="n",legend=breed_name, col=cols, pch=18)

for (k in 1:n_breeds) {
points(ld_table1[7:(6+n_chr),k]~chr, col=cols[k], pch=18, type="p")
}

dev.off()
