suppressMessages(library(ff))
#---------------------------------------------
#Script: compute Persistence of phase between breeds from Ldbycorrelation.R output
#by J.P.S. and Y.B. Summer 2011.
#Use as you see fit, please keep this header in every file using part or all of this code.
#----------------------------------------------

x<-getwd();
x <- paste(x,"input_information.R",sep="/")
x
source(x)
n_chr<-length(chr)
n_breeds<-length(breed_fold)

if(length(target)!=0) {
   dist_all<-c(target, seq(0,5000000,50000))  } else {
   dist_all<-distances}

n_int<-dist_all[2:length(dist_all)]>dist_all[1:(length(dist_all)-1)]

row_n<-paste(dist_all[1:(length(dist_all)-1)]/1000000, "-", dist_all[2:length(dist_all)]/1000000, sep='')[n_int]

col_n<-paste(t(combn(breed_name, 2))[,1], t(combn(breed_name, 2))[,2], sep='-')
n_col<-(factorial(n_breeds)/(factorial(2)*factorial(n_breeds-2)))

persist<-matrix(99, ncol=(n_col+1), nrow=sum(n_int))
colnames(persist)<-c(col_n, "n_snp")
rownames(persist)<-row_n

phase<-matrix(99, ncol=n_col, nrow=sum(n_int))
colnames(phase)<-col_n
rownames(phase)<-row_n


snp_in<-read.table(paste("./", breed_fold[1],"/", snp_file, sep=''), header=F)
for (k in 2:n_breeds) {
    snp_i<-read.table(paste("./", breed_fold[k],"/", snp_file, sep=''), header=F)
    snp_in<-merge(snp_in, snp_i)
    }

nrow(snp_in)          #22340


map<-read.table(map_path, header=F)


for (k in 1:n_breeds) {

     suppressMessages(file_nam<-paste("./", breed_fold[k], "/", file_names, sep=''))
     suppressWarnings(suppressMessages(ffload(file_nam))) #takes from a few secs to a couple minutes, much less that read.ffdf

    suppressMessages(idx<-(ld_long[,1]%in%snp_in[,1])&(ld_long[,2]%in%snp_in[,1])&(ld_long[,4]<=10000000))

    if(k==1) {
    suppressMessages(ld_all<-as.ffdf(ld_long[idx,c(1,2,4,3)]))
    } else {

    suppressMessages(ld_all<-as.ffdf(data.frame(ld_all, ld_long[idx,3])))
    }
    rm(ld_long)
    }

    suppressMessages(print(dim(ld_all)))

columns<-t(combn(c(4:(3+n_breeds)), 2))

for (w in 1:nrow(columns)) {
    per<-c()
    pha<-c()
    n_snp<-c()
    for (m in 1:(length(dist_all)-1)){
        if(dist_all[m+1]>dist_all[m])
        {
        idx<-(ld_all[,3]>=dist_all[m])*(ld_all[,3]<dist_all[m+1])
        ld_m<-ld_all[idx==1,]
        per<-c(per, cor(ld_m[,columns[w,1]], ld_m[,columns[w,2]]))
        pha<-c(pha, (1-sum((((ld_m[,columns[w,1]]<=0))*((ld_m[,columns[w,2]]<=0)))+(((ld_m[,columns[w,1]]>0))*((ld_m[,columns[w,2]]>0))))/nrow(ld_m)))
        n_snp<-c(n_snp, sum(idx))
        } }
        persist[,w]<-per
        phase[,w]<-pha
          }
    persist[,(n_col+1)]<-n_snp

suppressMessages(write.table(persist, persist_out, col.names=T, row.names=T, quote=F, append=F, sep='\t'))

suppressMessages(write.table(phase, phase_out, col.names=T, row.names=T, quote=F, append=F, sep='\t'))

pdf(pdf_out)

    plot(0, axes=T, ylim=c(0,1.0), xlim=c(0,10000), col="white", type="b", pch=18, xlab="Marker distance in kb", ylab="Correlation of Phase")
    legend(x="topright", bty="n",legend=col_n, pch=18, col=cols)

    dist<-(seq(0,10000000,100000)[2:101])/1000-50

    for (l in 1:(n_breeds*(n_breeds-1)/2)) {
        points(persist[(length(target)):nrow(persist),l]~dist, col=cols[l], type="b", pch=18)
        }
dev.off()

init<-10000
fin<-300000
stp<-2500

ld<-ld_all[(ld_all[,3]<=fin & ld_all[,3]>init),]

combin<-columns
dist<-seq(init,fin,stp)

output<-matrix(99, ncol=7, nrow=(length(dist)-1))

for (i in 1:(length(dist)-1)) {
    index<-(ld[,3]>=dist[i])*(ld[,3]<dist[i+1])
 
    output[i,7]<-mean(ld[index==1, 3])/100000000
 
    for (k in 1:(n_breeds*(n_breeds-1)/2)) {
        output[i,k]<-cor(ld[index==1,combin[k,1]],ld[index==1,combin[k,2]])
        }
    }

t<-c()

for (k in 1:(n_breeds*(n_breeds-1)/2)) {
    y<-log(output[,k])
    model<-lm(y~output[,7])
    print(c(exp(coef(model)[1]), coef(model)[2]/(-2)))
    t<-append(t, coef(model)[2]/(-2), after=length(t))
       }

t_out<-data.frame(col_n, round(t, 0))
colnames(t_out)<-c("Populations", "Time since divergence")
write.table(t_out, "Timesincepopulationdiverged.txt", col.names=T, row.names=F, append=F, quote=F, sep='\t')

