# functions for the LD codes

# Function to estimate average LD for consequtive distances given by a vector
Calc_LD_average<-function(ld_matrix, distance_vector) {
    ld_out<-c()
    for (m in 1:(length(distance_vector)-1)){
        if(distance_vector[m+1]>distance_vector[m])
        {
        idx<-(ld_matrix[,4]>=distance_vector[m])*(ld_matrix[,4]<distance_vector[m+1])
        ld_out<-rbind(ld_out, c(mean(ld_matrix[idx==1,3]^2),sd(ld_matrix[idx==1,3]^2)))
        } }
    return(ld_out)
    }


# Function to estimate LD for adjacent markers, % r2> 0.2, % r2>0.3, number of markers used in a breed, average distance between adjacent markers
# will also estimate average LD between adjacent markers per chr
Calc_LD_table1<-function(ld_matrix, snps, snp, n_chr , get_chr) {
    ## calculate average LD for adjacent markers  - fill ld_table1

    idx<-(paste(ld_matrix[,2], ld_matrix[,1], sep='')%in%snps)+(paste(ld_matrix[,1], ld_matrix[,2], sep='')%in%snps)
    ld_table1<-c( mean(ld_matrix[idx==1,3]^2),sd(ld_matrix[idx==1,3]^2),
                  (sum(ld_matrix[idx==1,3]^2>0.2)/nrow(ld_matrix[idx==1,]))  ,
                  (sum(ld_matrix[idx==1,3]^2>0.3)/nrow(ld_matrix[idx==1,]))  ,
#                  (length(snps)+n_chr)                                                    ,
                   (length(snps)+1)                                                    ,
                  mean(ld_matrix[idx==1,4])                                             )
    if(get_chr==1) {
    ld_table1<-c(ld_table1,by(ld_matrix[idx==1,3]^2, ld_matrix[idx==1,5], mean),
                           by(ld_matrix[idx==1,3]^2, ld_matrix[idx==1,5], sd))
    }
    return(ld_table1)
    }

# Function that will assemble the SNP list for the sparse set of markers
Render_snps_sparse<-function(map, sparsing, snp, n_chr) {
    idx<-map[,1]%in%snp[,1]
    map_new<-map[idx,]
    for (i in 1:n_chr) {
        map_i<-map_new[map_new[,2]==i,]
        map_i<-map_i[order(map_i[,3]),]
        id_in<-seq(1, nrow(map_i), sparsing)
        snp_out<-map_i[id_in,1]
        if(i==1) {
        snps<-data.frame(snp_out[1:(length(snp_out)-1)],snp_out[2:(length(snp_out))])
        } else {
        snps<-rbind(snps, data.frame(snp_out[1:(length(snp_out)-1)],snp_out[2:(length(snp_out))]))
        }
        }
    return(snps)
    }







