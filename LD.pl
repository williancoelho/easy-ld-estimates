#!/usr/bin/perl -w

use strict;
use Cwd 'abs_path';
use File::Copy;

my $map_path = $ARGV[0];
my $chr = $ARGV[1];
my $distance_ld = $ARGV[2];
my $distance_pp = $ARGV[3];
my $run_sparse = $ARGV[4];
my $sparse_value = "";
my $counter=5;
if($run_sparse ne "FALSE"){
    $sparse_value = $ARGV[$counter];
    $counter++;
}
else{
    $sparse_value = "2,6,10";
}
my $out1 = $ARGV[$counter];
$counter++;
my $out2 = $ARGV[$counter];
$counter++;
my $out3 = $ARGV[$counter];
$counter++;
my $out4 = $ARGV[$counter];
$counter++;
my $out5 = $ARGV[$counter];
$counter++;
my $out6 = $ARGV[$counter];
$counter++;
my $out7 = $ARGV[$counter];
$counter++;
my $out8 = $ARGV[$counter];
$counter++;
my $out9 = $ARGV[$counter];
$counter++;
my $out10 = $ARGV[$counter];
$counter++;
my $out11 = $ARGV[$counter];
$counter++;


my $breed="\"".$ARGV[$counter]."\"";
$counter++;
while($ARGV[$counter] ne "end_breed_names"){
    $breed=$breed.",\"".$ARGV[$counter]."\"";
    $counter++;
}
$counter++;
my $breed_file="\"".$ARGV[$counter]."\"";
$counter++;
while($ARGV[$counter] ne "end_breed_files"){
    $breed_file=$breed_file.",\"".$ARGV[$counter]."\"";
    $counter++;
}
my $allbreeds=$breed;
my $allbreeds_files=$breed_file;
#print $allbreeds;
$counter++;



#--------------------------------------------------------------Fixing __dq__ problem in target_dist, sparsing and target
# This problem doesn't happen on mac, however when one runs this script on our galaxy server it happens. No idea why.

$distance_ld =~ tr/\__dq__//d; #remove quote
$distance_pp =~ tr/\__dq__//d; #remove quote
$sparse_value =~ tr/\__dq__//d; #remove quote

#--------------------------------------------------------------End fixing __dq__ problem




#--------------------------------------------------------------breeds
my @breed_array = split(/,/, $allbreeds);
my $breed_number = @breed_array;
my @breed_file_array = split(/,/,$allbreeds_files);

$counter=0;
my $LDbycorrelation="";

foreach my $breed(@breed_array)
{
    $breed =~ tr/\"//d; #remove quotes
    mkdir(abs_path()."/".$breed,0777) || die "can't do mkdir!";
    $LDbycorrelation=abs_path()."/".$breed."/LDbycorrelation.R";
    
    #descomentar print "starting ".$breed."... ";
    
    #---------------cria LDbyCorrelation
    my $aux ='
    
				suppressMessages(library(ff))
				
				map_path="'.$map_path.'"
    
				chr<-1:'.$chr.'
    
    ##Functions used in the code##
				refer<-function(phased) {
                    b<-sort(unique(phased))[[1]]
                    recoded<-ifelse(phased==b, 1, 0)
                    return(recoded)
                }
    
    filenam= '.$breed_file_array[$counter].'
    pha<-read.table(filenam, header=T)
    map<-read.table(map_path, header=F)
    
    
    for (i in chr) {
        snps_chr <- map[map[,2]==i,1]
        phased<-as.matrix(pha[pha[,2] %in% snps_chr,-c(1,2)])
        ph_mat<-t(apply(phased, 1, refer))
        r<-cor(t(ph_mat))
        rownames(r)<-pha[pha[,2] %in% snps_chr,2]
        colnames(r)<-rownames(r)
        filenam<-paste("'.abs_path().'/'.$breed.'/SNP.txt",sep="")
        write.table(rownames(r),  filenam, row.names=F, col.names=F, append=T, quote=F)
        idx<-match(rownames(r),map[,1]) #index for SNP position in map
        pos<-map[idx,3]                 #snp position
        for (m in 2:ncol(r)) {           #loop over colums, save upper triangle of matrix
            d<-abs(pos[1:(m-1)]-pos[m])  #distance between snp
            filenam<-paste("'.abs_path().'/'.$breed.'/LD.txt", sep="")
            out_vec<-r[1:(m-1),m]        #output vector: (upper half of) column m of correlation matrix
            write.table(cbind(rownames(r)[1:(m-1)], colnames(r)[m], out_vec, d , i), filenam, col.names=F, row.names=F, quote=F, sep="\t", append=T)
        }
        print(paste("DONE WITH CHROMOSOME:", i, sep=" "))
    }
    ld_long<-read.table.ffdf(file="'.abs_path().'/'.$breed.'/LD.txt",colClasses=c("factor","factor","numeric","integer","integer"))
    ffsave("ld_long",file="'.abs_path().'/'.$breed.'/LD.ff")';
    $counter++;
    
    open(my $fh1, '>', $LDbycorrelation) or die "Can't build file $!";
    print $fh1 "$aux";
    close $fh1;
    
    call_R_func();
    #descomentar print("Done.\n")
}

sub call_R_func {
				my $execute = `Rscript $LDbycorrelation`;
				
}
#descomentar print "Done all breeds!";
#-------------------------------------------------------------------------------------------------------------------------------End Breeds


#-------------------------------------------------------------------------------------------------------------------------------Input.R
my @colors = ("red","green","blue","black","gold","darkviolet","lightcoral","brown","greenyellow","orange4","pink","snow3","tan", "turquoise4", "purple3");
my $colors_1= '"'.$colors[0].'"';

for my $i ( 1..$breed_number-1 ) {  #create a vector of color depending on the number os breeds
    $colors_1 = $colors_1.',"'.$colors[$i].'"';
}

$breed_number = fact($breed_number)/(fact($breed_number-2)*fact(2));

sub fact {
    my ($n) = @_;
    my $prod = 1;
    $prod *= $n-- while $n > 0;
    return $prod;
}

my $colors_2= '"'.$colors[0].'"';

for my $i ( 1..$breed_number-1 ) {  #create a vector of color depending on the number os breeds
    $colors_2 = $colors_2.',"'.$colors[$i].'"';
}

# breed_number!/((breed_number - combinations)! * (combinations)!)

#descomentar print "\nCreating input file... ";

$distance_ld =~ s/X/\,/g;
$sparse_value =~ s/X/,/g;
$distance_pp =~ s/X/,/g;

my $aux ='

map_path<-"'.$map_path.'"

chr<-1:'.$chr.'

breed_fold<-c('.$allbreeds.')

breed_name<-c('.$allbreeds.')

cols<-c('.$colors_1.')

target_dist<-c('.$distance_ld.')
#vetor de parametros n fixo

plot_by_chr<-("Average_LD_by_chr.pdf")

plot_by_dist<-("Average_LD_by_dist.pdf")

run_sparse<-'.$run_sparse.'

sparsing<-c('.$sparse_value.')
#vetor de parametros n fixo

cols2<-c('.$colors_2.')


target<-c('.$distance_pp.')
#vetor de parametros n fixo

persist_out<-"PersistenceofPhase_by_distance.txt"
phase_out<-"PercentPhasePersistence.txt"

pdf_out<-"PersistenceofPhase.pdf"

file_names<-"LD.ff"
snp_file<-"SNP.txt"

';

open(my $fh2, '>',"input_information.R") or die "Não foi possível criar o arquivo $!";
print $fh2 "$aux";
close $fh2;

#descomentar print "Done.";

#--------------------------------------------------------------------------------------------------------------------------------End Input File

#descomentar print "\nstarting copying files... ";

copy ("../../../../../tools/Labegen/LDTool/ld_functions.R",abs_path());
copy ("../../../../../tools/Labegen/LDTool/ld_average.R",abs_path());
copy ("../../../../../tools/Labegen/LDTool/PersistenceofPhase.R",abs_path());
copy ("../../../../../tools/Labegen/LDTool/ld_EffectivePopSize.R",abs_path());


my $PersistenceofPhase = "";

$PersistenceofPhase = abs_path()."/PersistenceofPhase.R";
#print "\n".$PersistenceofPhase."\n";
my $ldAverage = "";

$ldAverage = abs_path()."/ld_average.R";
#print "\n".$ldAverage."\n";

my $ld_EffectivePopSize = "";

$ld_EffectivePopSize = abs_path()."/ld_EffectivePopSize.R";
#print "\n".$ld_EffectivePopSize."\n";

#print "\n diretorio corrente=".abs_path()."\n";
#descomentar print "Done.";
#--------------------------------------------------------------------------------------------------------------------------------End Copying


#descomentar print "\nstarting calculating ld average... ";

call_R_func_ldAverage();

sub call_R_func_ldAverage {
				my $execute = `R --no-save -q < $ldAverage`;
				
}

#descomentar print "Done.";


#descomentar print "\nstarting calculating persistance... ";

call_R_func_Persistence();

sub call_R_func_Persistence {
				my $execute = `R --no-save -q < $PersistenceofPhase`;
				
}

#descomentar print "Done.";


#descomentar print "\nstarting calculating ld EffectivePopSize... ";

call_R_func_ld_EffectivePopSize();

sub call_R_func_ld_EffectivePopSize {
				my $execute = `R --no-save -q < $ld_EffectivePopSize`;
				
}

#descomentar print "Done.";

#descomentar print "\nGenerating out files... ";

copy (abs_path()."/Average_LD_by_dist.pdf",$out1);
copy (abs_path()."/Average_LD_by_chr.pdf",$out2);
copy (abs_path()."/Effective_Pop_Size.pdf",$out3);
copy (abs_path()."/PersistenceofPhase.pdf",$out4);
copy (abs_path()."/LD_average.txt",$out5);
copy (abs_path()."/LD_sparse.txt",$out6);
copy (abs_path()."/LD_table1.txt",$out7);
copy (abs_path()."/PercentPhasePersistence.txt",$out8);
copy (abs_path()."/PersistenceofPhase_by_distance.txt",$out9);
copy (abs_path()."/Timesincepopulationdiverged.txt",$out10);
copy (abs_path()."/Ne_by_generation.txt",$out11);
#descomentar print "Done!";

print "\nDone computing ld estimates!!\n";