use strict;
use warnings;
################ Angle Beginning ################
sub Bending{

my ($AngleType2ID_hr,$AngleTypeArray_ar,$topology_ar,$AtomType_ar,$MoleculeID_ar,$ID2Type_hr,$MoleculeHash_hr) = @_;
##provide all bending sets

my $bendCount = 0;
my @bending; #[i,j,k] for atom id
my @bend2type;# [type ids for i, j, k]
my @bend2MolID;# molecule id

my @topology_ar = @$topology_ar;
my @AngleTypeArray_ar = @$AngleTypeArray_ar;
my %AngleType2ID_hr = %$AngleType2ID_hr;
#for (0..$#topology_ar)
#    {
#        #print "$topology_ar[$_]->[0] $topology_ar[$_]->[1] $topology_ar[$_]->[2]\n";
#    }

#[atomid -1][0]: bonded atom number,[atomid -1][1..]: the bonded atom id -1 

for my $i (0..$#{$topology_ar}){ #begin from the first atom to the last one
    my $atomi = $i + 1;#perl array id + 1, making atom id from 1
    my $bondedNoi = ${$topology_ar}[$i][0];#atom No. bonded with $atomi

    for my $j (1..$bondedNoi){#loop over atoms bonded to $atomi 
        my $atomj = ${$topology_ar}[$i][$j] + 1; #lmp atom id       
        my $bondedNoj = ${$topology_ar}[$atomj -1][0];#atom No. bonded with $atomj

            for my $k (1..$bondedNoj){ 
                my $atomk = ${$topology_ar}[$atomj-1][$k] + 1;        
                
                if( $atomi < $atomk){#preventing from growing back to $atomi
                    $bendCount++;
                    push @bending, [$atomi,$atomj,$atomk];
                    push @bend2type,[ $ID2Type_hr->{$atomi},
                                      $ID2Type_hr->{$atomj},
                                      $ID2Type_hr->{$atomk} ];
                    #???                  
                    push @bend2MolID,[ $MoleculeHash_hr->{$atomi},
                                      $MoleculeHash_hr->{$atomj},
                                      $MoleculeHash_hr->{$atomk} ];
                }
                else{
                    next;            
                }
            }
    }
}

#my %AngleTypeHash;
#my $AngleType = 0;
#for (0..$#AngleTypeArray_ar){
#        chomp;
#        $AngleType++;
#        chomp($AngleType);
#        my $temp = $AngleTypeArray[$_];
#        chomp($temp);
#        $AngleTypeHash{"$temp"}= "$AngleType";
#        my @temp = split("-",$temp);
#        chomp @temp,
#        $AngleTypeHash{"$temp[2]-$temp[1]-$temp[0]"}= "$AngleType";
#    }
#@AngleDone -> Bend ID, type, atomi, atomj,atomk
my @AngleDone;
my @bendType;
#@bending =[bendid -1][i,j,k]
for my $bid (0..$#bending){
    #sort
    
    my $temp = join("-",@{$bend2type[$bid]});
    chomp $temp;
         #if (keys %{{ map {$_, 1} @{$bend2MolID[$bid]} }} == 1){#one molecule ID
            push @AngleDone, [ $bid+1,$AngleTypeHash{$temp},
                $bending[$bid][0],$bending[$bid][1],$bending[$bid][2] ];
            push  @bendType,$AngleTypeHash{$temp};   
         
        # }
} 
my $totbendTypes = max(@bendType);

return (\@AngleDone);#(lammps data for bending)

}
############ Angle Done ############
1;