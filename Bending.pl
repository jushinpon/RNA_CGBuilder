use strict;
use warnings;
################ Angle Beginning ################
sub Bending
    {
        my ($topology_ar,$AtomType_ar,$MoleculeID_ar,$ID2Type_hr,$MoleculeHash_hr) = @_;
my $bendCount = 0;
my @bending; #atom id
my @bend2type;# type id
my @bend2MolID;# molecule id

my @topology_ar = @$topology_ar;
#for (0..$#topology_ar)
#    {
#        #print "$topology_ar[$_]->[0] $topology_ar[$_]->[1] $topology_ar[$_]->[2]\n";
#    }

for my $i (0..$#{$topology_ar}){ #begin from the first atom
    my $atomi = $i + 1;
    my $bondedNoi = ${$topology_ar}[$i][0];

    for my $j (1..$bondedNoi){ 
        my $atomj = ${$topology_ar}[$i][$j] + 1;        
        my $bondedNoj = ${$topology_ar}[$atomj -1][0];

            for my $k (1..$bondedNoj){ 
                my $atomk = ${$topology_ar}[$atomj-1][$k] + 1;        
                
                if( $atomi < $atomk){
                    $bendCount++;
                    push @bending, [$atomi,$atomj,$atomk];
                    push @bend2type,[ $ID2Type_hr->{$atomi},
                                      $ID2Type_hr->{$atomj},
                                      $ID2Type_hr->{$atomk} ];
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

#for (0..$#bending)
#    {
#       # print "@bending[$_]->[0] @bending[$_]->[1] @bending[$_]->[2]\n";
#    }
=Angle
1-2(P-S) => 1-2-1 1-2-3、1-2-3p、1-2-8、1-2-8p、2-1-2
2-3、2-8(S-CG、S-CU) => 2-3-4、2-3-5、2-3-6、2-3-9、2-8-4、2-8-6、2-8-7
3-4、3-6、8-4、8-6(CG-N6、CG-O6、CU-N6、CU-O6) =>  3-4-9、3-6-5、8-4-7、8-6-7
5-3、9-3、7-8(N2-CG、CA-CG、O2-CU) => 5-3-6、9-3-4、7-8-4、7-8-6
4-7(N6-O2) => 4-7-8
4-9(N6-CA) => 4-9-3
6-5(O6-N2) => 6-5-3
6-7(O6-O2) => 6-7-8
=cut
my @AngleTypeArray = qw/1-2-1 1-2-3 1-2-3p 1p-2-3 1-2-8 1-2-8p 1p-2-8 2-1-2 2-3-4 2-3-5 2-3-6 2-3-9 2-8-4 2-8-6 2-8-7 3-4-9 3-6-5 5-3-6 8-4-7 9-3-4 7-8-4 7-8-6 8-6-7 4-7-8 4-9-3 6-5-3 6-7-8/;
my @AnglePrime = qw/1-2-3 3-2-1 1-2-8 8-2-1/;
my %AnglePrime = map { $_ => 1 } @AnglePrime;

my %AngleTypeHash;
my $AngleType = 0;
for (0..$#AngleTypeArray)
    {
        chomp;
        $AngleType++;
        chomp($AngleType);
        my $temp = $AngleTypeArray[$_];
        chomp($temp);
        $AngleTypeHash{"$temp"}= "$AngleType";
        my @temp = split("-",$temp);
        chomp @temp,
        $AngleTypeHash{"$temp[2]-$temp[1]-$temp[0]"}= "$AngleType";
    }
#@AngleDone -> Bend ID, type, atomi, atomj,atomk
my @AngleDone;
my @bendType;

for my $bid (0..$#bending){
    my $temp = join("-",@{$bend2type[$bid]});
    chomp $temp;
    if(exists($AnglePrime{$temp})) { # could have the prime para 

         if (keys %{{ map {$_, 1} @{$bend2MolID[$bid]} }} == 1){#one molecule ID
            push @AngleDone, [ $bid+1,$AngleTypeHash{$temp},
                $bending[$bid][0],$bending[$bid][1],$bending[$bid][2] ];
            push  @bendType,$AngleTypeHash{$temp};   
         }
         else{# prime included
            my $prime = "$temp"."p";
            push @AngleDone, [ $bid+1,$AngleTypeHash{$prime},
                $bending[$bid][0],$bending[$bid][1],$bending[$bid][2] ];
            push  @bendType,,$AngleTypeHash{$prime};  
               # print "Use prime parameter $prime for "."$bid+1 ". "bend ID\n";
               # print "molecule IDs: @{$bend2MolID[$bid]}\n";        
         }
     }
     else{
         push @AngleDone, [ $bid+1,$AngleTypeHash{$temp},
                $bending[$bid][0],$bending[$bid][1],$bending[$bid][2] ];
         push  @bendType,$AngleTypeHash{$temp};   
                
     }
}

my $totbendTypes = max(@bendType);

return (\@AngleDone,\$totbendTypes);#(lammps data for angle, )

}
############ Angle Done ############
1;