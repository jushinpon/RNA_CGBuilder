use strict;
use warnings;
use feature qw(say);
use List::Util qw(min max sum);
################ Torsion Beginning ################
sub Torsion
{

  my ($topology_ar,$AtomType_ar,$MoleculeID_ar,$ID2Type_hr,$MoleculeHash_hr) = @_;
my @torsion;
my $torCount;
my @topology_ar = @$topology_ar;
my @bend2MolID;
my @Torsion2type;
my @torsion2MolID;
#for (0..$#topology_ar)
#    {
#       # print "$topology_ar[$_]->[0] $topology_ar[$_]->[1] $topology_ar[$_]->[2] $topology_ar[$_]->[3]\n";
#    }


for my $i (0..$#{$topology_ar}) #begin from the first atom
    { 
        my $atomi = $i + 1;
        my $bondedNoi = ${$topology_ar}[$i][0];
        for my $j (1..$bondedNoi)
            { 
                my $atomj = ${$topology_ar}[$i][$j] + 1;        
                my $bondedNoj = ${$topology_ar}[$atomj -1][0];
                for my $k (1..$bondedNoj)
                    {                        
                        my $atomk = ${$topology_ar}[$atomj -1][$k] + 1;

                        if($atomk != $atomi){# good to form a reasonable bending so far
                                my $bondedNok = ${$topology_ar}[$atomk -1][0];
                                for my $l (1..$bondedNok) 
                                    {
                                        my $atoml = ${$topology_ar}[$atomk-1][$l] + 1;   

                                        if( $atomj != $atoml && $atomi < $atoml) 
                                            {
                                                $torCount++;
                                                push @torsion, [$atomi,$atomj,$atomk,$atoml];
                                                push @Torsion2type, [ $ID2Type_hr->{$atomi},
                                                                      $ID2Type_hr->{$atomj},
                                                                      $ID2Type_hr->{$atomk},
                                                                      $ID2Type_hr->{$atoml},
                                                                    ];
                                                push @torsion2MolID,[ $MoleculeHash_hr->{$atomi},
                                                                      $MoleculeHash_hr->{$atomj},
                                                                      $MoleculeHash_hr->{$atomk},
                                                                      $MoleculeHash_hr->{$atoml},
                                                                    ];
                                            }
                                       else
                                           {
                                               next;            
                                           }
                                    }#loop l
                        }
                        else{# if atomk == atomi
                          next;# invalid for a torsion  
                        }
                    }# k loop
            }#j loop
    }#i loop


for (0..$#torsion)
    {
       # print "$_ $torsion[$_][0] $torsion[$_][1] $torsion[$_][2] $torsion[$_][3]\n";
    }
#print "Sleeping\n";
#sleep(100);


my @TorsionTypeArray = qw/1-2-3-4 1-2-3-4p 1p-2-3-4 1-2-3-9 1-2-3-9p 1p-2-3-9 1-2-3-5 1-2-3-5p 1p-2-3-5 1-2-3-6 1-2-3-6p 1p-2-3-6 1-2-8-4 1-2-8-4p 1p-2-8-4 1-2-8-6 1-2-8-6p 1p-2-8-6 1-2-8-7 1-2-8-7p 1p-2-8-7 1-2-1-2 2-1-2-3 2-1-2-3p 2p-1-2-3 2-1-2-8 2-1-2-8p 2p-1-2-8 2-1-2-1 2-3-4-9 2-3-5-6 2-3-6-5 2-3-9-4 2-8-4-7 2-8-6-7 2-8-7-4 2-8-7-6/;
my @TorsionPrime = qw/1-2-3-4 4-3-2-1 1-2-3-9 9-3-2-1 1-2-3-5 5-3-2-1 1-2-3-6 6-3-2-1 1-2-8-4 4-8-2-1 1-2-8-6 6-8-2-1 1-2-8-7 7-8-2-1 2-1-2-3 3-2-1-2 2-1-2-8 8-2-1-2/;
my %TorsionPrime = map { $_ => 1 } @TorsionPrime;

my %TorsionTypeHash;
my $TorsionType = 0;
for (0..$#TorsionTypeArray)
    {
        chomp;
        $TorsionType++;
        chomp($TorsionType);
        my $temp = $TorsionTypeArray[$_];
        chomp($temp);
        $TorsionTypeHash{"$temp"}= "$TorsionType";
        my @temp = split("-",$temp);
        chomp @temp,
        $TorsionTypeHash{"$temp[3]-$temp[2]-$temp[1]-$temp[0]"}= "$TorsionType"; # 4-3-2-1
        #$TorsionTypeHash{"$temp[3]-$temp[2]-$temp[1]-$temp[0]p"}= "$TorsionType"; # 4-3-2-1
       # say "$temp[3]-$temp[2]-$temp[1]-$temp[0]";
    }

my @TorsionDone;
my @TorType;

for my $bid (0..$#torsion)
    {
        my $temp = join("-",@{$Torsion2type[$bid]});
        chomp $temp;
       # say "$temp";
        if(exists($TorsionPrime{$temp})) # could have the prime para 
            { 
                if (keys %{{ map {$_, 1} @{$torsion2MolID[$bid]} }} == 1)#one molecule ID
                    {
                        push @TorsionDone , [ $bid+1,$TorsionTypeHash{$temp},
                            $torsion[$bid][0],$torsion[$bid][1],$torsion[$bid][2],$torsion[$bid][3]];
                        push  @TorType,$TorsionTypeHash{$temp};   
                    }
                else # prime included
                    {
                        my $prime = "$temp"."p";
                        #say $prime;
                        push @TorsionDone, [ $bid+1,$TorsionTypeHash{$prime},
                          $torsion[$bid][0],$torsion[$bid][1],$torsion[$bid][2],$torsion[$bid][3] ];
                           # print "$TorsionTypeHash{$prime} $torsion[$bid][0] $torsion[$bid][1] $torsion[$bid][2] $torsion[$bid][3]\n";
                       # push  @TorsionDone,$TorsionTypeHash{$prime}; 

                       # print "Use prime parameter $prime for "."$bid+1 ". "bend ID\n";
                       # print "molecule IDs: @{$bend2MolID[$bid]}\n";        
                    }
            }
        else
            {
                push @TorsionDone, [ $bid+1,$TorsionTypeHash{$temp},
                       $torsion[$bid][0],$torsion[$bid][1],$torsion[$bid][2],$torsion[$bid][3] ];
                       #say "$torsion[$bid][0] $torsion[$bid][1] $torsion[$bid][2] $torsion[$bid][3]";
                push  @TorType,$TorsionTypeHash{$temp};          
            }
    }
#say  Dumper @TorsionDone;
my @GetTorsionType;
  #  print $TorsionDone[0]->[1];
    for (0..$#TorsionDone)
        {
           #say "$TorsionDone[$_]->[0] $TorsionDone[$_]->[1] $TorsionDone[$_]->[2] $TorsionDone[$_]->[3] $TorsionDone[$_]->[4] $TorsionDone[$_]->[5]";
           push @GetTorsionType,"$TorsionDone[$_]->[1]";
        }

    return (\@TorsionDone,\@GetTorsionType);
################################ Torsion Done
}
1;