#use List::MoreUtils qw(mesh); 
use 5.010;
use strict;
use warnings;
#use diagnostics;
use List::Util qw(min max sum);
use List::MoreUtils 'pairwise';
use Data::Dumper;
use List::MoreUtils qw(uniq);
#use Parallel::ForkManager;
#use MCE::Shared;
open my $StructureOut , "> RNA.data";

#(typeArray,BondingArray,CoordinationArrays according to bead number in typeArray)
#RNA grows along y dimension
my @growthV = (0,-10.6,0); # the vector to grow your chain,could be variable
my @BaseA = ([2,3,4,9],["1-2","2-3","3-4","4-2"],[0,0,0],[3.74,0,0],[8.659,1.944,0],[8.659,-1.986,0]);
my @BaseU = ([2,8,4,7],["1-2","2-3","3-4","4-2"],[0,0,0],[3.61,0,0],[7.393,2.467,0],[7.393,-2.483,0]);
my @BaseC = ([2,3,6,5],["1-2","2-3","3-4","4-2"],[0,0,0],[3.74,0,0],[7.920,3.063,0],[7.920,-1.907,0]);
my @BaseG = ([2,8,6,7],["1-2","2-3","3-4","4-2"],[0,0,0],[3.61,0,0],[8.095,2.877,0],[8.095,-4.053,0]);
#my @linker = ([1],[-2.422,-5.3,0]); #P
my @linker = ([1],[0,0,0]); #P
my %blocks = (
    A => \@BaseA,
    U => \@BaseU,
    C => \@BaseC,
    G => \@BaseG,
    Linker => \@linker
);

#Define bond types in lammps (make inverse also ok)
my @BondTypeArray = qw/1-2 2-3 2-8 3-4 5-3 3-6 9-3 4-7 8-4 4-9 6-5 6-7 8-6 7-8/;
my %BondTypeHash;
my $BondType = 0;
for (0..$#BondTypeArray){
    chomp;
    $BondType++;
    chomp($BondType);
    my $temp = $BondTypeArray[$_];
    chomp($temp);
    $BondTypeHash{"$temp"}= "$BondType";
    my @temp = split("-",$temp);
    chomp @temp,
    $BondTypeHash{"$temp[1]-$temp[0]"}= "$BondType";
}

open my $data2 , "< Sequence.txt";
my @temp_array = <$data2>;
close($data2);
my @input_sequence = grep (($_!~m{^\s*$|^#}), @temp_array); # remove blank lines
for (@input_sequence){$_  =~ s/^\s+|\s+$//;}
my $radInc = 2.0*3.1415926/@input_sequence ;# the radian increment to rotate (base number + 1) 
my $radius = 10.6 * @input_sequence/(2.0*3.1415926);
my @orig = ($radius,$radius,0);
my $phase = 3.1415926;
#print (@input_sequence + 1.0);
#print "radinc $radInc.\n";
#sleep(100);

my @Atom; # All of atom's information
my $BondNO1; #Bonding atom ID_1.
my $BondNO2; #Bonding atom ID_2.
my $count = 0; #Total Bases counting.
my $AtomCount = 0; #Total atoms counting.
my $Coord;  #coordinate information in anonymous array. 
my $DisOfEachBase;
my $input_sequence; #Last atoms's ID.
my @BondNO;  #Bonding number.
my $BondCount = 0; #count the total bond pair
my %ID2Type;
my @MoleculeID;
my @AtomNO;  #element id -> atom ID
my @AtomType;##element id -> atom type
my @Coordinates;
my @currentPos;#out of the bracket
my $moleculeID;
my %MoleculeHash;#atom ID => molecule ID
my @AtomDone;
my @BondDone;
my @AngleDone;
my @TorsionDone;
my @GetAtomType;
my @GetBondType;
my @GetBendType;
my @GetTorsionType;
my @Bond;
my $Angleinput;
my $bonding;
my @topology;# bonded atom number, bonded IDs
require "./Bending.pl";
require "./Torsion.pl";

for my $baseSeq (0..$#input_sequence-1)
    {
      my $AngRot = $count * $radInc;
      my $currentRad = $phase + $AngRot;
      #@currentPos = map {$_ * $count} @growthV;
      $currentPos[0] = $orig[0] + $radius * cos($currentRad); 
      $currentPos[1] = $orig[1] + $radius * sin($currentRad); 
      $currentPos[2] = 0.0; 
      $count++;
     $moleculeID = $count;
     my $temp = $input_sequence[$baseSeq];
# bonding information 
  for my $bond(0..$#{$blocks{$temp}->[1]})  
            {
               $BondCount++;
               if($blocks{$temp}->[1]->[$bond] =~ m/(\d+)\-(\d+)/x)
                    {
                        $BondNO1 = $1+$AtomCount;
                        $BondNO2 = $2+$AtomCount;
                        push  @BondNO,[$BondCount,$BondNO1,$BondNO2];
                    } 
            }
############# Main atom information of BASES 
        for my $AtomNO (0..$#{$blocks{$temp}->[0]}) 
            {
                $AtomCount++;
                $Coord = $AtomNO+2;#according to format
                push @AtomNO,"$AtomCount";
                push @MoleculeID,"$moleculeID";
                push @AtomType,"$blocks{$input_sequence[$baseSeq]}->[0]->[$AtomNO]";
                my @orginalCoord = @{$blocks{$temp}->[$Coord]};
                my @temp;
                $temp[0] = $orginalCoord[0]*cos($AngRot) - $orginalCoord[1]*sin($AngRot);
                $temp[1] = $orginalCoord[0]*sin($AngRot) + $orginalCoord[1]*cos($AngRot);
                $temp[2] = 0.0;
                my @sum = pairwise { $a + $b } @currentPos, @temp;
                push @Coordinates,[@sum];
            }
#linker atom information
############# P (Link) has been shifted by coordinates
      $currentRad = $currentRad + 0.5*$radInc;#based on previous one
      #@currentPos = map {$_ * $count} @growthV;
      $currentPos[0] = $orig[0] + ($radius + 2.422) * cos($currentRad); 
      $currentPos[1] = $orig[1] + ($radius + 2.422) * sin($currentRad); 
      $currentPos[2] = 0.0;
        
        $AtomCount++;
        push @AtomNO,"$AtomCount";
        push @MoleculeID,"$moleculeID";
        push @AtomType,"$blocks{Linker}->[0]->[0]";
        my @sum = pairwise { $a + $b } @currentPos, @{$blocks{Linker}->[1]};
        push @Coordinates,[@sum];
# linker bonding information 
        $BondCount++;
        $BondNO1 = $AtomCount - @{$blocks{$temp}->[0]};
        $BondNO2 = $AtomCount;
        push  @BondNO,[$BondCount,$BondNO1,$BondNO2];
# linker bonding information for next base
        $BondCount++;
        $BondNO1 = $AtomCount;
        $BondNO2 = $AtomCount + 1;
        push  @BondNO,[$BondCount,$BondNO1,$BondNO2];        
    }

############ Last Base
# bonding information for the last base
my $temp = $input_sequence[-1];
    for my $bond(0..$#{$blocks{$temp}->[1]})  
            {
               $BondCount++;
               if($blocks{$temp}->[1]->[$bond] =~ m/(\d+)\-(\d+)/x)
                    {
                        $BondNO1 = $1+$AtomCount;
                        $BondNO2 = $2+$AtomCount;
                        push  @BondNO,[$BondCount,$BondNO1,$BondNO2];
                    } 
            }

      my $AngRot = $count * $radInc;
      my $currentRad = $phase + $AngRot;
      #@currentPos = map {$_ * $count} @growthV;
      $currentPos[0] = $orig[0] + $radius * cos($currentRad); 
      $currentPos[1] = $orig[1] + $radius * sin($currentRad); 
      $currentPos[2] = 0.0; 
 
#@currentPos = map {$_ * $moleculeID} @growthV; 

$moleculeID++;#the same value of count++
for my $AtomNO(0..$#{$blocks{$temp}->[0]}) 
    {
        $AtomCount++;
                $Coord = $AtomNO+2;#according to format
                push @AtomNO,"$AtomCount";
                push @MoleculeID,"$moleculeID";
                push @AtomType,"$blocks{$input_sequence[-1]}->[0]->[$AtomNO]";
                my @orginalCoord = @{$blocks{$temp}->[$Coord]};
                my @temp;
                $temp[0] = $orginalCoord[0]*cos($AngRot) - $orginalCoord[1]*sin($AngRot);
                $temp[1] = $orginalCoord[0]*sin($AngRot) + $orginalCoord[1]*cos($AngRot);
                $temp[2] = 0.0;
                my @sum = pairwise { $a + $b } @currentPos, @temp;
                push @Coordinates,[@sum];
    }

#print Dumper (\@BondNO); #!!!!!Check Point!!!!!
for (0..$#AtomNO)
    {
        push @AtomDone,"\n$AtomNO[$_] $MoleculeID[$_] $AtomType[$_] $Coordinates[$_]->[0] $Coordinates[$_]->[1] $Coordinates[$_]->[2]";
        push @GetAtomType,"$AtomType[$_]";#?????
    }

############ Assemble 
for (0..$#AtomNO)
    {
                $ID2Type{$AtomNO[$_]} = "$AtomType[$_]";#for lammps ID
                $MoleculeHash{$AtomNO[$_]} = $MoleculeID[$_];#For identifying that the atom comes from another neighbor residue or not.
                #print "$AtomNO[$_] $MoleculeID[$_] $AtomType[$_] $Coordinates[$_]\n";
    }
################ Atom Information Done 
my @BondType;# the same data with @BondDone
my @TotalBondType;
foreach my $bondrf (@BondNO) #bond id and pairing id 
    { 
        #print "$bondrf->[1],$bondrf->[2]\n";
        my $BondType = $BondTypeHash{"$ID2Type{$bondrf->[1]}-$ID2Type{$bondrf->[2]}"};
        push @GetBondType,"$BondType";
        push @BondDone,"\n$bondrf->[0] $BondType $bondrf->[1] $bondrf->[2]";
        push @BondType,[$bondrf->[0],$BondType,$bondrf->[1],$bondrf->[2]];
        push @TotalBondType,"$BondType"; #For Get Total Bond Types
    }
#print Dumper (\@BondType); #!!!!!Check Point!!!!!
################## Bonding Done
## topology analysis

for (0..$#BondType) 
    {
        #print "$BondType[$_][0] $BondType[$_][1] $BondType[$_][2]\n";
    }


for my $bID (0..$#BondType)
    {
                my $atomi = $BondType[$bID][2] - 1;#Perl array ID from 0
                my $atomj = $BondType[$bID][3] - 1;
                #print "$atomi $atomj\n";
                $topology[$atomi][0]++;
                push @{$topology[$atomi]},$atomj;
                $topology[$atomj][0]++;
                push @{$topology[$atomj]},$atomi;                
    }
    #@topology is one more than the atom number!!
    #!! from atom 1 to atom n, $topology[0] is not used
#print Dumper (\@topology);
#for my $temp (0..$#topology){
    #print "temp: $temp\n";
    #print "$topology[$temp]->[0] $topology[$temp]->[1] $topology[$temp]->[2]\n";
#print "@{$topology[$temp]}\n";
#}

#
#print "sleeping\n";
#sleep(100);
#print @Bonding;

################# Bending Call by value
my ($AngleDone,$GetBendType) = &Bending(\@topology,\@AtomType,\@MoleculeID,\%ID2Type,\%MoleculeHash);
my @CatchAngleDone = @$AngleDone;
$GetBendType = $$GetBendType;


my @BendType;
for (0..$#CatchAngleDone)
    {
        push @AngleDone,"\n$CatchAngleDone[$_]->[0] $CatchAngleDone[$_]->[1] $CatchAngleDone[$_]->[2] $CatchAngleDone[$_]->[3] $CatchAngleDone[$_]->[4]";
        push  @BendType,[$CatchAngleDone[$_]->[2],$CatchAngleDone[$_]->[3],$CatchAngleDone[$_]->[4]]; # Bending bead ID
       # print "$CatchAngleDone[$_]->[0],$CatchAngleDone[$_]->[1],$CatchAngleDone[$_]->[2],$CatchAngleDone[$_]->[3],$CatchAngleDone[$_]->[4]\n";
    }




################# Torsion Call by value
my ($TorsionDone,$GetTorsionType) =  & Torsion(\@topology,\@AtomType,\@MoleculeID,\%ID2Type,\%MoleculeHash);
my @CatchTorsionDone = @$TorsionDone;
@GetTorsionType = @$GetTorsionType;
for (0..$#CatchTorsionDone)
    {
        push @TorsionDone,"\n$CatchTorsionDone[$_]->[0] $CatchTorsionDone[$_]->[1] $CatchTorsionDone[$_]->[2] $CatchTorsionDone[$_]->[3] $CatchTorsionDone[$_]->[4] $CatchTorsionDone[$_]->[5]";
        #say "$CatchTorsionDone[$_]->[0] $CatchTorsionDone[$_]->[1] $CatchTorsionDone[$_]->[2] $CatchTorsionDone[$_]->[3] $CatchTorsionDone[$_]->[4] $CatchTorsionDone[$_]->[5]";
    }




#@GetTorsionType = @$GetTorsionType;
#@TorsionDone = @$TorsionDone;

############### Output Structure
& OutputSturcture;


sub OutputSturcture
    {
my $atoms = $#AtomDone+1;
my $bonds = $#BondDone+1;
my $bends = $#AngleDone+1;
my $torsions = $#TorsionDone+1;
my $atomtypes = max @GetAtomType;
my $bondtypes = max @GetBondType;
my $bendtypes = $GetBendType;
my $torsiontypes = max @GetTorsionType;

print $StructureOut "LAPPMS data file\n\n$atoms atoms\n$atomtypes atom types\n$bonds bonds\n$bondtypes bond types \n$bends angles\n$bendtypes angle types\n$torsions dihedrals\n$torsiontypes dihedral types\n\nAtoms\n";

print $StructureOut "@AtomDone\n\nBonds\n@BondDone\n\nAngles\n@AngleDone\n\nDihedrals\n@TorsionDone";
    }






