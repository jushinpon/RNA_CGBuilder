#use List::MoreUtils qw(mesh); 
use 5.010;
use strict; 
use warnings;
#use diagnostics;
use List::Util qw(min max sum);
use List::MoreUtils 'pairwise';
use Data::Dumper;
#use List::MoreUtils qw(uniq);

require "./Bending.pl";

open my $StructureOut , "> ps_b_p2vp.data";#the topology file 
my $beadtypeNo = 3;#the total bead types
my @AngleTypeArray;#all angle types
my %AngleType2ID;#angle type (i-j-k) to angle type ID (1,2,3..)
my $Ang_counter = 0;
for my $i (1..$beadtypeNo){
    for my $j (1..$beadtypeNo){
        for my $k (1..$beadtypeNo){
            if($i <= $j){
                $Ang_counter++;
                my $temp = "$i-$j-$k";
                push @AngleTypeArray,$temp;
                $AngleType2ID{"$i-$j-$k"} = $Ang_counter;
            }
        }
    }
}
#print Dumper (\@AngleTypeArray); #!!!!!Check Point!!!!!
#die;
#RNA grows along y dimension
my @growthV = (0,2,0); # the vector to grow your chain,could be variable
#current types-> SCY:1, STY:2, P4:3 1-1 1-2 1-3 2-2 2-3 3-3
#{local typeArray,local BondingArray(topology by atom id),CoordinationArrays according to bead number in typeArray}
my @block_A = ([2,2,2],["1-2","2-3","1-3"],[0,0,0],[0,-1,-1],[0,1,-1]);#PS: STY only
my @block_B = ([3,2,2],["1-2","2-3","1-3"],[0,0,0],[0,-1,-1],[0,1,-1]);#P2VP: P4 and two STY
#my @linker = ([1],[-2.422,-5.3,0]); #P
my @linkerA = ([1],[0,0,0]); #SCY
my @linkerB = ([3],[0,0,0]); #P4
my %type2Charge = ( #hash to convert type to its charge
   1 => 0.0,
   2 => 0.0,
   3 => 0.0
);
#the first one the is monomer number
my @moInfo = (
    [248,\@block_A,\@linkerA],
    [195,\@block_B,\@linkerB]
    );# monomer No for each type, and related informaton

#Define bond types in lammps (make inverse also ok)
my @BondTypeArray = qw/1-1 1-2 1-3 2-2 2-3 3-3/;#smaller type ID first for sorting (id to type),id starts from 0
my %BondTypeHash;#bond type to lmp bond ID
my $BondType = 0;
for (0..$#BondTypeArray){
    chomp;
    $BondType++;#lmp from 1
    chomp($BondType);
    my $temp = $BondTypeArray[$_];
    chomp($temp);
    $BondTypeHash{"$temp"}= "$BondType";#convert pair (bond type) to lmp bond ID
 }

my @Atom; # All of atom's information
my $BondNO1; #Bonding atom ID_1.
my $BondNO2; #Bonding atom ID_2.
my $monocount = 0; #Total grown monomer
my $AtomCount = 0; #Total atoms counting and ID.
my $BondCount = 0; #count the total bond pair
#[$BondCount,$BondNO1,$BondNO2]
my @AllBondInfo;  #Bonding information for all bonds [bid,atom1,atom2,lmp type].
my $Coord;  #coordinate information in anonymous array. 

my $DisOfEachBase;
my $input_sequence; #Last atoms's ID.
my %ID2Type;
my @MoleculeID;# for different monomer types
my @Atomid;  #element id -> atom ID
my @AtomType;##element id -> atom type
my @Coordinates;
my @currentPos = (0,0,0);#current grown position (growthV unit)
my $moleculeID = 0;# for different monomer types
my %MoleculeHash;#atom ID => molecule ID
my @Atom4data;#format for the data file
my @Bond4data;#format for the data file
my @Angle4data;# information of all angles in lmp data format (id type i j k)

my @GetAtomType;
my @BondType4bid;# bond types for all bond ids
my @GetBendType;
my @Bond;
my $Angleinput;
my $bonding;
# @topology: 2-d array used for finding all bending or torsion network
#[atomid -1][0]: bonded atom number,[atomid -1][1..]: the bonded atom id -1 
my @topology;# bonded atom number, bonded IDs
#require "./Torsion.pl";#no torsion term for MARTIN applied to polymer

for my $moIn (@moInfo){#loop for assembling monomers. Loop over different monomer types
    $moleculeID++;
    my $mNo = $moIn ->[0];#monomer No of each type
    #[2,2,2]
    my $atomNo4mo = @{$moIn ->[1]->[0]};#atom No of each monomer block
    print "\$mNo: $mNo\n";
    #["1-2","2-3","1-3"]
    my @lbonds = @{$moIn ->[1]->[1]};#local bond by atomid (not type id)
#    print "\@lbonds: @lbonds\n";
#die;
#$AtomCount = 0 and @currentPos == 0 for the first loop
#beginning to assemble the segment of a monomer type
    for my $mono (1.. $mNo){#loop over the monomer sequence of a monmer type
       $monocount++;
        for my $b (0 .. $#lbonds){#loop over local monomer bond ID
            $BondCount++;
            my @temp = split("-",$lbonds[$b]);
            $BondNO1 = $temp[0] + $AtomCount;#convert local atom id to global atom id
            $BondNO2 = $temp[1] + $AtomCount;#convert local atom id to global atom id
            #the following is for lmp data file for bonds iterm        
            push  @AllBondInfo,[$BondCount,$BondNO1,$BondNO2];#for data file
        }
    ############### setting atom information 
        for my $id (0.. $atomNo4mo -1){#loop over block local atom id
             $AtomCount++;
             $Coord = $id + 2;#according to block array format, @block_A
             push @Atomid,"$AtomCount";#????????
             push @MoleculeID,"$moleculeID";#for different monomer types
             push @AtomType,"$moIn->[1]->[0]->[$id]";#"$blocks{$input_sequence[$baseSeq]}->[0]->[$AtomNO]";
             #print "id, Coord: $id, $Coord\n";
             #print "$moIn->[1]->[$Coord]\n";
            
             my @orginalCoord = @{$moIn->[1]->[$Coord]};

             my @sum = pairwise { $a + $b } @currentPos, @orginalCoord;
             push @Coordinates,[@sum];
        }     
    #        linker atom information
    #############  (Link) has been shifted by coordinates
          @currentPos = map {$_ * $monocount} @growthV;

            $AtomCount++;
            push @Atomid,"$AtomCount";
            push @MoleculeID,"$moleculeID";
            push @AtomType,"$moIn->[2]->[0]->[0]";
            my @sum = pairwise { $a + $b } @currentPos,  @{$moIn ->[2]->[1]};
            push @Coordinates,[@sum];
    # linker bonding information 
            $BondCount++;
            $BondNO1 = $AtomCount - @{$moIn ->[1]->[0]};#bonded with the first atom in @block_X
            $BondNO2 = $AtomCount;
            push  @AllBondInfo,[$BondCount,$BondNO1,$BondNO2];
    # linker bonding information for next base
            $BondCount++;
            $BondNO1 = $AtomCount;
            $BondNO2 = $AtomCount + 1;#bonded by the first atom of next monomer
            push  @AllBondInfo,[$BondCount,$BondNO1,$BondNO2];        

    }
}      

# 
pop @AllBondInfo;#take the last bonding information out (after all loops), no further new monomer

for (0..$#Atomid){
    my $temp = $AtomType[$_];
    my $charge = $type2Charge{"$temp"};
    chomp ($Atomid[$_],$MoleculeID[$_],$AtomType[$_],$charge,
     $Coordinates[$_]->[0],$Coordinates[$_]->[1],$Coordinates[$_]->[2]);
    push @Atom4data,[$Atomid[$_],$MoleculeID[$_],$AtomType[$_],$charge,
     $Coordinates[$_]->[0],$Coordinates[$_]->[1],$Coordinates[$_]->[2]];
    #push @GetAtomType,"$AtomType[$_]";#?????
}
#print Dumper (\@Atom4data); #!!!!!Check Point!!!!!

#
############# Assemble 
for (0..$#Atomid){
    $ID2Type{$Atomid[$_]} = "$AtomType[$_]";#for lammps ID
    #$MoleculeHash{$AtomNO[$_]} = $MoleculeID[$_];#For identifying that the atom comes from another neighbor residue or not.
    #print "$AtomNO[$_] $MoleculeID[$_] $AtomType[$_] $Coordinates[$_]\n";
}
################# Atom Information Done 
## get bond types for the lmp data file
my @BondType;# the same data with @BondDone
my @TotalBondType;
#[$BondCount,$BondNO1,$BondNO2]
foreach my $bondrf (@AllBondInfo){ #bond id and pairing id 
     
    chomp ($ID2Type{$bondrf->[1]},$ID2Type{$bondrf->[2]});
    my @temp = ($ID2Type{$bondrf->[1]},$ID2Type{$bondrf->[2]});
    my @temp1 = sort {$b < $a} @temp;
    my $pair = join("-", @temp1);
    #print "\$pair $pair\n";
    #BondTypeHash:pair to lmp bond ID
    my $BondType = $BondTypeHash{"$pair"};
    
    #print "$bondrf->[1],$bondrf->[2]\n";
    #print "$ID2Type{$bondrf->[1]},$ID2Type{$bondrf->[2]}\n";
    #print "\$BondType: $BondType\n";
    push @BondType4bid,"$BondType";#??????
    chomp ($bondrf->[0],$BondType,$bondrf->[1],$bondrf->[2]);
    push @Bond4data,[$bondrf->[0],$BondType,$bondrf->[1],$bondrf->[2]];
#    push @BondType,[$bondrf->[0],$BondType,$bondrf->[1],$bondrf->[2]];
    push @TotalBondType,"$BondType"; #For Get Total Bond Types?????
}

#print Dumper (\@Bond4data); #!!!!!Check Point!!!!!

for my $bID (0..$#Bond4data){
    my $atomi = $Bond4data[$bID][2] - 1;#Perl array ID from 0
    my $atomj = $Bond4data[$bID][3] - 1;
    #print "$atomi $atomj\n";
    $topology[$atomi][0]++;#bonded atom counter for atomi
    push @{$topology[$atomi]},$atomj;
    $topology[$atomj][0]++;#bonded atom counter for atomj
    push @{$topology[$atomj]},$atomi;                
}
    #@topology is one more than the atom number!!
    #!! from atom 1 to atom n, $topology[0] is not used
#print Dumper (\@topology);
#
################## Bending Call by value
my $Angle4data = &Bending(\@AngleTypeArray,\%AngleType2ID,\@topology,\@AtomType,\@MoleculeID,\%ID2Type,\%MoleculeHash);
# a little stupid, but more clear for subroutine
@Angle4data = @$Angle4data;
#$GetBendType = $$GetBendType;
#
#
#my @BendType;
#for (0..$#CatchAngleDone)
#    {
#        push @AngleDone,"\n$CatchAngleDone[$_]->[0] $CatchAngleDone[$_]->[1] $CatchAngleDone[$_]->[2] $CatchAngleDone[$_]->[3] $CatchAngleDone[$_]->[4]";
#        push  @BendType,[$CatchAngleDone[$_]->[2],$CatchAngleDone[$_]->[3],$CatchAngleDone[$_]->[4]]; # Bending bead ID
#       # print "$CatchAngleDone[$_]->[0],$CatchAngleDone[$_]->[1],$CatchAngleDone[$_]->[2],$CatchAngleDone[$_]->[3],$CatchAngleDone[$_]->[4]\n";
#    }
#
################ Output Structure
#& OutputSturcture;
#
#
#sub OutputSturcture
#    {
#my $atoms = $#AtomDone+1;
#my $bonds = $#BondDone+1;
#my $bends = $#AngleDone+1;
#my $torsions = $#TorsionDone+1;
#my $atomtypes = max @GetAtomType;
#my $bondtypes = max @GetBondType;
#my $bendtypes = $GetBendType;
#my $torsiontypes = max @GetTorsionType;
#
#print $StructureOut "LAPPMS data file\n\n$atoms atoms\n$atomtypes atom types\n$bonds bonds\n$bondtypes bond types \n$bends angles\n$bendtypes angle types\n$torsions dihedrals\n$torsiontypes dihedral types\n\nAtoms\n";
#
#print $StructureOut "@AtomDone\n\nBonds\n@BondDone\n\nAngles\n@AngleDone\n\nDihedrals\n@TorsionDone";
#    }






