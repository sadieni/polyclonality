#!/usr/bin/perl -w
# Simulation of random colour marking of crypts and then patch counting
# Multiple runs and patch size collation

#open(OUTPUT, ">$ARGV[0]");

print "Marked\tSingle\tTwo\tThree\tFour\tFive\tSix\tSeven\n";
$total_crypts = 10000;

$blue = 280;

$iteration = 1;

while ($iteration <= 1000){
# Create the unmarked crypt field
@CRYPTS = ();
$crypts = 1;
while ($crypts <= $total_crypts){
	push(@CRYPTS,"•");
	$crypts = $crypts +1;
}
	
# Randomly mark some of these crypts 

$count = 1;
while ($count <= $blue){
	$marked_crypt = int(rand($total_crypts));
	if 	($CRYPTS[$marked_crypt] eq "•"){
		$CRYPTS[$marked_crypt] = "B";
		$count = $count +1;
	}
}

#print "\n\n";
#foreach $val (@CRYPTS){
#	print "$val";
#}

#print "\n\n";


# Count the overall marked crypt frequency and store all marked crypt IDs in array (colour_position)
$index = 0;
@MARKED_CRYPTS =();
$blue = 0;
$yellow = 0;
$red = 0;
while($index < scalar@CRYPTS){

	$position = $index +1;
	
	if ($CRYPTS[$index] eq "B"){$blue = $blue +1; push(@MARKED_CRYPTS,"B_$position")}
	if ($CRYPTS[$index] eq "Y"){$yellow = $yellow +1; push(@MARKED_CRYPTS,"Y_$position")}
	if ($CRYPTS[$index] eq "R"){$red = $red +1; push(@MARKED_CRYPTS,"R_$position")}

	$index = $index +1;
}

$all_marked = $blue + $yellow + $red;
#print "Overall marked crypt count = $all_marked\tBlue = $blue\tYellow = $yellow\tRed = $red\n\n";


#foreach $val (@MARKED_CRYPTS){
#	print "$val\t";
#}

# Count up patches of marked crypts
# Currently assuming 100 crypts per row

@PATCH=();

$single = 0;
$two = 0;
$three = 0;
$four = 0;
$five = 0;
$six = 0;
$seven = 0;

#print "\n\n";

while(scalar@MARKED_CRYPTS > 0){


	$main_index = 0;
	while($main_index < scalar@MARKED_CRYPTS){

		push(@PATCH,$MARKED_CRYPTS[$main_index]);
	
		
		splice(@MARKED_CRYPTS, $main_index, 1);
		$main_index = $main_index -1;
	
	
		foreach $val1 (@PATCH){
			($colour1,$position1) = split(/_/,$val1);
			foreach $val2 (@MARKED_CRYPTS){
				($colour2,$position2) = split(/_/,$val2);
				if(($position2 == $position1 +1) || ($position2 == $position1 +99) || ($position2 == $position1 +100) ||
					($position2 == $position1 -1) || ($position2 == $position1 -99) || ($position2 == $position1 -100)){
					push(@PATCH,$val2);
				
				
					$index = 0;
					until ($MARKED_CRYPTS[$index] eq "$val2"){$index++};
					splice(@MARKED_CRYPTS, $index, 1);
				
				}
			
			}
		
		}
#		print "Patch @PATCH\n";
		$patch_size = scalar@PATCH;
#


#		print OUTPUT "$patch_size\n";
		
		if($patch_size == 1){$single++};
		if($patch_size == 2){$two++};
		if($patch_size == 3){$three++};
		if($patch_size == 4){$four++};
		if($patch_size == 5){$five++};
		if($patch_size == 6){$six++};
		if($patch_size == 7){$seven++};

		@PATCH=();
		$main_index = $main_index +1;
	}
}

print "$blue\t$single\t$two\t$three\t$four\t$five\t$six\t$seven\n";
$iteration++;
}
exit;