

#!/usr/bin/perl
# this is to generate xw.al from xw.all.join.cull
my @a;


while(<>){
		chomp;
		if (/QNAME/){
				print "$_\t1\n";
		}else{
				print "$_\t2\n";
		}
}		
