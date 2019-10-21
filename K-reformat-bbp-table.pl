#!/usr/bin/perl
       
# Author: Ksenia Krasileva                                                                                                                                
# Contact: kseniak [at] berkeley.edu

# Feed a tabular file containg gene ids, their pathway assignments and blast hits to the script                                                          

use strict;
use warnings;
use Getopt::Long;

my $usage = "Usage: perl K-reformat-bbp-table.pl <options>\n                                                                                             
Necessary options for the script to run:                                                                                                                 
-i|--input filename                                                                                                                                      
-o|--output filename                                                                                                                                     
";
my ($infile, $outfile);

GetOptions(
       'i|input:s'    => \$infile,
       'o|output:s'    => \$outfile,
       );
die $usage unless ( defined($infile) );

my %table;

open (my $FILEHANDLE, "<", $infile) or die "cannot open input file: $infile";

my $FILEOUT;

if (defined($outfile)) {

    open ($FILEOUT, ">", $outfile);
}

while (my $line = <$FILEHANDLE>){

chomp $line;

my ($gene, $pathway, $hit) = split (/\s/, $line);

push (@{$table{$gene}{'pathways'}}, $pathway);

push (@{$table{$gene}{'hits'}}, $hit);

}

close $FILEHANDLE;

my $key;

foreach $key ( keys %table ) {

my @pathways = &uniq ( @{$table{$key}{'pathways'}} );
my @hits = &uniq ( @{$table{$key}{'hits'}} );

    if (defined($outfile)){
        print $FILEOUT $key, "\t", join(", ", @pathways), "\t", join(", ", @hits), "\n";

    }

    else {
        print $key, "\t", join(", ", @pathways), "\t", join(", ", @hits), "\n";
    }

}

sub uniq{

    my @non_unique_array = @_;
    my %seen =() ;
    my @unique_array = grep { ! $seen{$_}++ } @non_unique_array ;
    return @unique_array;
}



__END__  

