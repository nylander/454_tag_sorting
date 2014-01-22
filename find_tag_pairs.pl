#!/usr/bin/perl 

=pod

=head2

            FILE:  find_tag_pairs.pl

           USAGE:  See ./find_tag_pairs.pl --help

                   ./find_tag_pairs.pl  -t=tags.tab  dat.fas
                   ./find_tag_pairs.pl  -t=tags.tab  -o=out.fas  dat.fas
                   ./find_tag_pairs.pl  -t=tags.tab  -r=rejects.fas  -o=out.fas  dat.fas
                   ./find_tag_pairs.pl  -t=tags.tab  --mask  -o=out.fas  dat.fas
                   ./find_tag_pairs.pl  -t=tags.tab  -o=out.fas --norejectfile --nosplit --nomask  dat.fas

                   ./find_tag_pairs.pl  -t=tags.tab  dat.fas  >  dat.retag.fas
                   ./find_tag_pairs.pl  -t=tags.tab  dat.fas.gz  |  gzip  >  dat.retag.fas.gz
                   ./find_tag_pairs.pl  -t=tags.tab  dat.fas.gz  |  zip  >  dat.retag.fas.zip

         OPTIONS:  -t=<file>, --tagfile=<file>    Specify <file> with tags.
                   -o=<file>, --outfile=<file>    Specify outfile (optional).
                   -r=<file>, --rejectfile=<file> Specify reject file (optional).
                   -nor, --norejectfile           Do not print reject file.
                   -m, --mask                     Mask (lower case) the matched parts of the sequence.
                   -nom, --nomask                 Do not mask. Default.
                   -s, --split                    Split found sequences into tag-specific files. Default.
                   -nos, --nosplit                Do not print separate files.
                   -h, --help                     Prints usage.
                   -v, --verbose                  Be verbose.
                   --noverbose                    Be quiet.

     DESCRIPTION:  Reads fasta file, primarily from a 454 sequencing run where samples have "paired tags" (see below), and
                   searches the sequences for tag pairs (specified in the tag-file, see below) in the end of the sequences.
                   If both tags in the pair is found, the sequence is captured, labelled and printed to STDOUT and/or to file.

                   Note: the match needs to be exact, i.e., all sequences with any mismatches will be rejected and put in
                   the rejects file.

                   Input SEQUENCE FILE:

                   Expected input sequence (Update 09/03/2013 10:03:50 AM: The lengths of the sequences are taken from the tags.tab file): 

                    MID-tag 1 (11 bp)  Tag 1 (6 bp)  Forward primer (25 bp)  The Sequence  Reverse primer (20 bp) Tag 2 (6 bp)  MID-tag 2 (11 bp)
                   |-----------------|-------------|-----------------------|--------------|----------------------|------------|-----------------|


                   Sequence length of all tags and primers: 11+6+25+X+20+6+11 = 79+X

                   The string ' mytag=<tag name>' will also be added to the fasta header. For example:
                    >H4XH63U01BXA6A rank=0176223 x=672.0 y=1344.0 length=383 mytag=MID4_L4_M6

                   Note: in this example, the length in the original label 'length=<nr>' will not be changed when trimming the output sequence. 
                   

                   Input TAGFILE:

                   The file with tags needs to be formatted this way (tab- or white-space separated columns). The header
                   ("Name_primer_tag" etc) is optional, but if present, should look like this. Pay close attention to
                   how the sequences are arranged and concatenated.

                   Name_primer_tag Name_primer_tag_mid Fwd_mid Fwd_tag Rev_mid Rev_tag Fwd_mid_Fwd_tag Rev_tag_Rev_mid
                   LCO_1_MLepR1_1 MID1_L1_M1 ACACGACGACT TACAAG AGTCGTGGTGT CTTGTA ACACGACGACTTACAAG CTTGTAAGTCGTGGTGT
                   LCO_1_MLepR1_2 MID1_L1_M2 ACACGACGACT TACAAG AGTCGTGGTGT GGTGTA ACACGACGACTTACAAG GGTGTAAGTCGTGGTGT



    REQUIREMENTS:  ---

            BUGS:  Prints zero length sequences in STDOUT stream (IE4ERIL01ATX89, IE4ERIL01BF6K8, IE4ERIL01C0YQU, IE4ERIL01CAV50, IE4ERIL01CSZOX, IE4ERIL01DNOAW)

           NOTES:  ---

          AUTHOR:  Johan A. A. Nylander (JN), <johan.nylander @ bils.se>

         COMPANY:  BILS

         VERSION:  1.0

         CREATED:  03/04/2011 07:45:16 PM CET

        REVISION:  09/03/2013 04:10:20 PM

            TODO:  Solve zero-length sequence bug in output.

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw( min );

exec("perldoc", $0) unless (@ARGV);


## Globals
my $VERBOSE          = 1;  # 0: no verbose
my $mask             = 0;  # 0: no masking of matched tag mid combinations
my $split            = 1;  # 0: no printing of separate fasta files
my %count_hash       = (); # key: name, value: n seqs
my %comb_hash        = (); # key:name, value: comb of first and second tags
my %comb_length_hash = (); # key:name, value: hash with length of first, and length of second
my %comb_print_hash  = (); # key:name, value: string with the two tags
my %split_hash       = (); # key: comb name, value: fasta entry
my %tag_order_hash        = (); # key:name, value: order in input file
my @comblengths      = ();
my $nseqs            = 0;
my $npairs           = 0;
my $nrejects         = 0;
my $ntags            = 0;
my $ntagsfound       = 0;
my $mincomblength;
my $term             = $/;
my $reportfile;
my $reportfileending = '.report.txt';
my $tagfile;
my $outfile;
my $norejectfile;
my $rejectfile;
my $rejectfileending = '.rejects.fas';


## Read args
my $r = GetOptions("tagfile=s"    => \$tagfile,
                   "verbose!"     => \$VERBOSE,
                   "outfile=s"    => \$outfile,
                   "norejectfile" => \$norejectfile,
                   "rejectfile=s" => \$rejectfile,
                   "mask!"        => \$mask,
                   "split!"       => \$split,
                   "help"         => sub { exec("perldoc", $0); exit(0); },
                  );


## Read the tagfile
## TODO: Warn if non-unique name primer tag mid (column 2)?
## TODO: Find the shortest combination of tags and mids, to be used later in reject step
READTAGFILE:
    die "need a tag file. See $0 --help\n" unless $tagfile;
    open my $TAGFILE, "<", $tagfile or die "could not open $tagfile : $! \n";
    print STDERR "Reading $tagfile\n" if ($VERBOSE);
    my $i = 0;
    while(<$TAGFILE>) {
        chomp;
        next if (/^Name/);
        my ($name_primer_tag, $name_primer_tag_mid, $Fwd_mid, $Fwd_tag, $Rev_mid, $Rev_tag, $Fwd_mid_Fwd_tag, $Rev_tag_Rev_mid) = split /\s+/, $_;
        my $name = $name_primer_tag_mid; # $name = $name_primer_tag . '-' . $name_primer_tag_mid
        my $Fwd =  uc($Fwd_mid) . uc($Fwd_tag);
        my $Rev =  uc($Rev_mid) . uc($Rev_tag);
        $comb_length_hash{$name}{'first'} = length($Fwd_mid_Fwd_tag);
        $comb_length_hash{$name}{'last'} = length($Rev_tag_Rev_mid);
        $comb_hash{$name} = uc($Fwd_mid_Fwd_tag) . uc($Rev_tag_Rev_mid);
        push @comblengths, length($comb_hash{$name});
        $comb_print_hash{$name} = "$Fwd\t$Rev";
        $tag_order_hash{$name} = $i;
        $i++;
    }

$mincomblength = min(@comblengths);


## Open fasta file (allow compression)
my $infilename = shift(@ARGV) or die "Error. No input fasta file.\n";
my $INFH;
if ( ($infilename =~ /\.zip$/) or ($infilename =~ /\.gz$/) ) {
    open $INFH, "gunzip -c $infilename |", or die " could not open $infilename for reading : $! \n";
}
else {
    open $INFH, "<", $infilename, or die " could not open $infilename for reading : $! \n";
}


## Open outfile if needed
if ($outfile) {
    if ( -e $outfile ) {
        die "$outfile already exists, will not overwrite. Quitting.\n";
    }
    else {
        open OUTFILE, ">", $outfile or die "could not open $outfile for writing : $! \n";
    }
}


## Open reject file unless norejectfile
if ($norejectfile) {
    ## Do not open rejectfile
}
elsif ($rejectfile) {
    if ( -e $rejectfile ) {
        die "$rejectfile already exists, will not overwrite. Quitting.\n";
    }
    else {
        open REJECTFILE, ">", $rejectfile or die "could not open $rejectfile for writing : $! \n";
    }
}
else {
    $rejectfile = $infilename . $rejectfileending;
    if ( -e $rejectfile ) {
        die "$rejectfile already exists, will not overwrite. Quitting.\n";
    }
    else {
        open REJECTFILE, ">", $rejectfile or die "could not open $rejectfile for writing : $! \n";
    }
}


## Read sequences from the infile
READSEQUENCEFILE:
    print STDERR "Reading $infilename (please wait...)\n" if ($VERBOSE);
    $/ = '>';
    while(<$INFH>) {
        chomp;
        next if ($_ eq '');
        my ($id, @sequencelines) = split /\n/;
        $id = '>' . $id;
        my $sequence = '';
        foreach my $line (@sequencelines) {
            $sequence .= $line;
        }
        $nseqs++;
        #next if (length($sequence) <= $mincomblength); # Saves time or not?
        $sequence = uc($sequence);

        ## Search for tags in the first N (mid+tag) and last N characters in the sequence. Accepts no mismatches!
        my $found = 0;
        foreach my $key (keys %comb_hash) {
            my $length_first = $comb_length_hash{$key}{'first'};  # Set the length here based on $comb_length_hash{$key}
            my $length_last = $comb_length_hash{$key}{'last'};    # Set the length here based on $comb_length_hash{$key}
            my $first_part = substr($sequence, 0, $length_first); # ($sequence, 0, N)
            my $last_part = substr($sequence, -$length_last);     # ($sequence, -N)
            my $first_and_last_part = $first_part . $last_part;
            
            ## Check for match: If the comination of first and last parts of the seq matches combo given in tab file
            if ($first_and_last_part eq $comb_hash{$key}) {
                $found = 1;
                if ($mask) {
                    $sequence =~ s/^($first_part)/lc($1)/e;
                    $sequence =~ s/($last_part)$/lc($1)/e;
                }
                else {
                    $sequence =~ s/^$first_part//;
                    $sequence =~ s/$last_part$//;
                }
                my $fasta = $id . " mytag=" . $key . "\n" . $sequence . "\n";
                push @{$split_hash{$key}}, $fasta; # For splitting into files, push to tag-specific arrays in hash

                ## Print sequence
                if ($outfile) {
                    print OUTFILE $id, " mytag=", $key, "\n", $sequence, "\n";
                }
                else {
                    print STDOUT $id, " mytag=", $key, "\n", $sequence, "\n";
                }
                $count_hash{$key}++;
                last;
            }
        }
        if ($found) {
            $npairs++;
            $found = 0;
        }
        else {
            $nrejects++;
            print REJECTFILE $id, "\n", $sequence, "\n" unless $norejectfile; ## Print the rejects
        }
    }
    close($INFH);
    print STDERR "Done reading input\n" if ($VERBOSE);

## Reset separator and close possible outfile
$/ = $term;
if($outfile) {
    close(OUTFILE);
}
if($rejectfile) {
    close(REJECTFILE);
}


## Create splitfiles
if ($split) {
    foreach my $key (keys %split_hash) {
        my $combofile = $key . '.fas';
        open my $CF, ">", $combofile or die "$! \n";
        print $CF @{$split_hash{$key}};
        close($CF);
    }
}


## Create report file
WRITEREPORTFILE:
    $reportfile  = $infilename . $reportfileending;
    die "$reportfile already exists. Will not overwrite. Quitting\n" if ( -e $reportfile);
    open my $REPORTFILE, ">", $reportfile or die "could not open $reportfile : $! \n";
    print $REPORTFILE "Name\tFwdMID+FwdTAG\tRevTAG+RevMID\tComb_Fwd+Rev\tcount\n";
    foreach my $key (sort keys %tag_order_hash) {
        print $REPORTFILE $key, "\t", $comb_print_hash{$key}, "\t", $comb_hash{$key}, "\t";
        if( defined $count_hash{$key} ) {
            print $REPORTFILE $count_hash{$key}, "\n";}
        else {
            print $REPORTFILE "0", "\n";
        }
    }
    close($REPORTFILE);

$ntags = keys %tag_order_hash;  # Number of tag pairs in tab file
$ntagsfound = keys %count_hash; # Number of tag-pairs found in infile

## Print some info and quit
if ( -e $reportfile ) {
    print STDERR "Created $reportfile\n" if ($VERBOSE);
}
if ($outfile) {
    if ( -e $outfile ) {
        print STDERR "Created $outfile\n" if ($VERBOSE);
    }
}
if ($rejectfile) {
    if ( -e $rejectfile ) {
        print STDERR "Created $rejectfile\n" if ($VERBOSE);
    }
}
if ($split) {
    print STDERR "Created separate fasta files for tag-pairs found.\n" if ($VERBOSE);
}
if($VERBOSE) {
    print STDERR "Summary:\n";
    print STDERR "  $ntags tag pairs to look for in\n"; 
    print STDERR "  $nseqs input sequences. And script found\n";
    print STDERR "  $npairs paired-tagged sequences coming from\n";
    print STDERR "  $ntagsfound tag pairs.\n";
    print STDERR "  $nrejects sequences were rejected.\n";
}
print STDERR "Done with all.\n" if ($VERBOSE);

exit(0);


#===  FUNCTION  ================================================================
#         NAME:  RevComp
#      VERSION:  03/04/2011 07:51:44 PM CET
#  DESCRIPTION:  Reverse complements DNA sequence. Note only ACGT are handled.
#   PARAMETERS:  string
#      RETURNS:  string
#         TODO:  ???
#===============================================================================
sub RevComp {

    my $seq = shift(@_);
    my $rseq = uc(reverse($seq));
    $rseq =~ tr/ACGT/TGCA/;

    return($rseq);

} # end of RecComp

__END__
