454_tag_sorting
===============

Script for finding paired tags ("Binladen tags") in 454 data

NOTE: Work in progress!
-----------------------


Documentation for find_tag_pairs.pl
-----------------------------------


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



