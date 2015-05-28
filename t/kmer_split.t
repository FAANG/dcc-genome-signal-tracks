#!/usr/bin/env perl
use strict;
use Test::More;
use FindBin qw($Bin);
use File::Temp qw/ tempdir /;
use lib "$Bin/../lib";

use Bio::GenomeSignalTracks::Util::FastaKmerWriter;
use Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit;
use File::Path qw(remove_tree make_path);
use autodie;

my $test_data_dir = "$Bin/data";
my $test_out_dir  = tempdir(CLEANUP => 1);

make_path($test_out_dir);
test1();
test2();
test3();
done_testing();


sub test1 {
    my $test_out_fn = "$test_out_dir/out1.kmers";
    open( my $in_fh,  '<', "$test_data_dir/small.fa" );
    open( my $out_fh, '>', $test_out_fn );

    my $splitter = Bio::GenomeSignalTracks::Util::FastaKmerWriter->new(
        in_fh     => $in_fh,
        kmer_size => 2,
        out_fh    => $out_fh
    );

    $splitter->split();

    close($in_fh);
    close($out_fh);

    my @kmers = slurp_kmers($test_out_fn);

    is_deeply(
        \@kmers,
        [ 'AB', 'BC', 'CD', 'DE' ],
        "Splitter Kmers match expectation"
    );
    unlink($test_out_fn);
}

sub test2 {
    my $test_out_fn = "$test_out_dir/out2.kmers";
    open( my $in_fh,  '<', "$test_data_dir/three_seqs.fa" );
    open( my $out_fh, '>', $test_out_fn );

    my $splitter = Bio::GenomeSignalTracks::Util::FastaKmerWriter->new(
        in_fh     => $in_fh,
        kmer_size => 5,
        out_fh    => $out_fh
    );

    $splitter->split();

    close($in_fh);
    close($out_fh);

    my @kmers = slurp_kmers($test_out_fn);

    my @expercted_kmers = qw(12345
      23456
      34567
      45678
      56789
      67890
      ABCDE
      FFFFF
      FFFFF
      FFFFF
      FFFFF
      FFFFF
      FFFFF
    );
    is_deeply( \@kmers, \@expercted_kmers,
        "Splitter Kmers from multi fasta match expectation" );

    unlink($test_out_fn);
}

sub test3 {

    open( my $in_fh, '<', "$test_data_dir/three_seqs.fa" );

    #open( my $out_fh, '>', $test_out_fn );

    tie( *OUT_FH, 'Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit',
        $test_out_dir, 5, 0, '.kmers', 'splitXXXXX' );

    my @out_files;

    # Register code to listen to file creation
    ( tied *OUT_FH )->add_file_creation_listeners(
        sub {
            my ( $tied_object, $filename ) = @_;

            push @out_files, $filename;
        }
    );

    my $splitter = Bio::GenomeSignalTracks::Util::FastaKmerWriter->new(
        in_fh     => $in_fh,
        kmer_size => 5,
        out_fh    => \*OUT_FH,
    );

    $splitter->split();

    close($in_fh);
    close(*OUT_FH);

    my @kmers;

    for my $of (@out_files) {
        push @kmers, slurp_kmers($of);
        unlink($of);
    }

    my @expercted_kmers = qw(
      12345
      23456
      34567
      45678
      56789
      67890
      ABCDE
      FFFFF
      FFFFF
      FFFFF
      FFFFF
      FFFFF
      FFFFF
    );

    is_deeply( \@kmers, \@expercted_kmers,
        "Splitter Kmers from multi fasta match expectation" );
}


sub slurp_kmers {
    my ($file) = @_;
    my @kmers;
    open( my $results_in_fh, '<', $file );
    while (<$results_in_fh>) {
        chomp;
        push @kmers, $_;
    }
    close $results_in_fh;
    return @kmers;
}
