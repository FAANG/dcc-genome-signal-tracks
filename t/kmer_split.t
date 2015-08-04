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
my $test_out_dir = tempdir( CLEANUP => 1 );
my $seq_2_start_pos = 18;
$seq_2_start_pos--;

make_path($test_out_dir);

simple_split();
multi_fasta_split();
multi_fasta_multiplex_out();
multi_fasta_seek_out();
multi_fasta_seek_limit_seqs_out();

done_testing();

sub simple_split {
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

sub multi_fasta_split {
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

    my @expected_kmers = qw(12345
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
    is_deeply( \@kmers, \@expected_kmers,
        "Splitter Kmers from multi fasta match expectation" );

    unlink($test_out_fn);
}

sub multi_fasta_multiplex_out {

    open( my $in_fh, '<', "$test_data_dir/three_seqs.fa" );

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

    my @expected_kmers = qw(
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

    is_deeply( \@kmers, \@expected_kmers,
        "Splitter Kmers from multi fasta match expectation" );
}

sub multi_fasta_seek_out {
    my $test_out_fn = "$test_out_dir/out4.kmers";
    open( my $in_fh,  '<', "$test_data_dir/three_seqs.fa" );
    open( my $out_fh, '>', $test_out_fn );

    seek( $in_fh, $seq_2_start_pos, 0 );    # seek to start of second seq

    my $splitter = Bio::GenomeSignalTracks::Util::FastaKmerWriter->new(
        in_fh     => $in_fh,
        kmer_size => 5,
        out_fh    => $out_fh

    );

    $splitter->split();

    close($in_fh);
    close($out_fh);

    my @kmers = slurp_kmers($test_out_fn);

    my @expected_kmers = qw(
      ABCDE
      FFFFF
      FFFFF
      FFFFF
      FFFFF
      FFFFF
      FFFFF
    );
    is_deeply( \@kmers, \@expected_kmers,
        "Splitter Kmers after seek in multi fasta match expectation" );

    unlink($test_out_fn);
}

sub multi_fasta_seek_limit_seqs_out {
    my $test_out_fn = "$test_out_dir/out4.kmers";
    open( my $in_fh,  '<', "$test_data_dir/three_seqs.fa" );
    open( my $out_fh, '>', $test_out_fn );

    seek( $in_fh, $seq_2_start_pos, 0 );    # seek to start of second seq

    my $splitter = Bio::GenomeSignalTracks::Util::FastaKmerWriter->new(
        in_fh        => $in_fh,
        kmer_size    => 5,
        out_fh       => $out_fh,
        max_num_seqs => 1,
    );

    $splitter->split();

    close($in_fh);
    close($out_fh);

    my @kmers = slurp_kmers($test_out_fn);

    my @expected_kmers = qw(
      ABCDE
    );
    is_deeply( \@kmers, \@expected_kmers,
        "Splitter Kmers after seek with max num of seqs to read in multi fasta match expectation" );

    unlink($test_out_fn);
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
