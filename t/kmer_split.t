#!/usr/bin/env perl
use strict;
use Test::More;
use FindBin qw($Bin);
use File::Temp qw/ tempdir /;
use lib "$Bin/../lib";
use PerlIO::gzip;
use Bio::GenomeSignalTracks::Util::FastaKmerWriter;
use Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit;
use File::Path qw(remove_tree make_path);
use autodie;

use Data::Dumper;

my $test_data_dir   = "$Bin/data";
my $test_out_dir    = tempdir( CLEANUP => 1 );
my $seq_2_start_pos = 18;
$seq_2_start_pos--;

my @multi_kmer_expected = qw(
  >a:1-5
  12345
  >a:2-6
  23456
  >a:3-7
  34567
  >a:4-8
  45678
  >a:5-9
  56789
  >a:6-10
  67890
  >b:1-5
  ABCDE
  >c:1-5
  FFFFF
  >c:2-6
  FFFFF
  >c:3-7
  FFFFF
  >c:4-8
  FFFFF
  >c:5-9
  FFFFF
  >c:6-10
  FFFFF
);

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
        out_fh    => $out_fh,
    );

    $splitter->split();

    close($in_fh);
    close($out_fh);

    my @kmers    = slurp_kmers($test_out_fn);
    my $expected = [
        '>seqname:1-2', 'AB', '>seqname:2-3', 'BC',
        '>seqname:3-4', 'CD', '>seqname:4-5', 'DE'
    ];

    is_deeply( \@kmers, $expected, "Splitter Kmers match expectation" );
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

    my @expected_kmers = @multi_kmer_expected;
    is_deeply( \@kmers, \@expected_kmers,
        "Splitter Kmers from multi fasta match expectation" );

    unlink($test_out_fn);
}

sub multi_fasta_multiplex_out {

    open( my $in_fh, '<', "$test_data_dir/three_seqs.fa" );

    my $file_count = 0;
    
    my $file_namer = sub {
      $file_count++;
      return "tmp_${file_count}.kmers.gz";
    };

    tie( *OUT_FH, 'Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit',
        $test_out_dir, 5, 1, $file_namer );

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

    my @expected_kmers = @multi_kmer_expected;

    is_deeply( \@kmers, \@expected_kmers,
        "Splitter Kmers from multi fasta match expectation" );
}

sub multi_fasta_seek_out {
    my $test_out_fn = "$test_out_dir/out4.kmers";
    open( my $in_fh,  '<', "$test_data_dir/three_seqs.fa" );
    open( my $out_fh, '>', $test_out_fn );

    seek( $in_fh, $seq_2_start_pos, 0 );    # seek to start of second seq

    my $splitter = Bio::GenomeSignalTracks::Util::FastaKmerWriter->new(
        in_fh          => $in_fh,
        kmer_size      => 5,
        out_fh         => $out_fh,
        first_seq_name => 'b',

    );

    $splitter->split();

    close($in_fh);
    close($out_fh);

    my @kmers = slurp_kmers($test_out_fn);

    my @expected_kmers = @multi_kmer_expected[ 12 .. $#multi_kmer_expected ];
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
        in_fh          => $in_fh,
        kmer_size      => 5,
        out_fh         => $out_fh,
        max_num_seqs   => 1,
        first_seq_name => 'b',
    );

    $splitter->split();

    close($in_fh);
    close($out_fh);

    my @kmers = slurp_kmers($test_out_fn);

    my @expected_kmers = @multi_kmer_expected[ 12 .. 13 ];
    is_deeply( \@kmers, \@expected_kmers,
"Splitter Kmers after seek with max num of seqs to read in multi fasta match expectation"
    );

    unlink($test_out_fn);
}

sub slurp_kmers {
    my ($file) = @_;
    my @kmers;
    my $op = '<';
    if ($file =~ m/\.gz$/) {
      $op .= ':gzip';
    }
    
    open( my $results_in_fh, $op, $file );
    while (<$results_in_fh>) {
        chomp;
        push @kmers, $_;
    }
    close $results_in_fh;
    return @kmers;
}
