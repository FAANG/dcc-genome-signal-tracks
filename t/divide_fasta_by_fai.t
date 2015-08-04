#!/usr/bin/env perl
use strict;

use FindBin qw($Bin);
use File::Temp qw/ tempdir /;
use lib "$Bin/../lib";
use autodie;
use Data::Dumper;
use Test::More;
use Bio::EnsEMBL::Hive::Utils::Test qw(standaloneJob);

my $test_data_dir = "$Bin/data";
my $test_fai      = "$test_data_dir/Galgal4.fa.fai";

standaloneJob(
    'Bio::GenomeSignalTracks::Process::DivideFastaByFai',
    {
        'fai'               => $test_fai,
        'target_base_pairs' => 500_000_000
        ,    #large limit so we don't have to many chunks to worry about

    },
    [
        [
            'DATAFLOW',
            {
                num_seqs_to_read => 13,
                seq_start_pos    => 57,
            },
            2
        ],
        [
            'DATAFLOW',
            {
                num_seqs_to_read => 18,
                seq_start_pos    => 503626595,
            },
            2
        ],
        [
            'DATAFLOW',
            {
                num_seqs_to_read => 15901,
                seq_start_pos    => 936036306,
            },
            2
        ]

    ],
);

done_testing();