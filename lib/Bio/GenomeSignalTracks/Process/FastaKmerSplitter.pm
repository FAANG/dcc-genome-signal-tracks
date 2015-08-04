package Bio::GenomeSignalTracks::Process::FastaKmerSplitter;
use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');

use autodie;
use PerlIO::gzip;
use Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit;
use Bio::GenomeSignalTracks::Util::FastaKmerWriter;
use Time::HiRes;
use Time::Stopwatch;

sub fetch_input {
    my ($self) = @_;

    my $kmer_size        = $self->param_required('kmer_size');
    my $fasta_file       = $self->param_required('fasta_file');
    my $output_dir       = $self->param_required('output_dir');
    my $split_limit      = $self->param_required('split_limit');
    my $fan_branch_code  = $self->param('fan_branch_code');
    my $gzip_output      = $self->param('gzip');
    my $seq_start_pos    = $self->param('seq_start_pos');
    my $num_seqs_to_read = $self->param('num_seqs_to_read');
}

sub param_defaults {
    return { fan_branch_code => 2, };
}

sub write_output {
    my ($self) = @_;

    my $kmer_size        = $self->param_required('kmer_size');
    my $fasta_file       = $self->param_required('fasta_file');
    my $output_dir       = $self->param_required('output_dir');
    my $split_limit      = $self->param_required('split_limit');
    my $fan_branch_code  = $self->param('fan_branch_code');
    my $gzip_output      = $self->param('gzip');
    my $seq_start_pos    = $self->param('seq_start_pos');
    my $num_seqs_to_read = $self->param('num_seqs_to_read');

    my $in_fh;

    $self->dbc and $self->dbc->disconnect_when_inactive(1);

    if ( $fasta_file =~ m/\.gz$/ ) {
        open( $in_fh, '<:gzip', $fasta_file );
    }
    else {
        open( $in_fh, '<', $fasta_file );
    }
    tie my $timer, 'Time::Stopwatch';

    if ($seq_start_pos) {
        my $adjusted_seek_start = $seq_start_pos - 2;
        seek( $in_fh, $adjusted_seek_start, 0 );
        print STDERR "Seek to $seq_start_pos: $timer\n";
    }

    tie( *OUT_FH, 'Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit',
        $output_dir, $split_limit, $gzip_output, '.kmers',
        "${kmer_size}_XXXXX" );

    my $current_fn;

    # Register code to listen to file creation
    ( tied *OUT_FH )->add_file_creation_listeners(
        sub {
            my ( $tied_object, $filename ) = @_;
            print STDERR "New file - $filename: $timer\n";
            if ($current_fn) {
                print STDERR "Emmitted new job for $current_fn\n";
                $self->dataflow_output_id( { kmer_file => $current_fn },
                    $fan_branch_code );
            }
            $current_fn = $filename;
        }
    );

    my $splitter = Bio::GenomeSignalTracks::Util::FastaKmerWriter->new(
        in_fh     => $in_fh,
        kmer_size => $kmer_size,
        out_fh    => \*OUT_FH,
    );

    print STDERR "Start split\n";
    $splitter->split();
    print STDERR "End split\n";
    close($in_fh);
    close(*OUT_FH);

    $self->dbc and $self->dbc->disconnect_when_inactive(0);
    if ( $current_fn  ) {
        $self->dataflow_output_id( { kmer_file => $current_fn },
            $fan_branch_code );
        print STDERR "Emmitted last job for $current_fn\n";
    }

}

1;
