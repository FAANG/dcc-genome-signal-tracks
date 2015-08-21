package Bio::GenomeSignalTracks::Process::FastaKmerSplitter;
use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');

use autodie;

#use PerlIO::gzip;
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
    my $first_seq_name   = $self->param('first_seq_name');
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
    my $first_seq_name   = $self->param('first_seq_name');

    my $in_fh;

    $self->dbc and $self->dbc->disconnect_when_inactive(1);

    if ( $fasta_file =~ m/\.gz$/ ) {
        open( $in_fh, '-|', 'gzip', '-dc', $fasta_file );
    }
    else {
        open( $in_fh, '<', $fasta_file );
    }
    tie my $timer, 'Time::Stopwatch';

    if ($seq_start_pos) {
        $seq_start_pos = $seq_start_pos - 2;
        seek( $in_fh, $seq_start_pos, 0 );
    }

    my $splitter;
    my $counter;
    my $suffix = '.kmers';

    if ($gzip_output) {
        $suffix .= '.gz';
    }

    my $filenamer_callback = sub {
        my $seq_name = $splitter->current_seq_name();
        my $byte_pos = tell($in_fh);

        my $count = $counter++;

        my $file_name =
          "k${kmer_size}_${seq_name}_${count}${suffix}";
        return $file_name;
    };

    tie( *OUT_FH, 'Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit',
        $output_dir, $split_limit, $gzip_output, $filenamer_callback );

    $splitter = Bio::GenomeSignalTracks::Util::FastaKmerWriter->new(
        in_fh          => $in_fh,
        kmer_size      => $kmer_size,
        out_fh         => \*OUT_FH,
        first_seq_name => $first_seq_name,
        max_num_seqs   => $num_seqs_to_read,
    );

    $splitter->split();
    close($in_fh);
    close(*OUT_FH);

    if ($num_seqs_to_read && $num_seqs_to_read != $splitter->seqs_read()){
      my $start_pos = $seq_start_pos || 0;
      die "Expected to read $num_seqs_to_read seqs from $fasta_file (starting at $start_pos), but read ".$splitter->seqs_read();
    }

    $self->dbc and $self->dbc->disconnect_when_inactive(0);

    for ( ( tied *OUT_FH )->get_filenames ) {
        $self->dataflow_output_id( { kmer_file => $_ }, $fan_branch_code );
    }

}

1;
