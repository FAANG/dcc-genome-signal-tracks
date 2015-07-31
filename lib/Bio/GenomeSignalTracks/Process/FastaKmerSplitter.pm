package Bio::GenomeSignalTracks::Process::FastaKmerSplitter;
use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');

use autodie;
use PerlIO::gzip;
use Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit;
use Bio::GenomeSignalTracks::Util::FastaKmerWriter;

sub fetch_input {
    my ($self) = @_;

    my $kmer_size   = $self->param_required('kmer_size');
    my $fasta_file  = $self->param_required('fasta_file');
    my $output_dir  = $self->param_required('output_dir');
    my $split_limit = $self->param_required('split_limit');
    my $output_branch_id = $self->param_required('output_branch_id');
    my $gzip_output = $self->param('gzip');
}

sub write_output {
    my ($self) = @_;

    my $kmer_size   = $self->param_required('kmer_size');
    my $fasta_file  = $self->param_required('fasta_file');
    my $output_dir  = $self->param_required('output_dir');
    my $split_limit = $self->param_required('split_limit');
    my $output_branch_id = $self->param_required('output_branch_id');
    my $gzip_output = $self->param('gzip');

    my $in_fh;

    $self->dbc and $self->dbc->disconnect_when_inactive(1);

    if ( $fasta_file =~ m/\.gz$/ ) {
        open( $in_fh, '<:gzip', $fasta_file );
    }
    else {
        open( $in_fh, '<', $fasta_file );
    }

    tie( *OUT_FH, 'Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit',
        $output_dir, $split_limit, $gzip_output, '.kmers',
        "${kmer_size}_XXXXX" );

    my $current_fn;

    # Register code to listen to file creation
    ( tied *OUT_FH )->add_file_creation_listeners(
        sub {
            my ( $tied_object, $filename ) = @_;
            if ($current_fn) {
                $self->dataflow_output_id( { kmer_file => $current_fn }, $output_branch_id );
            }
            $current_fn = $filename;
        }
    );

    my $splitter = Bio::GenomeSignalTracks::Util::FastaKmerWriter->new(
        in_fh     => $in_fh,
        kmer_size => $kmer_size,
        out_fh    => \*OUT_FH,
    );

    $splitter->split();

    close($in_fh);
    close(*OUT_FH);
    
    $self->dbc and $self->dbc->disconnect_when_inactive(0);
    if ($current_fn && *OUT_FH->{num_writes}) {
      $self->dataflow_output_id( { kmer_file => $current_fn }, $output_branch_id ) 
    }

}

1;