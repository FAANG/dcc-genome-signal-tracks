package Bio::GenomeSignalTracks::Process::FastaKmerSplitter;
use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');

use autodie;
use IO::Uncompress::AnyUncompress qw($AnyUncompressError);
use Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit;
use Bio::GenomeSignalTracks::Util::SlurpFastaKmerWriter;
use Bio::GenomeSignalTracks::Util::FastaKmerWriter;

sub fetch_input {
    my ($self) = @_;

    my $kmer_size   = $self->param_required('kmer_size');
    my $fasta_file  = $self->param_required('fasta_file');
    my $output_dir  = $self->param_required('output_dir');
    my $split_limit = $self->param_required('split_limit');
    my $gzip_output = $self->param('gzip');
    my $slurp       = $self->param('slurp');
}

sub write_output {
    my ($self) = @_;

    my $kmer_size   = $self->param_required('kmer_size');
    my $fasta_file  = $self->param_required('fasta_file');
    my $output_dir  = $self->param_required('output_dir');
    my $split_limit = $self->param_required('split_limit');
    my $gzip_output = $self->param('gzip');
    my $slurp       = $self->param('slurp');

    my $in_fh;

    if ( $fasta_file =~ m/\.gz$/ ) {
        $in_fh = IO::Uncompress::AnyUncompress->new($fasta_file)
          or die "anyuncompress failed: $AnyUncompressError\n";
    }
    else {
        open( $in_fh, '<', $fasta_file );
    }

    tie(
        *OUT_FH,
        'Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit',
        $output_dir,
        $split_limit,
        $gzip_output,
        '.kmers',
        "split_${kmer_size}_XXXXX"
    );

    my $current_fn;

    # Register code to listen to file creation
    ( tied *OUT_FH )->add_file_creation_listeners(
        sub {
            my ( $tied_object, $filename ) = @_;
            if ($current_fn) {
                $self->dataflow_output_id( { kmer_file => $current_fn }, 1 );
            }

            $current_fn = $filename;
        }
    );

    my $splitter;

    if ($slurp) {
        $splitter = Bio::GenomeSignalTracks::Util::SlurpFastaKmerWriter->new(
            in_fh     => $in_fh,
            kmer_size => $kmer_size,
            out_fh    => \*OUT_FH,
        );
    }
    else {
        $splitter = Bio::GenomeSignalTracks::Util::FastaKmerWriter->new(
            in_fh     => $in_fh,
            kmer_size => $kmer_size,
            out_fh    => \*OUT_FH,
        );
    }

    $splitter->split();

    close($in_fh);
    close(*OUT_FH);
    $self->dataflow_output_id( { kmer_file => $current_fn }, 1 );

}

1;

