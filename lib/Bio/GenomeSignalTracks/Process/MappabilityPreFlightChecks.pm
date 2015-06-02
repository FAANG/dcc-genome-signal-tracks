package Bio::GenomeSignalTracks::Process::MappabilityPreFlightChecks;
use strict;
use warnings;

use autodie;
use base ('Bio::EnsEMBL::Hive::Process');

use Scalar::Util::Numeric qw(isint);
use Try::Tiny;

sub fetch_input {
    my ($self) = @_;

    # do the executables exist, and are they actually executable?
    my %executables = (
        samtools         => $self->param_required('samtools'),
        bowtie           => $self->param_required('bowtie'),
        bedtools         => $self->param_required('bedtools'),
        bedGraphToBigWig => $self->param_required('bedGraphToBigWig'),
    );
    for my $executable ( values %executables ) {
        $self->check_executable($executable);
    }

    # does the fasta dir exist and contain fasta files?
    $self->check_fasta(
        $self->param_required('fasta_dir'),
        $self->param_required('fasta_suffix')
    );

    # does the index exist?
    $self->check_index(
        $self->param_required('index_dir'),
        $self->param_required('index_name'),
    );

    # does the output dir exist
    $self->check_output_dir( $self->param_required('output_dir') );

    # does the chrom_list exist and contain two columns, second one numeric?
    $self->check_chrom_list( $self->param_required('chrom_list') );

    # are the kmer sizes positive integers
    $self->check_positive_integer_list( $self->param_required('kmer_sizes'),
        'kmer_sizes' );
}

sub check_positive_integer_list {
    my ( $self, $list_var, $param_name ) = @_;

    try {
        my $int_list = [ eval($list_var) ];
        for my $i (@$int_list) {
            if ( !isint($i) == 1 ) {
                die("$i isn't a positive integer");
            }
        }
    }
    catch {
        $self->throw("Error when converting $param_name to list of positive integers: $_");
    };

}

sub check_chrom_list {
    my ( $self, $chrom_list ) = @_;

    if (! -e $chrom_list ) {
        $self->throw("Chrom list does not exist: $chrom_list");
    }

    open( my $fh, '<', $chrom_list );

    my $counter;
    while (<$fh>) {
        chomp;
        my @line = split /\t/;

        if ( scalar(@line) == 1 ) {
            $self->throw("Chrom list line $. has too few columns");
        }
        if ( scalar(@line) > 1 ) {
            $counter++;
            my $label  = $line[0];
            my $length = $line[1];

            if ( isint($length) != 1 ) {
                $self->throw(
"Chrom list. Entry for $label does not have a valid length: $length"
                );
            }
        }
    }

    if ( $counter < 1 ) {
        $self->throw("Chrom list does not contain any lines");
    }

    close($fh);
}

sub check_output_dir {
    my ( $self, $dir ) = @_;

    if ( !-d $dir ) {
        $self->throw("Output dir is not a directory: $dir");
    }

    if ( !-w $dir ) {
        $self->throw("Output dir is not writeable: $dir");
    }
}

sub check_index {
    my ( $self, $index_dir, $index_name ) = @_;

    if ( !-d $index_dir ) {
        $self->throw("Index dir is not a directory: $index_dir");
    }

    my $index_pattern = "$index_dir/$index_name.*.ebwt";

    my @index_files = glob($index_pattern);

    if ( scalar(@index_files) < 1 ) {
        $self->throw(
            "Could not find any index files matching pattern $index_pattern");
    }
}

sub check_executable {
    my ( $self, $executable ) = @_;

    if ( !-x $executable ) {
        $self->throw("$executable is not executable");
    }
}

sub check_fasta {
    my ( $self, $dir, $suffix ) = @_;

    my @files = glob("$dir/*.$suffix");

    if ( scalar(@files) < 1 ) {
        $self->throw(
            "Could not find any files in directory $dir matching *.$suffix");
    }
}

1;
