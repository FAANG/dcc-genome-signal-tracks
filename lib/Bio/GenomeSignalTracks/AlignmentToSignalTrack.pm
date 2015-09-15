package Bio::GenomeSignalTracks::AlignmentToSignalTrack;

use strict;
use Moose;
use Carp;
use File::Temp qw(tempdir);
use File::Copy qw(move);
use File::Basename qw(fileparse);
use IPC::System::Simple qw(system capture);
use POSIX qw(mkfifo);
use autodie;
use Parallel::ForkManager;
use Moose::Util::TypeConstraints;

subtype 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable',
  as 'Str',
  where {
    constraint => sub { -r $_ },
    message    => 'File path must exist and be readable'
  };

subtype 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Executable',
  as 'Str',
  where {
    constraint => sub { -x $_ },
    message    => 'File path must exist and be executable'
  };

#executable paths, all requried
has 'samtools_path' => (
    is => 'rw',
    isa =>
      'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Executable'
);
has 'wiggletools_path' => (
    is => 'rw',
    isa =>
      'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Executable'
);
has 'bedtools_path' => (
    is => 'rw',
    isa =>
      'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Executable'
);

# chrom sizes bed file, will be created from the bam if not otherwise present
has 'chrom_sizes_bed_file_path' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable'
);

#
has 'mappability_file_path' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable'
);
has 'mappability_auc' => ( is => 'rw', isa => 'Num' );
has 'mappability_auc_file_path' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable',
);

# required

has 'fragment_length' => ( is => 'rw', isa => 'Int' );
has 'alignment_file_path' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable'
);
has 'output_file_path' => ( is => 'rw', isa => 'Str' );

# optional
has 'read_count' => ( is => 'rw', isa => 'Int' )
  ;    #will be calculated with samtools if necessary
has 'read_length' => ( is => 'rw', isa => 'Int' )
  ;    # will calculate from the first read if not specified
has 'working_dir' => ( is => 'rw', isa => 'Maybe[Str]' )
  ;    #will use default temp folder if not specificed
has 'smooth_length' => ( is => 'rw', isa => 'Int' )
  ;    #will default to fragment length
has 'exclude_regions_file_path' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable'
);

has 'sort'         => ( is => 'rw', isa => 'Bool', default => 1 );
has 'verbose'      => ( is => 'rw', isa => 'Bool', default => 1 );
has 'cleanup_temp' => ( is => 'rw', isa => 'Bool', default => 1 );
has 'output_precision_dp' => ( is => 'rw', isa => 'Maybe[Int]', );

has 'output_type' =>
  ( is => 'rw', isa => enum( [qw[ wig bg ]] ), default => 'wig' );
no Moose::Util::TypeConstraints;

#todo - check that you have sufficient information to run the process up front
#todo - allow wig output
#todo - check that you can write to output location early
#todo - mark required information in moose declarations
#todo - filter output to set precision, tidy up any -0.000 values
#todo - can you do this in a single pass with fifos?
#todo - can you apply this to paired end data?

sub generate_track {
    my ($self) = @_;

    my $tmp_dir = $self->temp_dir;
    $self->log( 'Temp dir:', $tmp_dir );
    $self->_prep($tmp_dir);
    $self->log('Prep complete');

    #    my $fwd_output        = $tmp_dir . '/fwd.bed';
    #    my $rev_output        = $tmp_dir . '/rev.bed';

    my $fwd_scaled_output = $tmp_dir . '/fwd.scaled.shifted.bg';
    my $rev_scaled_output = $tmp_dir . '/rev.scaled.shifted.bg';

    my $fifo_mode = 0600;
    mkfifo( $fwd_scaled_output, $fifo_mode );
    mkfifo( $rev_scaled_output, $fifo_mode );

    my $combined_output = $self->_temp_output_location;

    my $pm = new Parallel::ForkManager(2);

#start the inputs into the FIFOs, they will block until somthing reads from them
    for (
        [ '-F16', $fwd_scaled_output, '+' ],
        [ '-f16', $rev_scaled_output, '-' ],
      )
    {
        my $pid = my $pid = fork; $pm->start and next;
        $self->_split_by_dir_to_bed(@$_);
        $pm->finish;    # Terminates the child process
    }

    # this now reads from the FIFOs, unblocking their inputs
    $self->log( 'Combine, extend and smooth all reads:', $combined_output );
    $self->_combine_extend_smooth( $combined_output, $fwd_scaled_output,
        $rev_scaled_output );

    $self->log( 'Move output to final location:', $self->output_file_path );
    move( $combined_output, $self->output_file_path );
    $self->log( 'Output in', $self->output_file_path );

    $self->_do_cleanup();
}

sub _temp_output_location {
    my ($self) = @_;

    my ( $name, $path ) = fileparse( $self->output_file_path );
    my $temp_output_location = $path . '.' . $name;

    $self->log( 'Using temp location for final output', $temp_output_location );

    return $temp_output_location;

}

sub _split_by_dir_to_bed {
    my ( $self, $samtools_filter, $target, $shift_op, ) = @_;

    my @cmd = (
        $self->samtools_path,
        'view -ub',
        '-F4',
        $samtools_filter,
        $self->alignment_file_path,
        '|',
        $self->bedtools_path,
        'bamtobed',
        '|',
        'awk \'BEGIN{FS="\t";OFS=FS}{print $1,$2,$3,1}\'',
    );

    if ( $self->exclude_regions_file_path ) {
        push @cmd, '|', $self->bedtools_path, 'intersect', '-sorted', '-v',
          '-a', 'stdin',
          '-b',
          $self->exclude_regions_file_path;
    }

    push @cmd, '|', 'sort -k1,1 -k2,2n' if ( $self->sort );

    push @cmd, '|', $self->wiggletools_path, 'write_bg -', 'ratio strict', '-',
      'scale', $self->_expectation_correction_factor,
      $self->mappability_file_path;

    my $shift_size = $self->_shift_size;

    push @cmd,
      '|',
      'awk',
      '-F $\'\t\'',
"'BEGIN{OFS=FS}{\$2 = \$2 $shift_op $shift_size; \$3 = \$3 $shift_op $shift_size; if(\$2<0){\$2=0};print}'";

    push @cmd, '>', $target;

    my $cmd = join( ' ', @cmd );

    $self->do_system($cmd);
}

sub do_system {
    my ( $self, $cmd ) = @_;
    $self->log( "executing", $cmd );
    system($cmd);
}

sub _combine_extend_smooth {
    my ( $self, $target, @files ) = @_;

    my $output_cmd = 'write_bg';
    if ( $self->output_type eq 'wig' ) {
        $output_cmd = 'write';
    }

    my $cmd = join( ' ',
        $self->wiggletools_path,          $output_cmd,
        '-',                              'mult strict',
        $self->chrom_sizes_bed_file_path, 'smooth',
        $self->smooth_length,             'extend',
        $self->_extend_length,            'sum',
        @files, );

    open( my $in_fh,  '-|', $cmd );
    open( my $out_fh, '>',  $target );

    while (<$in_fh>) {
        if (   m/^fixedStep/
            || m/^variableStep/ )
        {
            print $out_fh $_;
            next;
        }
        chomp;
        my @vals = split /\t/;
        if ( scalar(@vals) == 4 ) {

            #bedgraph format
            $vals[3] = $self->post_process_score( $vals[3] );
        }
        if ( scalar(@vals) == 2 ) {

            #variableStep format
            $vals[1] = $self->post_process_score( $vals[1] );
        }
        if ( scalar(@vals) == 1 ) {

            #fixedStep format
            $vals[0] = $self->post_process_score( $vals[0] );
        }

        print $out_fh join( "\t", @vals ) . "\n";
    }

    close($in_fh);
    close($out_fh);
}

sub post_process_score {
    my ( $self, $score ) = @_;

    if ( $score =~ m/^-0.00000/ ) {
        $score = 0;
    }
    if ( $self->output_precision_dp ) {
        $score = sprintf '%.' . $self->output_precision_dp . 'f', $score;
    }
    return $score;
}

sub _create_chrom_bed_file {
    my ( $self, $target ) = @_;

    my @cmd = (
        $self->samtools_path, 'view',
        '-H',                 $self->alignment_file_path,
        '|',                  'grep @SQ',
        '|',                  'cut -f2,3',
        '|',                  "sed 's/[SL]N://g'",
        '|',                  'awk',
        '-F $\'\t\'',         "'BEGIN{OFS=FS}{print \$1,0,\$2,1}'",
    );

    push @cmd, '|', 'sort -k1,1 -k2,2n' if ( $self->sort );
    push @cmd, '>', $target;
    my $cmd = join( ' ', @cmd );

    $self->do_system($cmd);
}

sub temp_dir {
    my ($self) = @_;

    my %opts = ( CLEANUP => $self->cleanup_temp );

    if ( $self->working_dir ) {
        $opts{DIR} = $self->working_dir;
    }

    if ( !$self->cleanup_temp ) {
        $self->log("Temp dir will not be cleaned up at program exit");
    }

    return tempdir(%opts);
}

sub _extend_length {
    my ($self) = @_;

    return int( ( $self->fragment_length - $self->read_length ) / 2 );
}

sub _shift_size {
    my ($self) = @_;

    return int( $self->fragment_length / 2 );
}

sub _do_cleanup {
    my ( $self, $tmp_dir ) = @_;

    if ( $self->cleanup_temp ) $self->do_system( 'rm -rf ' . $tmp_dir );
}

sub _expectation_correction_factor {
    my ($self) = @_;

    my $rc = $self->read_count;
    my $rl = $self->read_length;
    my $m  = $self->mappability_auc;

    my $ec = ( ( $rc * $rl ) / $m );

    $self->log(
"read_count: $rc read_length: $rl total_mappability: $m correction_factor: $ec"
    );

    return $ec;
}

sub _load_bed_chrom_lengths {
    my ($self) = @_;

    my %c;
    open my $fh, '<', $self->chrom_sizes_bed_file_path;
    while (<$fh>) {
        next if m/^#/;
        chomp;
        my ( $seq, $start, $end ) = split /\t/;

        $c{$seq} = $end;
    }
    close $fh;

    return \%c;
}

sub _prep {
    my ( $self, $tmp_dir ) = @_;

    if ( $self->mappability_auc_file_path && !$self->mappability_auc ) {
        $self->log(
            'Getting mappability AUC from file:',
            $self->mappability_auc_file_path
        );
        my $mappability =
          capture( 'head -1 ' . $self->mappability_auc_file_path );
        chomp $mappability;
        $self->mappability_auc($mappability);
    }

    if ( !$self->mappability_auc && $self->mappability_file_path ) {
        $self->log( 'Calculating mappability AUC from mappability file:',
            $self->mappability_file_path );
        my $mappability = capture(
            $self->wiggletools_path . ' AUC ' . $self->mappability_file_path );
        chomp $mappability;
        $self->mappability_file_path($mappability);
    }

    if ( !$self->read_count ) {
        $self->log( 'Getting read count from alignment',
            $self->alignment_file_path );
        my $read_count =
          capture( $self->samtools_path
              . ' view -c -F4 '
              . $self->alignment_file_path );
        chomp $read_count;
        $self->read_count($read_count);
    }

    if ( !$self->smooth_length ) {
        $self->log("Setting smooth length to fragment length");
        $self->smooth_length( $self->fragment_length );
    }

    if ( !$self->chrom_sizes_bed_file_path ) {
        my $chroms = $tmp_dir . '/chroms.bed';
        $self->log( 'Creating chroms file:', $chroms );
        $self->_create_chrom_bed_file($chroms);
        $self->chrom_sizes_bed_file_path($chroms);
    }

    if ( !$self->read_length ) {
        $self->log( 'Getting read length from first read in alignment:',
            $self->alignment_file_path );
        my $cmd = join( ' ',
            $self->samtools_path, 'view', $self->alignment_file_path, '|',
            'head -1', '|', 'cut -f10' );
        my $first_read = capture($cmd);

        chomp($first_read);
        $self->log( 'First read length:',
            length($first_read), '(', $first_read, ')' );
        $self->read_length( length($first_read) );
    }
}

sub log {
    my ( $self, @msg ) = @_;

    print STDERR join( ' ', @msg ) . "\n" if ( $self->verbose );
}

1;
