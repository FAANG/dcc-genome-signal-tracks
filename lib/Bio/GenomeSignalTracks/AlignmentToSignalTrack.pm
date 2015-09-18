
=head1 LICENSE

   Copyright 2015 EMBL - European Bioinformatics Institute
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

=cut

package Bio::GenomeSignalTracks::AlignmentToSignalTrack;

=head1 NAME

AlignmentToSignalTrack - Convert an aligmment to a normalized signal track

=head1 SYNOPSIS

    use Bio::GenomeSignalTracks::AlignmentToSignalTrack;

    my $a2s = Bio::GenomeSignalTracks::AlignmentToSignalTrack->new(
      samtools_path => '/path/to/samtools',
      bedtools_path => '/path/to/bedtools',
      wiggletools_path => '/path/to/wiggletools',
      mappability_file_path => '/path/to/mappability.bw',
      exclude_regions_file_path => '/path/to/regions_to_ignore.bed',
      alignment_file_path => '/path/to/alignment.bam',
      fragment_length => 200,
      output_file_path => '/path/to/output.wig',
    );

    $a2s->generate_track();

=head1 Description

    This module converts alignments to signal tracks, normalized by read count
    and mappability. It relies on the samtools, bedtools and wiggletools

    It is inspired by Anshul Kundaje's wiggler/align2rawsignal software
      * https://code.google.com/p/align2rawsignal/
      * https://sites.google.com/site/anshulkundaje/projects/wiggler

    The input is a BAM file and a mappability track.  

=cut

use strict;
use Moose;
use Carp;
use File::Path qw(remove_tree);
use File::Temp qw(tempdir);
use File::Copy qw(move);
use File::Basename qw(fileparse);
use IPC::System::Simple qw(system capture);
use POSIX qw(mkfifo);
use autodie;
use Moose::Util::TypeConstraints;
use namespace::autoclean;

subtype 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::PositiveInt',
  as 'Int',
  where { $_ > 0 },
  message { "The number you provided, $_, was not a positive number" };

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

=head1 Attributes

=head2 samtools_path

    path to the samtools executable. Required

=cut

has 'samtools_path' => (
    is => 'rw',
    isa =>
      'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Executable',
    required => 1,
);

=head2 wiggletools_path

    path to the wiggletools executable. Required. 

=cut

has 'wiggletools_path' => (
    is => 'rw',
    isa =>
      'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Executable',
    required => 1,
);

=head2 bedtools_path

    path to the bedtools executable. Required. 

=cut

has 'bedtools_path' => (
    is => 'rw',
    isa =>
      'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Executable',
    required => 1,
);

=head2 _chrom_sizes_bed_file_path

    path to a bed file listing all chromosomes in the alignment header. This is
    used to prevent the output exceeding the chromsome length. This will be 
    created from the alignment file header, there is no need to provide it.

=cut

has '_chrom_sizes_bed_file_path' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable'
);

=head2 mappability_file_path

    path to a signal track listing mappability scores. This should be in a 
    format wiggletools can read (we use bigwig). The signal values will be
    divided by this value (hard to map regions have their score 
    increased). We have tested this with a score range of 0-1, tracks 
    produced with https://github.com/FAANG/reference_data_builder or
    https://github.com/FAANG/genome-signal-tracks

=cut

has 'mappability_file_path' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable',
    required => 1,
);

=head2 mappability_auc

    total mappable area in the mappability file. This is used to normalise the
    signal. If not provided, it will be read from mappability_auc_file_path, 
    or calculated from the mappabilty files with wiggletools AUC.

=cut

has 'mappability_auc' => ( is => 'rw', isa => 'Num' );

=head2 mappability_auc_file_path

    Source for the mappability_auc value. Should be a single numeric value in
    the first line. These files are routinely produced by the pipeline
    https://github.com/FAANG/reference_data_builder

=cut

has 'mappability_auc_file_path' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable',
);

=head2 fragment_length

    Reads are shifted by half the fragment length when producing signal data.
    Required.

=cut

has 'fragment_length' => (
    is       => 'rw',
    isa      => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::PositiveInt',
    required => 1,
);

=head2 alignment_file_path

    Path to the source alignment file. Must be coordinate sorted, with a header. Tested with BAM, but SAM or CRAM should also work, depending on your version of samtools. Required. 

=cut

has 'alignment_file_path' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable',
    required => 1,
);

=head2 output_file_path

    Target output destination. Required. 

=cut

has 'output_file_path' => ( is => 'rw', isa => 'Str', required => 1, );

=head2 read_count

    Number of reads in alignment. Used when normalising the signal. 
    This will be calculated from the alignment if not provided.  

=cut

has 'read_count' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::PositiveInt'
);

=head2 read_length

    Read length is used when normalising the signal. Will be calculated  
    from the first read if not provided.

=cut

has 'read_length' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::PositiveInt'
);

=head2 working_dir

    The module creates a temp directory for intermediate files, using 
    File::Temp. If set, the temp directory will be created within the
    working_dir.

=cut

has 'working_dir' => ( is => 'rw', isa => 'Maybe[Str]' );

=head2 smooth_length

    Signal track output is smoothed by this value, using wiggletools. If not
    provided, the fragment_length is used.

=cut

has 'smooth_length' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::PositiveInt'
);

=head2 exclude_regions_file_path

    Path to a bed file containing regions to exclude. To emulate 
    align2rawsignal, use the set of regions where mappability is less than 0.25
    
    These files are created by the pipeline https://github.com/FAANG/reference_data_builder

=cut

has 'exclude_regions_file_path' => (
    is  => 'rw',
    isa => 'Bio::GenomeSignalTracks::AlignmentToSignalTrack::FilePath_Readable'
);

=head2 force_lexographic_sort

    wiggletools requires input to be lexicographically sorted. We can guarantee
    this by running everything through sort -k1,1 -k2,2n, and will attempt to 
    do this if the alignment chromosomes are not in the required order. 

    If force_lexographic_sort is true, we skip trying to work it out and just
    sort everything.    

=cut

has 'force_lexographic_sort' => ( is => 'rw', isa => 'Bool', default => 1, );

=head2 output_precision_dp
  
    wiggletools output is post-processed to remove negative 0 values, and to
    control the precision of the output. If set, this controls how many decimal
    places to include in the output.
  
=cut  

has 'output_precision_dp' => (
    is => 'rw',
    isa =>
      'Maybe[Bio::GenomeSignalTracks::AlignmentToSignalTrack::PositiveInt]',
    default => 2
);

=head2 output_type
  
    output format. can be wig (wiggle) or bg (bedgraph). Required, defaults 
    to wig.
  
=cut  

has 'output_type' => (
    is       => 'rw',
    isa      => enum( [qw(wig bg)] ),
    default  => 'wig',
    required => 1,
);

=head2 intermediate_files
  
    Should intermediate files be named pipes (fifo) or temporary (tmp) files? 
    Named pipes should use less temp space where a sort is not required, with
    increased RAM use, as more processes run concurrently.

    Required, default is auto - will use tmp if a lexographic sort is required.
  
=cut  

has 'intermediate_files' => (
    is       => 'rw',
    isa      => enum( [qw(tmp fifo auto)] ),
    default  => 'tmp',
    required => 1,
);

=head2 verbose
  
    When true, the module will spam STDERR with logging information.

    Required, defaults to true.
  
=cut  

has 'verbose' => ( is => 'rw', isa => 'Bool', default => 1 );

=head2 cleanup_temp
  
    Is the temp directory deleted when the program exits?

    Required, defaults to true.
  
=cut

has 'cleanup_temp' => ( is => 'rw', isa => 'Bool', default => 1 );

no Moose::Util::TypeConstraints;

#todo - check that you have sufficient information to run the process up front
#todo - check that you can write to output location early
#todo - can you apply this to paired end data?
#todo - documentation
#todo - do we need to force sort collate locale to C ?

=head1 Methods

=head2 generate_track
  
  Create a normalized signal track, based on the aligment file.
  
=cut

sub generate_track {
    my ($self) = @_;

    $self->_check_requirements;

    my $tmp     = File::Temp->newdir( $self->_temp_dir_options );
    my $tmp_dir = $tmp->dirname;

    $self->log( 'Temp dir:', $tmp_dir );
    $self->_prep($tmp_dir);
    $self->log('Prep complete');

    my $fwd_scaled_output = $tmp_dir . '/fwd.scaled.shifted.bg';
    my $rev_scaled_output = $tmp_dir . '/rev.scaled.shifted.bg';

    my $combined_output = $tmp_dir . '/combined.out';

    my $split_options = [
        [ '-F16', $fwd_scaled_output, '+' ],
        [ '-f16', $rev_scaled_output, '-' ],
    ];

    if ( $self->intermediate_files eq 'fifo' ) {
        $self->_run_with_fifos(
            $fwd_scaled_output, $rev_scaled_output,
            $split_options,     $combined_output
        );
    }

    if ( $self->intermediate_files eq 'tmp' ) {
        $self->_run_with_temp_files(
            $fwd_scaled_output, $rev_scaled_output,
            $split_options,     $combined_output
        );
    }

    $self->log( 'Move output to final location:', $self->output_file_path );
    move( $combined_output, $self->output_file_path );
    $self->log( 'Output in', $self->output_file_path );
}

=head2 _check_requirements

 croak if we don't have enough information to proceed.
 
 only used for things that can't be managed through Moose settings.

=cut

sub _check_requirements {
    my ($self) = @_;

    #can we write to the output file? - autodie will croak for us if this fails.
    open( my $fh, '>', $self->output_file_path );
    close($fh);

    if ( $self->mappability_file_path ) {
        my $suffix = $self->_suffix( $self->mappability_file_path );
        croak(
"Mappability file must be wiggle (wig), bedGraph (bg) or bigWig file, but suffix is "
              . $suffix )
          if ( $suffix ne 'wig' && $suffix ne 'bg' && $suffix ne 'bw' );
    }

    if ( $self->exclude_regions_file_path ) {
        my $suffix = $self->_suffix( $self->exclude_regions_file_path );
        croak( "Excluded regions file must be in bed format but suffix is "
              . $suffix )
          if ( $suffix ne 'bed' );
    }

    if ( $self->alignment_file_path ) {
        my $suffix = $self->_suffix( $self->alignment_file_path );
        croak( "Alignment file should be .bam, .cram or .sam, but suffix is "
              . $suffix )
          if ( $suffix ne 'bam' && $suffix ne 'cram' && $suffix ne 'sam' );
    }
}

=head2 _suffix  

    returns the suffix from a file name

=cut

sub _suffix {
    my ( $self, $fn ) = @_;

    my ($suffix) = reverse split /\./, $fn;
    
    return $suffix;
}

=head2 _is_sort_required

  We don't need to sort the aligment if the file is coordinate sorted and the
  reference sequences are already in lexographic order.
  
  params - sort order (scalar), list of reference sequence names (array ref)

=cut

sub _is_sort_required {
    my ( $self, $so, $sqs ) = @_;

    my $lex_order = $self->_is_list_in_lex_order($sqs) || 0;
    $self->log( "Sort order:",               $so );
    $self->log( "SQs in lexographic order:", $lex_order );
    if ( $so eq 'coordinate' && $lex_order ) {
        $self->log("Sort is not necessary");
        return 0;
    }
    else {
        $self->log("Sort is necessary");
        return 1;
    }

}

=head2 _is_list_in_lex_order

  Returns true if the first arg (an array ref) is already in the order produced by sort {$a cmp $b}
  
=cut

sub _is_list_in_lex_order {
    my ( $self, $sqs ) = @_;

    my @sorted_copy = sort { $a cmp $b } @$sqs;

    my $order_match = 1;
    for ( my $i = 0 ; $order_match && $i < scalar(@sorted_copy) ; $i++ ) {
        $order_match = ( $sorted_copy[$i] eq $sqs->[$i] );
    }

    return $order_match;
}

=head2 _parse_header

    Parse the alignment header for the sort order and reference sequence information

    Returns (sort order, array ref to sequence names, hash ref to sequence names and their lengths)

=cut

sub _parse_header {
    my ($self) = @_;

    my $so;
    my @sqs;
    my %sq_lengths;
    my $cmd = join( ' ',
        $self->samtools_path, 'view', '-H', $self->alignment_file_path );
    open( my $aligment_header_fh, '-|', $cmd );
    while (<$aligment_header_fh>) {
        chomp;
        my @tokens = split /\t/;
        if ( scalar(@tokens) && $tokens[0] eq '@HD' ) {
            ($so) = grep { /^SO:/ } @tokens;
            $so =~ s/^SO://;
        }
        if ( scalar(@tokens) && $tokens[0] eq '@SQ' ) {
            my ($sn) = grep { /^SN:/ } @tokens;
            my ($ln) = grep { /^LN:/ } @tokens;
            $sn =~ s/^SN://;
            $ln =~ s/^LN://;

            push @sqs, $sn;
            $sq_lengths{$sn} = $ln;
        }
    }

    close($aligment_header_fh);

    return ( $so, \@sqs, \%sq_lengths );
}

=head2 _run_with_fifos

  Generate the track using named pipe intermediate files
  
=cut

sub _run_with_fifos {
    my ( $self, $fwd_scaled_output, $rev_scaled_output, $split_options,
        $combined_output )
      = @_;

    $self->log("Passing intermediate data through fifos");
    my $fifo_mode = 0600;
    mkfifo( $fwd_scaled_output, $fifo_mode );
    mkfifo( $rev_scaled_output, $fifo_mode );

    for (@$split_options) {
        my $pid = fork;

        if ( !defined $pid ) {
            die "Cannot fork: $!";
        }
        elsif ( $pid == 0 ) {

   #write intermediate data to fifos - process will block until there's a reader
            $self->_split_by_dir_to_bed(@$_);
            exit 0;
        }
    }

    # this now reads from the FIFOs, unblocking their inputs
    $self->_combine_extend_smooth( $combined_output, $fwd_scaled_output,
        $rev_scaled_output );
}

=head2 _run_with_temp_files

  Generate the track using temporary intermediate files
  
=cut

sub _run_with_temp_files {
    my ( $self, $fwd_scaled_output, $rev_scaled_output, $split_options,
        $combined_output )
      = @_;

    $self->log("Passing intermediate data through tmp files");

    for (@$split_options) {
        $self->log( "Split alignment by dir", @$_ );
        $self->_split_by_dir_to_bed(@$_);
    }
    $self->log( 'Combine, extend and smooth all reads:', $combined_output );
    $self->_combine_extend_smooth( $combined_output, $fwd_scaled_output,
        $rev_scaled_output );
}

=head2 _split_by_dir_to_bed

    Take reads from one strand of the alignment, generate normalized signal for
    it, shifted based on the fragment length.

    Depending on the settings, it can also exclude reads based on the file exclude_regions_file_path.

    Data will be sorted lexographically if force_lexographic_sort is set.

=cut

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

    push @cmd, '|', 'sort -k1,1 -k2,2n' if ( $self->force_lexographic_sort );

    push @cmd, '|', $self->wiggletools_path, 'write_bg -', 'ratio strict',
      '-',
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

    $self->_do_system($cmd);
}

=head2 _do_system 
  
    Run a command through system, wrapped up to use bash and the pipefail option

=cut

sub _do_system {
    my ( $self, $cmd ) = @_;
    $self->log( "executing", $cmd );

    system( 'bash', '-o', 'pipefail', '-c', $cmd );
}

=head2 _combine_extend_smooth 
  
    Combine data from multiple files, while extending and smoothing it.
    Data is post-processed 

    params - target file name, list of files to combine
  
=cut

sub _combine_extend_smooth {
    my ( $self, $target, @files ) = @_;

    my $output_cmd = 'write';
    if ( $self->output_type eq 'bg' ) {
        $output_cmd = 'write_bg';
    }

    my $cmd = join( ' ',
        $self->wiggletools_path,           $output_cmd,
        '-',                               'mult strict',
        $self->_chrom_sizes_bed_file_path, 'smooth',
        $self->smooth_length,              'extend',
        $self->_extend_length,             'sum',
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
            $vals[3] = $self->_post_process_score( $vals[3] );
        }
        if ( scalar(@vals) == 2 ) {

            #variableStep format
            $vals[1] = $self->_post_process_score( $vals[1] );
        }
        if ( scalar(@vals) == 1 ) {

            #fixedStep format
            $vals[0] = $self->_post_process_score( $vals[0] );
        }

        print $out_fh join( "\t", @vals ) . "\n";
    }

    close($in_fh);
    close($out_fh);
}

=head2 _post_process_score

    Replace -0.00000 with 0.
    If output_precision_dp is set, round the scores down to that many decimal points

=cut

sub _post_process_score {
    my ( $self, $score ) = @_;

    if ( $score =~ m/^-0.00000/ ) {
        $score = 0;
    }
    if ( defined $self->output_precision_dp ) {
        $score = sprintf '%.' . $self->output_precision_dp . 'f', $score;
    }
    return $score;
}

=head2 _create_chrom_bed_file
  
    Create a bed file, based on a hash listing reference sequences and their
    lengths

=cut

sub _create_chrom_bed_file {
    my ( $self, $target, $sq_lengths ) = @_;

    open( my $fh, '>', $target );

    for my $sq ( sort keys %$sq_lengths ) {
        my $length = $sq_lengths->{$sq};

        print $fh join( "\t", $sq, 0, $length, 1 ) . "\n";
    }

    close $fh;
}

=head2 _temp_dir_options
  
    Provide options hash for creating a temp dir

=cut

sub _temp_dir_options {
    my ($self) = @_;

    my %opts = ( CLEANUP => $self->cleanup_temp );

    if ( $self->working_dir ) {
        $opts{DIR} = $self->working_dir;
    }

    if ( !$self->cleanup_temp ) {
        $self->log("Temp dir will not be cleaned up at program exit");
    }

    return %opts;
}

=head2 _extend_length
  
    By how many base pairs should reads be extended? 
    Calculated as (fragment_length - read_length)/2, rounded down.

=cut

sub _extend_length {
    my ($self) = @_;

    return int( ( $self->fragment_length - $self->read_length ) / 2 );
}

=head2 _shift_size
  
    By how many base pairs should reads be shifted? 
    Calculated as fragment_length/2, rounded down.

=cut

sub _shift_size {
    my ($self) = @_;

    return int( $self->fragment_length / 2 );
}

=head2 _expectation_correction_factor
  
    Factor for normalizing signal, based on how much of the genome is mappable, how many reads you have, and how long they are.

=cut

sub _expectation_correction_factor {
    my ($self) = @_;

    my $rc = $self->read_count;
    my $rl = $self->read_length;
    my $m  = $self->mappability_auc;

    my $ec = ( ( $rc * $rl ) / $m );

    $self->log(
        "read_count:",        $rc, "read_length:",       $rl,
        "total_mappability:", $m,  "correction_factor:", $ec
    );

    return $ec;
}

=head2 _prep
  
    Prepare for generating a track

=cut

sub _prep {
    my ( $self, $tmp_dir ) = @_;

    my ( $sort_order, $sqs, $sq_lengths ) = $self->_parse_header;

    my $sort_required = $self->_is_sort_required( $sort_order, $sqs );

    if ($sort_required) {
        $self->log('Will sort input data');
        $self->force_lexographic_sort(1);
    }

    if ( $self->intermediate_files eq 'auto' ) {
        if ( $self->force_lexographic_sort ) {
            $self->intermediate_files('tmp');
        }
        else {
            $self->intermediate_files('fifo');
        }
        $self->log( 'Will use', $self->intermediate_files,
            'intermediate files' );
    }

    if ( $self->mappability_auc_file_path && !$self->mappability_auc ) {
        $self->_read_auc_from_file;
    }

    if ( !$self->mappability_auc && $self->mappability_file_path ) {
        $self->_calc_auc_from_mappability_file;
    }

    if ( !$self->read_count ) {
        $self->_get_read_count;
    }

    if ( !$self->smooth_length ) {
        $self->log("Setting smooth length to fragment length");
        $self->smooth_length( $self->fragment_length );
    }

    if ( !$self->read_length ) {
        $self->_get_read_length;
    }

    my $chroms = $tmp_dir . '/chroms.bed';
    $self->log( 'Creating chroms file:', $chroms );
    $self->_create_chrom_bed_file( $chroms, $sq_lengths );
    $self->_chrom_sizes_bed_file_path($chroms);
}

=head2 _read_auc_from_file 

  get the mappabilty AUC from a value stored in a file

=cut

sub _read_auc_from_file {
    my ($self) = @_;
    $self->log(
        'Getting mappability AUC from file:',
        $self->mappability_auc_file_path
    );
    my $mappability = capture( 'head -1 ' . $self->mappability_auc_file_path );
    chomp $mappability;
    $self->mappability_auc($mappability);
}

=head2 _calc_auc_from_mappability_file 

    Calculate the mappability AUC with wiggletools.

=cut

sub _calc_auc_from_mappability_file {
    my ($self) = @_;
    $self->log( 'Calculating mappability AUC from mappability file:',
        $self->mappability_file_path );
    my $mappability = capture(
        $self->wiggletools_path . ' AUC ' . $self->mappability_file_path );
    chomp $mappability;
    $self->mappability_file_path($mappability);
}

=head2 _get_read_count 

    Get the number of mapped reads from the alignment with samtools

=cut

sub _get_read_count {
    my ($self) = @_;

    $self->log( 'Getting read count from alignment',
        $self->alignment_file_path );
    my $read_count =
      capture(
        $self->samtools_path . ' view -c -F4 ' . $self->alignment_file_path );
    chomp $read_count;
    $self->read_count($read_count);
}

=head2 _get_read_length

    Get the first read with samtools and return it's length

=cut

sub _get_read_length {
    my ($self) = @_;
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

=head2 log

    If verbose is true, join the arguments and print to STDERR

=cut

sub log {
    my ( $self, @msg ) = @_;

    print STDERR join( ' ', @msg ) . "\n" if ( $self->verbose );
}

__PACKAGE__->meta->make_immutable;
1;
