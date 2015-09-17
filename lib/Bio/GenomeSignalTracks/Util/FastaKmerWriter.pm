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
package Bio::GenomeSignalTracks::Util::FastaKmerWriter;

use strict;
use Moose;
use namespace::autoclean;
use Moose::Util::TypeConstraints;

has 'kmer_size' => ( is => 'rw', isa => 'Int' );
has 'out_fh'    => ( is => 'rw', isa => 'FileHandle' );

has 'in_fh'          => ( is => 'rw', isa => 'FileHandle' );
has 'max_num_seqs'   => ( is => 'rw', isa => 'Maybe[Int]' );
has 'first_seq_name' => ( is => 'rw', isa => 'Maybe[Str]' );
has 'output_type' => (is => 'rw', isa => enum([qw(FASTA RAW)]), default => 'FASTA' );
has 'current_seq_name' => (is => 'rw', isa => 'Maybe[Str]');
has 'seqs_read' => (is => 'rw', isa => 'Int', default => 0);

no Moose::Util::TypeConstraints;


sub split {
    my ($self) = @_;

    my $in_fh          = $self->in_fh();
    my $out_fh         = $self->out_fh();
    my $kmer_size      = $self->kmer_size();
    my $max_num_seqs   = $self->max_num_seqs();
    my $first_seq_name = $self->first_seq_name();
    my $output_type = $self->output_type();
    
    my $seq_count = 0;

    local $/ = '>';
    while (<$in_fh>) {
        chomp;
        if ($_) {
            my $seq    = $_;
            my $nl_idx = index( $seq, "\n" );
            my $ws_idx = index( $seq, ' ' );
            my $name;

            if ( $seq_count == 0 && $first_seq_name ) {
                $name = $first_seq_name;
            }
            else {
                $name = substr( $seq, 0, $nl_idx );
                ($name) = split /\s/, $name;
            }
            $self->current_seq_name($name);

            $seq =~ s/\n//g;
            my $limit = length($seq) - $kmer_size + 1;

            for ( my $i = $nl_idx ; $i < $limit ; $i++ ) {
              my $kmer = substr( $seq, $i, $kmer_size );
              if ($output_type eq 'FASTA'){
                my $start = $i + 1 - $nl_idx;
                my $end   = $i + $kmer_size - $nl_idx;

                print $out_fh ">$name:$start-$end\n$kmer\n";
              }
              else {
                print $out_fh $kmer."\n";
              }
               
            }

            $seq_count += 1;
            $self->seqs_read($seq_count);
            return if ( $max_num_seqs && $seq_count >= $max_num_seqs );

        }
    }
}

__PACKAGE__->meta->make_immutable;
1;

