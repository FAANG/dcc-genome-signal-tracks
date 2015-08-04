package Bio::GenomeSignalTracks::Util::FastaKmerWriter;

use strict;
use Moose;

has 'kmer_size' => ( is => 'rw', isa => 'Int' );
has 'out_fh'    => ( is => 'rw', isa => 'FileHandle' );

has 'in_fh'   => ( is => 'rw', isa => 'FileHandle' );
has 'max_num_seqs' => (is => 'rw', isa => 'Maybe[Int]');

sub split {
    my ($self) = @_;

    my $in_fh     = $self->in_fh();
    my $out_fh    = $self->out_fh();
    my $kmer_size = $self->kmer_size();
    my $max_num_seqs = $self->max_num_seqs();
    
    my $seq_count = 0;
    
    local $/ = '>';
    while (<$in_fh>) {
        chomp;
        if ($_) {
            my $seq    = $_;
            my $nl_idx = index( $seq, "\n" );
            

            $seq =~ s/\n//g;
            my $limit = length($seq) - $kmer_size + 1;

            for ( my $i = $nl_idx ; $i < $limit ; $i++ ) {
                print $out_fh substr( $seq, $i, $kmer_size ) ."\n";
            }
            
            $seq_count++;
            return if ($max_num_seqs && $seq_count >= $max_num_seqs);
        }
    }
}

1;

