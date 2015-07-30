package Bio::GenomeSignalTracks::Util::FastaKmerWriter;

use strict;
use Moose;

has 'kmer_size' => ( is => 'rw', isa => 'Int' );
has 'out_fh'    => ( is => 'rw', isa => 'FileHandle' );

has 'in_fh'   => ( is => 'rw', isa => 'FileHandle' );
has 'seqname' => ( is => 'rw', isa => 'Maybe[Str]' );

sub split {
    my ($self) = @_;

    my $in_fh     = $self->in_fh();
    my $out_fh    = $self->out_fh();
    my $kmer_size = $self->kmer_size();

    local $/ = '>';
    while (<$in_fh>) {
        chomp;
        if ($_) {
            my $name   = $1;
            my $seq    = $_;
            my $nl_idx = index( $seq, "\n" );

            $self->seqname( substr( $seq, 0, $nl_idx ) );
            $seq =~ s/\n//g;
            my $limit = length($seq) - $kmer_size + 1;

            for ( my $i = $nl_idx ; $i < $limit ; $i++ ) {
                print $out_fh substr( $seq, $i, $kmer_size ) ."\n";
            }
        }
    }
}

1;

