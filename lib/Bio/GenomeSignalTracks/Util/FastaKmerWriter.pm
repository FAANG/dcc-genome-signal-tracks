package Bio::GenomeSignalTracks::Util::FastaKmerWriter;

use strict;
use Moose;

has 'kmer_size' => ( is => 'rw', isa => 'Int' );
has 'out_fh'    => ( is => 'rw', isa => 'FileHandle' );

has 'in_fh'        => ( is => 'rw', isa => 'FileHandle' );
has 'state'        => ( is => 'rw', isa => 'Int', default => 0 );
has 'seqname'      => ( is => 'rw', isa => 'Maybe[Str]' );
has 'next_seqname' => ( is => 'rw', isa => 'Maybe[Str]' );

sub split {
    my ($self) = @_;

    my $out_fh    = $self->out_fh();
    my $kmer_size = $self->kmer_size();

    while ( $self->next_seq() ) {
        my $buffer;

        while ( my $line = $self->next_line() ) {
            $buffer .= $line;

            if (length($buffer) < $kmer_size) {
              next;
            }
            
            my $limit = length($buffer) - $kmer_size + 1;

            my $i = 0;
            for ( ; $i < $limit ; $i++ ) {
                my $kmer = substr( $buffer, $i, $kmer_size );
                print $out_fh $kmer . $/;
            }
            substr( $buffer, 0, $limit ) = "";
        }
    }
}

sub next_seq {
    my ($self) = @_;

    #state 0 - awaiting seq
    #state 1 - found next seq
    #state 2 - reading seq
    #state 3 - done

    my $state = $self->state();

    while ( $self->state() == 0 || $self->state() == 2 ) {
        $self->next_line();
    }
    if ( $self->state() == 1 ) {
        $self->seqname( $self->next_seqname() );
        $self->next_seqname(undef);
        return 1;
    }
    if ( $self->state() == 3 ) {
        return undef;
    }
}

sub next_line {
    my ($self) = @_;

    my $fh = $self->in_fh();
    if ( my $line = <$fh> ) {
        chomp $line;

        if ( substr( $line, 0, 1 ) eq '>' ) {
            $self->state(1);
            $self->next_seqname( substr( $line, 1 ) );
            return undef;
        }
        else {
            return $line;
        }
    }
    else {
        $self->state(3);
        return undef;
    }
}

1;
