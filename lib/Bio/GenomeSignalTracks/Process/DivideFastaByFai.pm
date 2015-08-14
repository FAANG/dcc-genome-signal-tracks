package Bio::GenomeSignalTracks::Process::DivideFastaByFai;
use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use autodie;

sub fetch_input {
    my ($self) = @_;

    my $fai               = $self->param_required('fai');
    my $target_base_pairs = $self->param('target_base_pairs');
    my $fan_branch_code   = $self->param_required('fan_branch_code');

    print STDERR Dumper( [ $fai, $target_base_pairs, $fan_branch_code ] );

    my @fai_entries;
    my @fai_columns = qw(seq_name length offset line_length line_byte_length);

    open( my $fh, '<', $fai );
    while (<$fh>) {
        chomp;
        my %fai_entry;
        @fai_entry{@fai_columns} = split /\t/;

        push @fai_entries, \%fai_entry;
    }
    close($fh);

    $self->param( 'fai_entries', \@fai_entries );
}

sub param_defaults {
    return {
        fan_branch_code   => 2,
        target_base_pairs => 5_000_000,
    };
}

sub write_output {
    my ($self) = @_;

    my $target_base_pairs = $self->param_required('target_base_pairs');
    my $fai_entries       = $self->param_required('fai_entries');
    my $fan_branch_code   = $self->param_required('fan_branch_code');

    my @fai_bundles;
    my $current_bundle;

    for my $fai (@$fai_entries) {
      
        if (  !$current_bundle
            || $fai->{length} + $current_bundle->{total_length} >
            $target_base_pairs )
        {
            $current_bundle = {
                start_pos    => $fai->{offset},
                seq_count    => 1,
                total_length => $fai->{length},
            };
            push @fai_bundles, $current_bundle;
            next;
        }

        $current_bundle->{seq_count}++;
        $current_bundle->{total_length} =
          $fai->{length} + $current_bundle->{total_length};
    }

    for my $bundle (@fai_bundles) {
        my $seq_start_pos = $self->dataflow_output_id(
            {
                seq_start_pos    => $bundle->{start_pos},
                num_seqs_to_read => $bundle->{seq_count}
            },
            $fan_branch_code
        );
    }

}

1;
