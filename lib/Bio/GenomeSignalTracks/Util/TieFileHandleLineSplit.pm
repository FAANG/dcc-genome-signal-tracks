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
package Bio::GenomeSignalTracks::Util::TieFileHandleLineSplit;

=head1 NAME

Tie::FileHandle::LineSplit - Filehandle tie that captures, splits and stores output into files in a given path.

=head1 VERSION

Version 0.1

=cut

$VERSION = 0.1;

=head1 DESCRIPTION

This module, when tied to a filehandle, will capture and store all that
is output to that handle. You should then select a path to store files and a
size to split files.

=head1 SYNOPSIS

 # $path should exist or the current process have enough priv. for creation.
 # $size should be > 0.
 tie *HANDLE, 'Tie::FileHandle::LineSplit', $path, $size;

 # Register code to listen to file creation
 (tied *HANDLE)->add_file_creation_listeners( sub {
	my ( $tied_object, $filename) = @_;
	print "Created $filename with size: " . -s $filename . "\n";
 } );

 # Will create int( $many_times_size / $split_size) files of size $split_size.
 # Will call each listener int( $many_times_size / $split_size) times.
 # Buffers will hold $many_times_size % $split_size outstanding bytes.
 (tied *HANDLE)->print( ' ' x $many_times_size );

 # Write all outstanding output from buffers to files.
 # The last file created can be smaller than split_size
 (tied *HANDLE)->write_buffers;

 # Get generated filenames to the moment
 (tied *HANDLE)->get_filenames();

=cut

use 5.10.0;
use strict;
use warnings;

use vars qw(@ISA $VERSION);
use base qw(Tie::FileHandle::Base);

use File::Path;
use IO::Compress::Gzip;
use File::Temp;
use Carp;

# Tie::FileHandle implementation
# Usage: tie *HANDLE, 'Tie::FileHandle::Split', $path, $split_size
sub TIEHANDLE {
    my ( $class, $path, $split_size, $gzip, $filenamer_callback ) = @_;

    my $self = {
        class              => $class,
        path               => $path,
        split_size         => $split_size,
        current_fh         => '',
        num_writes         => 0,
        filenames          => [],
        listeners          => {},
        gzip               => $gzip,
        filenamer_callback => $filenamer_callback,
    };

    File::Path::make_path( $self->{path} ) unless -d $self->{path};

    bless $self, $class;
}

# Tie::FileHandle implementation
# Print to the selected handle
sub PRINT {
    my ( $self, $data ) = @_;

    if ( !$self->{current_fh} ) {
        $self->_create_file();
    }
    my $fh = $self->{current_fh};

    print $fh $data;

    $self->{num_writes}++;
    if ( $self->{num_writes} >= $self->{split_size} ) {
        close( $self->{current_fh} );
        $self->{current_fh} = undef;
        $self->{num_writes} = 0;
    }
}

sub CLOSE {
    my ($self) = @_;

    if ( $self->{current_fh} ) { return close( $self->{current_fh} ); }
    return 1;

}

sub _create_file {
    my ($self) = @_;

    if ($self->{current_fh}) {
      close( $self->{current_fh} );
    }

    my $filenamer_callback = $self->{filenamer_callback};
    my $file_name  = $filenamer_callback->();
    
    my $dir = $self->{path};
    
    my $path = "$dir/$file_name";
    my $op = $self->{gzip} ? '>:gzip' : '>';
    
    open(my $fh, $op, $path);  

    $self->{current_fh} = $fh;

    # Call listeners
    foreach my $listener ( keys %{ $self->{listeners} } ) {
        &{ $self->{listeners}->{$listener} }( $self, $path );
    }

    push @{ $self->{filenames} }, $path;
}

=head1 METHODS

=head3 C<get_filenames>

C<get_filenames> returns a list of the files generates until the moment of the
call. It should be used to get the names of files and rename them to the
desired filenames. In a standard splitting operation C<get_filenames> is
called after outputting all data to the filehandle and calling C<write_buffers>.

=cut

# Returns filenames generated up to the moment the method is called
sub get_filenames {
    my ($self) = @_;

    return @{ $self->{filenames} } if defined $self->{filenames};
}

=head3 C<add_file_creation_listeners>

C<add_file_creation_listeners> adds methods to the list of listeners of the
file creation event. Methods should be code, array, arrayref or any
non-recursive structure resulting from them. Since methods are added to a HASH,
several elements pointing to the same piece of code will be added only once.
Code observing this event is called once per file created of the $split_size
size defined in the tie clause. When called the Tie::FileHandle::Split object
and the complete path to the newly created file is passed as parameter. The file
is of the specified C<$split_size> defined in the tie clause unless generated
from a C<write_buffers> call, has been closed and an effort has been made for it
to sync (untested).

=cut

sub add_file_creation_listeners {
    my ( $self, @listeners ) = @_;

    foreach my $listener (@listeners) {
        if ( ref($listener) eq 'CODE' ) {
            $self->{listeners}->{$listener} = $listener;
        }
        elsif ( ref($listener) eq 'ARRAY' ) {
            $self->add_file_creation_listeners(@$listener);
        }
        elsif ( ref($listener) eq 'ARRAYREF' ) {
            $self->add_file_creation_listeners($listener);
        }
        else {
            croak(  "Unsupported structure in add_file_creation_listeners. "
                  . "Can use any structure containing CODE, ARRAY and ARRAYREF. "
                  . "Looks like a "
                  . ref($listener) );
        }
    }
}

=head3 C<remove_file_creation_listeners>

C<remove_file_creation_listeners> removes a list of methods from the list of
listeners of the file creation event. Methods should be code, array, arrayref or
any non-recursive structure resulting from them.

=cut

sub remove_file_creation_listeners {
    my ( $self, @listeners ) = @_;

    foreach my $listener (@listeners) {
        if ( ref($listener) eq 'CODE' ) {
            delete $self->{listeners}->{$listener};
        }
        elsif ( ref($listener) eq 'ARRAY' ) {
            $self->remove_file_creation_listeners(@$listener);
        }
        elsif ( ref($listener) eq 'ARRAYREF' ) {
            $self->remove_file_creation_listeners($listener);
        }
        else {
            croak(  "Unsupported structure in add_file_creation_listeners. "
                  . "Can use any structure containing CODE, ARRAY and ARRAYREF. "
                  . "Looks like a "
                  . ref($listener) );
        }
    }
}

=head3 C<clear_file_creation_listeners>

C<clear_file_creation_listeners> removes all methods from the list of listeners
of the file creation event.

=cut

sub clear_file_creation_listeners {
    my ($self) = @_;

    $self->{listeners} = {};
}

sub _get_listeners {
    my ($self) = @_;

    # Behold! Dereferencing fixes incompatibility with pre 5.14 perl.
    # Both keys and each are affected if a hashref is passed.
    return map $_, keys %{ $self->{listeners} };
}

1;

=head1 TODO

=over 4

=item * Very untested for anything other than writing to the filehandle.

=item * write_buffers should sync to disk, untested and seeking advice.

=back

=head1 BUGS

No known bugs. Please report and suggest tests to gbarco@cpan.org.

=cut

=head1 AUTHORS AND COPYRIGHT

Written by Gonzalo Barco based on Tie::FileHandle::Buffer written by Robby Walker ( robwalker@cpan.org ).

Project repository can be found at https://github.com/gbarco/Tie-FileHandle-Split.

You may redistribute/modify/etc. this module under the same terms as Perl itself.

