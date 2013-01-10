package Physics::Photoemission::Bin;

use Moo;

use Physics::Photoemission::Constants qw/mass/;

=head1 NAME

Physics::Photoemission::Bin

=head1 SYNOPSIS

 my $bin = Physics::Photoemission::Bin->new($start,$end);
 while( --data-- ){
   $bin->add_slice( --data-- );
 }
 say for @{ $bin->result };

=head1 DESCRIPTION

This package implements the spatial bin objects used in the C<Physics::Photoemission> module/simulation. For each time slice, as electrons are generated and propagated forward they are placed into these spatial bins. This slice will have a certain number of electrons and some velocity distribution (currently implemented by mapping to an average velocity --- so be sure to use enough bins so that this approximation is ok). Once all the electrons are placed, calling the C<result()> method will compute the total number of electrons, their average momentum and their momentum uncertainty. The results are cached so that if C<result()> is called multiple times (without adding more slices), the results are only computed the first time.

=head1 METHODS

=head2 new( $start, $end )

The usual constructor method. For convenience (but recommended), the constructor can take a start and end position at creation. Returns a bin object.

=cut

sub BUILDARGS {
  my $class = shift;
  my ($begin, $end) = @_;
  return { begin => $begin, end => $end };
}

=head2 add_slice( $num, $v_begin, $v_end )

Adds data to the bin. This data is the number of electron in the slice and the velocities (lowest and highest) of those electrons. Currently the velocity is simply converted into an average and stored as such, therefore be sure to use enough bins to make the approximation that the these velocities are well represented by an average.

N.B. Adding slice data removes any cached results for the bin, see C<result()>.

=cut

# Method to store slice data in 'slices' array 
sub add_slice {
  my $self = shift;
  my ($num, $v_begin, $v_end) = @_;

  # data to store for each slice
  my $slice_data = [
    $num, ($v_begin + $v_end) / 2,
  ];

  # store data
  push @{ $self->slices }, $slice_data;
  # clear result to signal need to recalculate result set
  $self->_clear_result;
}

=head2 result( )

Method to retrieve results from this bin; results are calculated and cached if the cache is not already populated. The return is an array reference whose data is given as [ total electrons, average momentum of these electrons, uncertainty of these electrons].

=cut

# Method to return results and if needed call _calc_result to get those results
has 'result' => (
  is => 'lazy',
  clearer => '_clear_result',
);

=head1 ACCESSOR METHODS

Simple accessor methods (combined set/get) are provided for 

=over 

=item * begin

=item * end

=back

 $bin->begin($new_position);
 $beginning_postition = $bin->begin();

Note that changing these values after storing data will not update the data. To state differently, the bin's data is independent of its position. 

=cut

has 'begin' => (
  is => 'rw',
  required => 1,
);

has 'end' => (
  is => 'rw',
  required => 1,
);

has 'slices' => (
  is => 'ro',
  default => sub { [] },
);

#-----------------------
# Private Methods
#-----------------------

# Method to compute final state of bin
sub _build_result {
  my $self = shift;

  # initialize holders
  my $total_num = 0;
  my $avg_momentum = 0;
  my $momentum_2 = 0;
  foreach my $slice (@{ $self->slices }) { # loop over all stored slices
    # get values in slice
    my $num = $slice->[0];
    my $v_avg = $slice->[1];
    # add values to holders
    $total_num += $num;
    $avg_momentum += $num * $v_avg;
    $momentum_2 += $num * ($v_avg ** 2);
  }

  # divide to find average momentum
  $avg_momentum = mass * $avg_momentum / $total_num;
  # calculate the momentum uncertainty in the bin
  my $momentum_uncertainty = sqrt( abs( (mass ** 2) * $momentum_2 / $total_num - ($avg_momentum ** 2))); #DeltaP^2 = eta

  return [$total_num, $avg_momentum, $momentum_uncertainty];
}

1;
