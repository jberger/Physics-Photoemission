package Bin;

use strict;
use warnings;

use constant {
  mass => 9.1E-31
};

sub new {
  my $class = shift;
  my ($begin, $end) = @_;

  my $self = {
    slices => [], # array to hold slice data (as ordered hashrefs)
    begin => undef, # location of front of bin
    end => undef, # location of end of bin
    result => undef, # holder for array of results, or undef for result needing to be (re)calculated
  };

  bless($self, $class);

  # set locations
  $self->begin($begin);
  $self->end($end);

  return $self;
}

### Accessor Methods ###

sub begin {
  my $self = shift;

  $self->{'begin'} = shift if @_;

  return $self->{'begin'};
}

sub end {
  my $self = shift;

  $self->{'end'} = shift if @_;

  return $self->{'end'};
}

### Methods ###

# Method to store slice data in 'slices' array 
sub add_slice {
  my $self = shift;
  my ($num, $v_begin, $v_end) = @_;

  # data to store for each slice
  my %slice_data = (
    num => $num,
    v_avg => ($v_begin + $v_end) / 2,
    #v_begin => $v_begin,
    #v_end => $v_end,
  );

  # store data
  push @{ $self->{'slices'} }, \%slice_data;
  # set $self->{'result'} to undef to signal need to recalculate result set
  $self->{'result'} = undef;
}

# Method to return results and if needed call _calc_result to get those results
sub result {
  my $self = shift;

  unless (defined $self->{'result'}) {
    $self->_calc_result();
  }

  return @{$self->{'result'}};
}

### Private Methods ###

# Method to compute final state of bin
sub _calc_result {
  my $self = shift;

  # initialize holders
  my $total_num = 0;
  my $avg_momentum = 0;
  my $momentum_2 = 0;
  foreach my $slice (@{ $self->{'slices'} }) { # loop over all stored slices
    # get values in slice
    my $num = $slice->{'num'};
    my $v_avg = $slice->{'v_avg'};
    # add values to holders
    $total_num += $num;
    $avg_momentum += $num * $v_avg;
    $momentum_2 += $num * ($v_avg ** 2);
  }

  # divide to find average momentum
  $avg_momentum = mass * $avg_momentum / $total_num;
  # calculate the momentum uncertainty in the bin
  my $momentum_uncertainty = sqrt( abs( (mass ** 2) * $momentum_2 / $total_num - ($avg_momentum ** 2))); #DeltaP^2 = eta

  $self->{'result'} = [$total_num, $avg_momentum, $momentum_uncertainty];
}

1;
