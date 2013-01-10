package Physics::Photoemission;

use Moo;
BEGIN { with 'Physics::Photoemission::Apparatus' }

use PDL;
use PDL::GSL::INTEG;
use PDL::Fit::LM;

use Term::ProgressBar;

use Physics::Photoemission::Constants ':all';
use Physics::Photoemission::Bin;

=head1 NAME

Physics::Photoemission

=head1 SYNOPSIS

 use Physics::Photoemission;
 my $model = Physics::Photoemission->new;
 $model->simulate;
 $model->process;
 print $model->report;

or

 perl -MPhysics::Photoemission -E 'run_sim'

same as above.

=head1 DESCRIPTION

This module models the initial longitudinal distribution of electrons photogenerated by a gaussian laser pulse. These parameters are of the form to be used as the initial conditions of the AG model first presented by Michalik and Sipe (J. App. Phys. [vol] [page]).

=head1 EXPORTED FUNCTION

=head2 run_sim( [ %opts | $opts_hashref ] )

Exported as a shortcut to run a basic simulation from the command line. See C<SYNOPSIS> for the description of the steps taken. Any parameters passed to C<run_sim()> are passed directly to the constructor C<new()>.

=cut

## Allow 'run_sim' to be used from the command line as
## `perl -MPhysics::Photoemission -E 'run_sim()'`
## N.B. parameters given to run_sim() are passed directly to new()
use parent 'Exporter';
our @EXPORT = qw/run_sim/;

sub run_sim {
  my $model = __PACKAGE__->new(@_);
  $model->simulate;
  $model->process;
  print $model->report;
}

=head1 METHODS

=head2 new( [ %opts | $opts_hashref ] )

=cut

has 'simple' => (
  is => 'rw',
  default => sub { 0 },
);

has 'verbose' => (
  is => 'rw',
  default => sub { 0 },
);

has 'output_filename' => (
  is => 'rw',
  default => sub { 'output.csv' },
);

has 'output_filehandle' => ( is => 'lazy' );
sub _build_output_filehandle {
  my $self = shift;
  my $file = $self->output_filename;
  open my $fh, '>', $file or die "Cannot open output filehandle";
  print $fh qq{"bin position","number of electrons","bin average momentum","bin momentum uncertainty"\n};
  return $fh;
}

has 'num_space_bins' => (
  is => 'rw',
  default => sub { 1000 },
);

has 'num_time_slices' => (
  is => 'rw',
  default => sub { 1000 },
);

has 'num_taus' => (
  is => 'rw',
  default => sub { 6 },
);

has 'epsrel' => (
  is => 'rw',
  default => sub { 0 },
);

has 'epsabs' => (
  is => 'rw',
  default => sub { 1e-6 },
);

has 'end_time' => ( is => 'lazy' );
sub _build_end_time { 
  my $self = shift;
  $self->num_taus * $self->tau;
}

has 'dmax' => ( is => 'lazy' );
sub _build_dmax {
  my $self = shift;
  $self->vmax * $self->end_time + ( $self->acc * ( $self->end_time ** 2)) / 2;
}

has 'binwidth' => ( is => 'lazy' );
sub _build_binwidth {
  my $self = shift;
  $self->dmax / $self->num_space_bins;
}

has 'slice_duration' => ( is => 'lazy' );
sub _build_slice_duration {
  my $self = shift;
  $self->end_time / $self->num_time_slices;
}

# vnorm for forming the distribution function num($v_l)
has 'vnorm' => ( is => 'lazy' );
sub _build_vnorm {
  my $self = shift;
  $self->V_Integral(0, $self->vmax);
}

has 'total' => (
  is => 'rw',
  default => sub { 0 },
);

has 'bins' => ( is => 'lazy' );
sub _build_bins {
  my $self = shift;
  my $binwidth = $self->binwidth;
  my $bins = [];
  # Create spatial bins
  foreach my $i ( 0 .. ($self->num_space_bins - 1) ) {
    push @$bins, Physics::Photoemission::Bin->new( 
        $binwidth * $i, $binwidth * ($i + 1) 
    );
  }
  return $bins;
}

has 'result' => ( is => 'lazy' );
sub _build_result {
  my $self = shift;
  my $result = {
    text => {
      conditions => undef,
      result     => undef,
    },
    values => {},
  };

  my $tau = $self->tau;
  my $WF = $self->WF;
  my $photon_energy = $self->photon_energy;
  my $DC_field = $self->DC_field;
  $result->{text}{conditions} = <<END;
### Conditions ###
tau = $tau
WF = $WF
photon_energy = $photon_energy
DC_field = $DC_field
END

  return $result;
}

=head2 simulate( )

=cut

### Main Simulation ###
sub simulate {
  my $self = shift;

  my $total = 0; # keep a total of electrons allocated
  my $num_time_slices = $self->num_time_slices - 1; # -1 accounts for perl index at zero
  my $slice_duration = $self->slice_duration;
  my $num_taus = $self->num_taus;
  my $num_electrons = $self->num_electrons;
  my $tau = $self->tau;

  my $progress = Term::ProgressBar->new({
    name => 'Simulating',
    count => $num_time_slices,
    ETA => 'linear',
  });
  my $num_before_slice = 0; # holder to keep the number to subtract off the erf() to get electrons in slice
  foreach my $j (0 .. $num_time_slices) { # loop over time slices

    # calculate time from start to slice generation
    my $elapsed_time = $j * $slice_duration;
    # calculate time from slice generation to end of simulation
    my $propagation_time = $self->end_time - $elapsed_time;

    # calculate number of electrons in slice
    my $num_in_slice = 
      $num_electrons * ( erf( ($elapsed_time - ($num_taus * $tau / 2)) / $tau ) + 1 ) / 2 - $num_before_slice;
    # add electrons in slice to total already photoemitted (num_in_slice)
    $num_before_slice += $num_in_slice;

    # simulate this slice
    $self->propagate_slice($num_in_slice, $propagation_time);

    # update progress
    $progress->update($j);
  }
  $progress->update($num_time_slices);
  print "\n";

  if ( $self->verbose ) {
    my $total = $self->total;
    my $percent = $total / $num_electrons * 100;
    print "Allocated a total of $total electrons ($percent\%)\n";
  }

  return $self->result;
}

=head2 process( )

=cut

### Collect results, calculate AG model parameters and write to file ###
sub process {
  my $self = shift;

  my $output = $self->output_filehandle;

  #initialize some variables
  my $av_x = 0;
  my $av_x2 = 0;

  my $num_electrons = $self->num_electrons;
  my $max_electrons = 0;
  my $x_peak = 0;
  my $v_peak = 0;

  my $bins = $self->bins;

  # loop over all spatial bins
  my $process_progress = Term::ProgressBar->new({
    name => 'Processing',
    count => scalar @$bins,
    ETA => 'linear',
  });
  my $iter = 0;
  foreach my $bin (@$bins) {
    # get results in bin (located at $x)
    my $x = $bin->end;
    my $result = $bin->result;

    my $electrons = $result->[0];

    # find peak position
    if ($electrons > $max_electrons) {
      $max_electrons = $electrons;
      $x_peak = $x;
      $v_peak = $result->[1] / mass;
    }

    # for finding sigma_z ...
    $av_x += $electrons * $x;
    $av_x2 += $electrons * ($x ** 2);

    # write output to file
    print $output join(',', $x, @$result) . "\n";
    $process_progress->update($iter++);
  }

  $process_progress->update(scalar @$bins);
  print "\n";

  # calculate sigma_z
  my $sigma = ($av_x2 / $num_electrons) - (($av_x / $num_electrons) ** 2);

  # call fit_momentum to calculate eta and gamma
  my ($eta, $gamma, $gamma_offset, $sigma_coeff) = $self->fit_routine($x_peak, $sigma);

  my $result = $self->result;
  $result->{values} = {
    sigma        => $sigma,
    sigma_coeff  => $sigma_coeff,
    eta          => $eta,
    gamma        => $gamma,
    gamma_offset => $gamma_offset,
    x_peak       => $x_peak,
    v_peak       => $v_peak,
  };

  return $result;
}

=head2 report( )

=cut

# print results
sub report {
  my $self = shift;
  my $result = $self->result;
  my %values = %{ $result->{values} };

  $result->{text}{result} = <<END;
### Results ###
Electron pulse peak position = $values{x_peak}
Electron pulse peak velocity = $values{v_peak}
sigma_z = $values{sigma}
eta_z = $values{eta}
gamma_z = $values{gamma}

gamma_offset = $values{gamma_offset}

Fit Equation for sigma:
$values{sigma_coeff} * exp(-(x - $values{x_peak})^2 / (2 * $values{sigma}))
Fit Equation for gamma:
($values{gamma} / $values{sigma}) * x + $values{gamma_offset}
END

  return $result->{text}{conditions} . $result->{text}{result};
}

### Functions ###

sub propagate_slice {
  my $self = shift;
  my ($num_in_slice, $propagation_time) = @_;

  foreach my $bin (@{ $self->bins }) { # loop over all bins
    my $begin = $bin->begin; # get bins start ...
    my $end = $bin->end; # ... and end locations
    my $acc = $self->acc;

    # calculate the initial velocities that reach the beginning and end of the spatial bin
    my $vi_begin = $begin / $propagation_time - $acc * $propagation_time / 2;
    my $vi_end = $end / $propagation_time - $acc * $propagation_time / 2;

    # calculate the final velocities that these electrons have as they reach these locations
    my $vf_begin = sqrt( ($vi_begin ** 2) + (2 * $acc * $begin) );
    my $vf_end = sqrt( ($vi_end ** 2) + (2 * $acc * $end) );

    # calculate the number of electrons from the slice that land in the spatial bin
    my $num = $self->number_by_vel($vi_begin, $vi_end) * $num_in_slice;
    # add this amount to the total electrons allocated
    $self->total( $self->total + $num );
    #print "In bin " . $begin . "-" . $end . " are " . $num . " electons\n" if $verbose;

    # store the calculated slice results in the bin object
    $bin->add_slice( $num, $vf_begin, $vf_end);
  }

}

sub number_by_vel {
  my $self = shift;
  my ($v_1, $v_2) = @_;
  my $vmax = $self->vmax;

  # set minimum lower velocity to 0
  $v_1 = ($v_1 > 0) ? $v_1 : 0;
  # set maximum upper velocity to vmax
  $v_2 = ($v_2 < $vmax) ? $v_2 : $vmax;

  my $return = 0;
  unless ($v_1 > $vmax or $v_2 < 0) { # as long as the velocities are in the proper range
    # calculate fraction of electrons from slice in the bin
    $return = 
      $self->simple
      ? (($v_2 / $vmax) ** 5) - (($v_1 / $vmax) ** 5)
      : $self->num($v_1, $v_2);
  }

  return $return;
}

### Numerical Integration Subroutines ###

# uses the velocity (and theta) integration of the transmission function but vnormalize 
# making it an effective distribution function
sub num {
  my $self = shift;
  my ($v_l, $v_u) = @_;

  return $self->V_Integral($v_l, $v_u) / $self->vnorm;
}

# 2D transmission over a 1D barrier in terms of external k_perpendicular
sub Trans {
  my $self = shift;

  my ($kz) = @_;
  my $sqrt = sqrt($self->kV ** 2 + $kz ** 2);

  my $numerator = 4 * $sqrt * $kz;
  my $denominator = ( $sqrt + $kz ) ** 2;

  return $numerator / $denominator;
}

# sub which does integration over v (and theta) as a function of the lower and upper limit velocity
sub V_Integral {
  my $self = shift;
  my ($vel_lower, $vel_upper) = @_;
  my $epsrel = $self->epsrel;
  my $epsabs = $self->epsabs;

  # anonymous sub which both does the theta integration and acts as the integrand for the velocity integration
  my $theta_integral = sub {
    my ($v) = @_;

    # anonymous sub which returns the integrand for the theta integration for a particular velocity
    my $theta_integrand = sub {
      my ($th) = @_;
      $self->Trans(mass * $v * cos($th) / hbar) * sin($th);
    };
    # end sub

    my $angle_lower = 0;
    my $angle_upper = pi / 2;

    # do the theta integration
    my ($res,$abserr,$ierr,$neval) = gslinteg_qng(
      $theta_integrand,$angle_lower,$angle_upper,$epsrel,$epsabs
    );

    return $res;
  };
  #end sub

  # do the velocity integration over the theta integral
  my ($res,$abserr,$ierr,$neval) = gslinteg_qng(
    $theta_integral,$vel_lower,$vel_upper,$epsrel,$epsabs
  );

  return $res;
}

#sub for finding eta and gamma
sub fit_routine {
  my $self = shift;
  my ($x0, $sigma) = @_;
  my $x_pm = sqrt( 2 * $sigma );

  # get only the bins in the desired range
  my @fit_bins = grep { ($_->begin >= ($x0 - $x_pm)) and ($_->end <= ($x0 + $x_pm)) } @{ $self->bins };
  my $x_pdl = pdl( map { $_->end } @fit_bins );
  my $N_pdl = pdl( map { ($_->result)[0] } @fit_bins );
  my $mom_pdl = pdl( map { ($_->result)[1] } @fit_bins );
  my $mom_uncert_2_pdl = pdl( map { (($_->result)[2]) ** 2 } @fit_bins );

  # coderef for linear fitting
  my $fit_linear = sub {
    my ($x,$par,$ym,$dyda) = @_;

    my ($m,$b) = map { $par->slice("($_)") } (0..1);

    $ym .= $m * $x + $b;

    my (@dy) = map { $dyda->slice(",($_)") } (0..1);

    $dy[0] .= $x;
    $dy[1] .= 1;
  };

  my ($mom_fit, $mom_fit_params, $mom_fit_cf, $mom_fit_iter) = lmfit(
    $x_pdl, 
    $mom_pdl, 
    ($N_pdl / maximum($N_pdl)) ** (-0.5), # Use a linear weight (rather than 1/w^2)
    $fit_linear,
    pdl(0,0)
  );
  my $gamma = $sigma * $mom_fit_params->index(0);
  my $gamma_offset = $mom_fit_params->index(1);

  my $eta = sum( $N_pdl * $mom_uncert_2_pdl ) / sum( $N_pdl );

  my $sigma_coeff = sum($N_pdl) * $self->binwidth / sqrt( 2 * pi * $sigma);

  return ($eta, $gamma, $gamma_offset, $sigma_coeff);
}

1;
