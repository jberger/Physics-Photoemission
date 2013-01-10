#!/usr/bin/env perl

use strict;
use warnings;

use PDL;
use PDL::GSL::INTEG;
use PDL::Fit::LM;

use Getopt::Long;
my $simple = 0;
my $output_filename = 'output';

my $tau = 100E-15; # HW1/eM temporal duration of the laser
my $num_electrons = 1E6; # total number of electrons
my $WF = 4.25; # Work fuction of the photocathode (Ta -> 4.25)
my $h_nu = 4.75; # Energy of the laser
my $DC_field = 1E6; # DC field strength
GetOptions(
  "simple" => \$simple,
  "filename=s" => \$output_filename,
  "tau=f" => \$tau,
  "num=i" => \$num_electrons,
  "WF=f" => \$WF,
  "hnu=f" => \$h_nu,
  "DC=f" => \$DC_field
);

use Term::ProgressBar;

#use Math::Libm 'erf';

my $verbose = 0;
my $update = 0;

use Bin;
my @bins;

### User supplied values ###
my $num_space_bins = 1000;
my $num_time_slices = 1000;
my $num_taus = 6; # number of taus to include

use constant {
  pi => 3.14159,
};
use constant {
  mass => 9.1E-31, # mass of the electron
  charge => 1.6E-19, # charge of the electron
  hbar => (6.63E-34) / (2 * pi),
};

### Calculated constants ###
my $end_time = $num_taus * $tau; # last time bin
my $Emax = charge * ($h_nu - $WF); # Energy of an electron with the full deltaE
my $acc = charge * $DC_field / mass; # acceleration in the gun

my $vmax = sqrt( 2 * $Emax / mass ); # maximum velocity of an electron
my $kV = sqrt( 2 * mass * $WF ) / hbar; # "momentum" of the barrier

print "Vmax = " . $vmax . "\n" if $verbose;

my $dmax = $vmax * $end_time + ( $acc * ( $end_time ** 2)) / 2;


my $conditions_text = <<END;
### Conditions ###
tau = $tau
WF = $WF
h_nu = $h_nu
DC_field = $DC_field

### Results ###
END

# Integration error limits
my $epsrel = 0;
my $epsabs = 1e-6;
# Norm of for forming the distribution function num($v_l)
my $norm = V_Integral(0, $vmax);

my $binwidth = $dmax / $num_space_bins;
my $slice_duration = $end_time / $num_time_slices;

### Create spatial bins ###
foreach my $i (0 .. ($num_space_bins - 1)) {
  push @bins, Bin->new( $binwidth * $i, $binwidth * ($i + 1) );
}

### Main Simulation ###
my $total = 0; # keep a total of electrons allocated

my $progress = Term::ProgressBar->new({
  name => 'Simulating',
  count => $num_time_slices - 1,
  ETA => 'linear',
});
my $num_before_slice = 0; # holder to keep the number to subtract off the erf() to get electrons in slice
foreach my $j (0 .. ($num_time_slices - 1)) { # loop over time slices

  # calculate time from start to slice generation
  my $elapsed_time = $j * $slice_duration;
  # calculate time from slice generation to end of simulation
  my $propagation_time = $end_time - $elapsed_time;

  # calculate number of electrons in slice
  my $num_in_slice = 
    $num_electrons * ( erf( ($elapsed_time - ($num_taus * $tau / 2)) / $tau ) + 1 ) / 2 - $num_before_slice;
  # add electrons in slice to total already photoemitted (num_in_slice)
  $num_before_slice += $num_in_slice;

  # simulate this slice
  propagate_slice($num_in_slice, $propagation_time, \@bins, \$total);

  # update progress
  $progress->update($j);
}
$progress->update($num_time_slices - 1);
print "\n";

print "Allocated a total of " . $total . " electrons (" . $total / $num_electrons * 100 . "\%)\n" if $verbose;

### Collect results, calculate AG model parameters and write to file ###

# open file for output
open(my $output, '>', $output_filename . '.csv');

#initialize some variables
my $av_x = 0;
my $av_x2 = 0;

my $max_electrons = 0;
my $x_max_electrons = 0;
my $v_max_electrons = 0;

# loop over all spatial bins
#TODO add progress bar here
my $process_progress = Term::ProgressBar->new({
  name => 'Processing',
  count => scalar @bins,
  ETA => 'linear',
});
my $iter = 0;
foreach my $bin (@bins) {
  # get results in bin (located at $x)
  my $x = $bin->end();
  my @result = $bin->result();

  my $electrons = $result[0];

  # find peak position
  if ($electrons > $max_electrons) {
    $max_electrons = $electrons;
    $x_max_electrons = $x;
    $v_max_electrons = $result[1] / mass;
  }

  # for finding sigma_z ...
  $av_x += $electrons * $x;
  $av_x2 += $electrons * ($x ** 2);

  # write output to file
  print $output join(',', $x, @result) . "\n";
  $process_progress->update($iter++);
}

$process_progress->update(scalar @bins);
print "\n";

# close file
close($output);

# calculate sigma_z
my $sigma = ($av_x2 / $num_electrons) - (($av_x / $num_electrons) ** 2);

# call fit_momentum to calculate eta and gamma
my ($eta, $gamma, $gamma_offset, $sigma_coeff) = fit_routine($x_max_electrons, $sigma);

# print results
open(my $result_handle, '>', $output_filename . '.rslt');

my $result_text = <<END;
Electron pulse peak position = $x_max_electrons
Electron pulse peak velocity = $v_max_electrons
sigma_z = $sigma
eta_z = $eta
gamma_z = $gamma

gamma_offset = $gamma_offset
Fit Equation for sigma:
$sigma_coeff * exp(-(x - $x_max_electrons)^2 / (2 * $sigma))
Fit Equation for gamma:
($gamma / $sigma) * x + $gamma_offset
END

print $result_text;

print $result_handle $conditions_text;
print $result_handle $result_text;

close($result_handle);

### Subroutines ###

sub propagate_slice {
  my ($num_in_slice, $propagation_time, $bin_array_ref, $total_ref) = @_;

  foreach my $bin (@{ $bin_array_ref }) { # loop over all bins
    my $begin = $bin->begin(); # get bins start ...
    my $end = $bin->end(); # ... and end locations

    # calculate the initial velocities that reach the beginning and end of the spatial bin
    my $vi_begin = $begin / $propagation_time - $acc * $propagation_time / 2;
    my $vi_end = $end / $propagation_time - $acc * $propagation_time / 2;

    # calculate the final velocities that these electrons have as they reach these locations
    my $vf_begin = sqrt( ($vi_begin ** 2) + (2 * $acc * $begin) );
    my $vf_end = sqrt( ($vi_end ** 2) + (2 * $acc * $end) );

    # calculate the number of electrons from the slice that land in the spatial bin
    my $num = number_by_vel($vi_begin, $vi_end) * $num_in_slice;
    # optionally add this amount to the total electrons allocated
    ${$total_ref} += $num if (defined $total_ref);
    #print "In bin " . $begin . "-" . $end . " are " . $num . " electons\n" if $verbose;

    # store the calculated slice results in the bin
    $bin->add_slice( $num, $vf_begin, $vf_end);
  }

}

sub number_by_vel {
  my ($v_1, $v_2) = @_;

  # set minimum lower velocity to 0
  $v_1 = ($v_1 > 0) ? $v_1 : 0;
  # set maximum upper velocity to vmax
  $v_2 = ($v_2 < $vmax) ? $v_2 : $vmax;

  my $return = 0;
  unless ($v_1 > $vmax or $v_2 < 0) { # as long as the velocities are in the proper range
    # calculate fraction of electrons from slice in the bin
    if ($simple) {
      $return = (($v_2 / $vmax) ** 5) - (($v_1 / $vmax) ** 5);
    } else {
      $return = num($v_1, $v_2);
    }
  }

  return $return;
}

### Numerical Integration Subroutines ###

# uses the velocity (and theta) integration of the transmission function but normalize 
# making it an effective distribution function
sub num {
  my ($v_l, $v_u) = @_;

  return V_Integral($v_l, $v_u) / $norm;
}

# 2D transmission over a 1D barrier in terms of external k_perpendicular
sub Trans {
  my ($kz) = @_;
  my $sqrt = sqrt($kV ** 2 + $kz ** 2);

  my $numerator = 4 * $sqrt * $kz;
  my $denominator = ( $sqrt + $kz ) ** 2;

  return $numerator / $denominator;
}

# sub which returns the integrand for the theta integration for a particular velocity
sub Generate_Theta_Integrand {
  my ($v) = @_;

  my $sub = sub {
    my ($th) = @_;
    Trans(mass * $v * cos($th) / hbar) * sin($th);
  };

  return $sub;
}

# sub which both does the theta integration and acts as the integrand for the velocity integration
sub Theta_Integral {
  my ($v) = @_;

  my $theta_integrand = Generate_Theta_Integrand($v);

  my $a = 0;
  my $b = pi / 2;

  my ($res,$abserr,$ierr,$neval) = gslinteg_qng($theta_integrand,$a,$b,$epsrel,$epsabs);

  return $res;
}

# sub which does integration over v (and theta) as a function of the lower and upper limit velocity
sub V_Integral {
  my ($v_l, $v_u) = @_;

  my ($res,$abserr,$ierr,$neval) = gslinteg_qng(\&Theta_Integral,$v_l,$v_u,$epsrel,$epsabs);

  return $res;
}

#sub for finding eta and gamma
sub fit_routine {
  my ($x0, $sigma) = @_;
  my $x_pm = sqrt( 2 * $sigma );

  # get only the bins in the desired range
  my @fit_bins = grep { ($_->begin() >= ($x0 - $x_pm)) and ($_->end() <= ($x0 + $x_pm)) } @bins;
  my $x_pdl = pdl( map { $_->end() } @fit_bins );
  my $N_pdl = pdl( map { ($_->result())[0] } @fit_bins );
  my $mom_pdl = pdl( map { ($_->result())[1] } @fit_bins );
  my $mom_uncert_2_pdl = pdl( map { (($_->result())[2]) ** 2 } @fit_bins );

  # coderef for linear fitting
  my $fit_linear = sub {
    my ($x,$par,$ym,$dyda) = @_;

    my ($m,$b) = map { $par->slice("($_)") } (0..1);

    $ym .= $m * $x + $b;

    my (@dy) = map {$dyda -> slice(",($_)") } (0..1);

    $dy[0] .= $x;
    $dy[1] .= 1;
  };

  my ($mom_fit, $mom_fit_params, $mom_fit_cf, $mom_fit_iter) = lmfit $x_pdl, $mom_pdl, 
    ($N_pdl / maximum($N_pdl)) ** (-0.5), # Use a linear weight (rather than 1/w^2)
    $fit_linear, pdl(0,0)
  ;
  my $gamma = $sigma * $mom_fit_params->index(0);
  my $gamma_offset = $mom_fit_params->index(1);

  my $eta = sum( $N_pdl * $mom_uncert_2_pdl ) / sum( $N_pdl );

  my $sigma_coeff = sum($N_pdl) * $binwidth / sqrt( 2 * pi * $sigma);

  return ($eta, $gamma, $gamma_offset, $sigma_coeff);
}
