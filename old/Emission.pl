#!/usr/bin/env perl

use strict;
use warnings;

use PDL;
use PDL::GSL::INTEG;

my $input = shift;

my $qe = 1.6E-19;
my $me = 9.1E-31;
my $Elaser = $qe * 4.75;
my $WF = $qe * 4.25;
my $pi = 3.14159;
my $hbar = 6.63E-34 / (2 * $pi);

my $kV = sqrt( 2 * $me * $WF ) / $hbar;
my $vmax = sqrt( 2 * ($Elaser - $WF) / $me );

my $epsrel = 0;
my $epsabs = 1e-6;

my $norm = V_Integral($vmax);

print num($input) . "\n";

### Subs ###

sub num {
  my ($v) = @_;

  return V_Integral($v) / $norm;
}

sub Trans {
  my ($kz) = @_;
  my $sqrt = sqrt($kV ** 2 + $kz ** 2);

  my $numerator = 4 * $sqrt * $kz;
  my $denominator = ( $sqrt + $kz ) ** 2;

  return $numerator / $denominator;
}

sub Generate_Theta_Integrand {
  my ($v) = @_;

  my $sub = sub {
    my ($th) = @_;
    Trans($me * $v * cos($th) / $hbar) * sin($th);
  };

  return $sub;
}

sub Theta_Integral {
  my ($v) = @_;

  my $theta_integrand = Generate_Theta_Integrand($v);

  my $a = 0;
  my $b = $pi / 2;

  my ($res,$abserr,$ierr,$neval) = gslinteg_qng($theta_integrand,$a,$b,$epsrel,$epsabs);

  return $res;
}

sub V_Integral {
  my ($v_l) = @_;

  my $a = 0;

  my ($res,$abserr,$ierr,$neval) = gslinteg_qng(\&Theta_Integral,$a,$v_l,$epsrel,$epsabs);

  return $res;
}
