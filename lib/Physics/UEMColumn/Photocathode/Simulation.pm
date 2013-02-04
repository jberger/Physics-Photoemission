package Physics::UEMColumn::Photocathode::Simulation;

=head1 NAME

Physics::UEMColumn::Photocathode::Simulation - Subclass of Physics::UEMColumn::Photocathode using Physics::Photoemission 

=head1 SYNOPSIS

  use strict;
  use warnings;

  use Physics::UEMColumn::Photocathode::Simulation 'SimPhotocathode' => {
    num_space_bins => 100,
    num_time_slices => 100,
  };

  my $photocathode = SimPhotocathode->new(
    work_function => '4.25 eV',
  );

  # note that $photocathode must have some access to an appropriate Column object

  my $pulse = $photocathode->generate_pulse( 1e8 );

=head1 DESCRIPTION

This class is a subclass of L<Physics::UEMColumn::Photocathode> which uses L<Physics::Photoemission> to generate the longitudinal initial conditions of the generated pulse.

Since it is as subclass, it implements all the attributes and methods of that class, with the following changes.

=cut

use Moo;
BEGIN{ extends 'Physics::UEMColumn::Photocathode' }

use Physics::Photoemission;
use Physics::Photoemission::Constants qw/charge/;

use Carp;

=head1 PACKAGE VARIABLES

=over

=item C<$Emission_Time>

The simulation allows a small amount of time to run off in generating the pulse. This package variable contains the amount of run off time for the last generation. One might want to localize this variable if it is to be used to correct the timing of a L<Physics::UEMColumn> simulation. The value is simply calculated by multiplying C<num_taus> and C<tau> from L<Physics::Photoemission>.

=cut

our $Emission_Time = 0;

=item C<$Anon_Class_Num>

This number is incremented for each anonymous subclass created (see L</IMPORTING>).

=back

=cut

our $Anon_Class_Num = 0;

=head1 ATTRIBUTES

=over

=item C<simulation>

Holder for the L<Physics::Photoemission> object. It is generated lazily when needed. Its C<tau>, C<work_function>, C<photon_energy> and C<field> are loaded from the information in the surrounding simulation. Additional overrides (see below) are merged with higher priority.

=cut

has 'simulation' => (
  is => 'lazy',
);
sub _build_simulation {
  my $self = shift;

  my $column = $self->column;
  my $laser = $column->laser;

  # sim expects energies in eV
  my $sim = Physics::Photoemission->new(
    tau           => $laser->duration,
    work_function => $self->work_function / charge,
    photon_energy => $laser->energy / charge,
    field         => $column->accelerator->field,

    %{ $self->simulation_overrides },
  );
}

=item C<simulation_overrides>

A hash reference merged with the automatic attributes to the C<simulation> object (above). These attributes are merged with higher priority than the automatic ones should you want to override those too.

=back

=cut

has 'simulation_overrides' => (
  is => 'lazy',
);
sub _build_simulation_overrides { +{} }

=head1 METHODS

=over

=item C<generate_pulse>

Overrides the method in L<Physics::UEMColumn::Photocathode>, though that method is still called in order to get the transverse conditions. This method is called by the L<Physics::UEMColumn> simulation to generate the pulse object (if necessary).

=back

=cut

sub generate_pulse {
  my ($self, $num) = @_;

  my $pulse = $self->SUPER::generate_pulse($num);

  print STDERR "Simulating pulse emission process\n";

  my $sim = $self->simulation;
  $sim->num_electrons( $num );

  $sim->simulate;
  $sim->process;
  my $result = $sim->result->{values};

  # The simulation of emission elapses some time which cannot otherwise be relayed to the column simulation.
  # Therefore this offset time is stored in a package variable for inspection if desired.
  # Future plans might include returning this in the pulse object itself.
  $Emission_Time = $sim->tau * $sim->num_taus;
  printf STDERR "Note that the pulse was created at an effective time %es beyond the simulation begin time.\n", $Emission_Time;

  $pulse->sigma_z( $result->{sigma} );
  $pulse->eta_z( $result->{eta} );
  $pulse->gamma_z( $result->{gamma} );

  $pulse->velocity( $result->{v_peak} );
  $pulse->location( $pulse->location + $result->{x_peak} );

  return $pulse;
}

=head1 IMPORTING

As with L<Physics::UEMColumn|Physics::UEMColumn/"IMPORTING"> this class can create an alias for easier use in your scripts. The calling convention is slightly different however since there is only one class (and its a subclass of one likely to be aliased). Any string passed as an import argument will be aliased to the name C<Physics::UEMColumn::Photocathode::Simulation>. In this documentation and examples the alias is by convention C<SimPhotocathode>.

 use Physics::UEMColumn alias => ':standard';
 use Physics::UEMColumn::Photocathode::Simulation 'SimPhotocathode';
 my $pulse = SimPhotocathode->new( ... );

As a further convenience, if a second arguement (a hash reference) is passed, that will be used as the default value of C<simulation_overrides>. This is accomplished by creating an "anonymous" subclass of this package and storing the name in the created alias.

 use Physics::UEMColumn alias => ':standard';
 use Physics::UEMColumn::Photocathode::Simulation 'SimPhotocathode' => {
   num_space_bins  => 100,
   num_time_slices => 100,
 };
 my $pulse = SimPhotocathode->new( ... );
 
 print ref $pulse; # Physics::UEMColumn::Photocathode::Simulation::Anon0

Finally for symmetry, an initial key "alias" is silently ignored in any case.

 use Physics::UEMColumn alias => ':standard';
 use Physics::UEMColumn::Photocathode::Simulation alias => 'SimPhotocathode';

=cut

# hic sunt dracones
sub import {
  my $class = shift;
  my $name = shift || return;

  # P::U aliasing uses an "alias" key, so this support is here for symmetry
  if ( $name eq 'alias' && $_[0] && ! ref $_[0] ) {
    $name = shift;
  }

  my $caller = caller;
  my $package = __PACKAGE__;

  # generate an "anonymous" subclass with new default simulation_overrides
  if ( my $overrides = shift ) {
    croak "Second import argument must be a hash reference"
      unless ref $overrides eq 'HASH';

    my $shallow_copy = { %$overrides };

    $package .= '::Anon' . $Anon_Class_Num++;

    no strict 'refs';
    @{ $package . '::ISA' } = __PACKAGE__;
    *{ $package . '::_build_simulation_overrides' } = sub { $shallow_copy };
  }

  no strict 'refs';
  *{ $caller . '::' . $name } = eval "sub () { '$package' }";
}

1;

=head1 SOURCE REPOSITORY

L<http://github.com/jberger/Physics-Photoemission>

=head1 AUTHOR

Joel Berger, E<lt>joel.a.berger@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Joel Berger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.


