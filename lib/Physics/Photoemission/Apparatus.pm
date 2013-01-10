package Physics::Photoemission::Apparatus;

use Moo::Role;
use Physics::Photoemission::Constants ':all';

=head1 NAME

Physics::Photoemission::Apparatus;

=head1 SYNOPSIS

 package Physics::Photoemission;
 use Moo;
 with 'Physics::Photoemission::Apparatus';

=head1 DESCRIPTION

This package implements the apparatus role used in the C<Physics::Photoemission> module/simulation. The role really is just a holder for experimental parameters of the simulation. This defines many physical quantities as attributes most are read-only except the few that are not reused anywhere or calculated from other quantities.

N.B. All units are MKS except energies which are in eV for convenience.

=head1 METHODS

=head2 new( [ $opts_hashref | %opts ] )

Constructor method which creates a new apparatus. The constructor can take a hash or hashref of parameters with the same name as the "get/set accessors" below.

=cut

# HW1/eM temporal duration of the laser
has 'tau' => (
  is => 'rw',
  default => sub { 10e-12 },
);

has 'num_electrons' => (
  is => 'rw',
  default => sub { 1e6 },
);

# Work fuction of the photocathode (Ta -> 4.25)
has 'WF' => (
  is => 'ro',
  default => sub { 4.25 },
);

# Photon energy of the laser
has 'photon_energy' => (
  is => 'ro',
  default => sub { 4.75 },
);

has 'DC_field' => (
  is => 'ro',
  default => sub { 1e6 },
);

# Energy of an electron with the full deltaE
has 'Emax' => ( is => 'lazy' );
sub _build_Emax {
  my $self = shift;
  charge * ( $self->photon_energy - $self->WF );
}

# acceleration in the gun
has 'acc' => ( is => 'lazy' );
sub _build_acc {
  my $self = shift;
  charge * $self->DC_field / mass;
}

# maximum velocity of an emitted electron
has 'vmax' => ( is => 'lazy' );
sub _build_vmax {
  my $self = shift;
  sqrt( 2 * $self->Emax / mass );
}

# "momentum" of the barrier
has 'kV' => ( is => 'lazy' );
sub _build_kV {
  my $self = shift;
  sqrt( 2 * mass * $self->WF ) / hbar;
}

=head1 GET/SET ACCESSOR METHODS

There are "get/set accessors" for the following "user-defined" quantities:

=over

=item * tau			HW1/eM laser duration (s)

=item * num_electrons		Total number of electrons

=item * WF			Work function of the photocathode (eV)

=item * photon_energy		Photon energy (eV)

=item * DC_field		Acceleration field strength (V/m)

=back

 $apparatus->num_electrons( 1e6 );
 my $num_electrons = $apparatus->num_electrons();

These names may also be used as the keys of the hash which is passed to the constructor C<new()>.

 my $apparatus = Physics::Photoemission::Apparatus( photon_energy => 5.0 );

=head1 GET-ONLY ACCESSOR METHODS

There are "get-only accessors" for the following derived quantities. These quantities are derived by the C<init()> method (which is implicitly called when accessing one of these quantities when it is undefined and is rerun when one of its dependencies is changed). The user probably does not need to call this method directly.

=over

=item * Emax		Maximum energy a photoemitted electron may have (eV)

=item * vmax		The velocity of an electron with C<Emax> energy (m/s)

=item * acc		Acceleration of an electron due to the electric field (m/s^2)

=item * kV		"Momentum" related to the work function (kg m/s)

=back

=cut



1;

