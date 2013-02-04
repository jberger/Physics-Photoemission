package Physics::UEMColumn::Photocathode::Simulation;

use Moo;
BEGIN{ extends 'Physics::UEMColumn::Photocathode' }

use Physics::Photoemission;
use Physics::Photoemission::Constants qw/charge/;

use Carp;

our $Emission_Time = 0;
our $Anon_Class_Num = 0;

has 'simulation_overrides' => (
  is => 'lazy',
);
sub _build_simulation_overrides { +{} }

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

# hic sunt dracones
sub import {
  my $class = shift;
  my $name = shift || return;

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

