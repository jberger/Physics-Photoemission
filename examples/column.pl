#!/usr/bin/env perl

use strict;
use warnings;

use Physics::UEMColumn alias => ':standard';
use Physics::UEMColumn::Auxiliary ':materials';
use Physics::UEMColumn::Photocathode::Simulation 'SimPhotocathode' => {
  num_space_bins  => 100,
  num_time_slices => 100,
};

my $laser = Laser->new(
  width    => '1 mm',
  duration => '1 ps',
  energy   => '4.75 eV',
);

my $acc = DCAccelerator->new(
  length  => '20 mm',
  voltage => '20 kilovolts',
);

my $column = Column->new(
  length       => '400 mm', 
  laser        => $laser,
  accelerator  => $acc,
  photocathode => SimPhotocathode->new(Ta),
);

my $sim = Physics::UEMColumn->new(
  column => $column,
  number => 1,
);

$sim->propagate;

