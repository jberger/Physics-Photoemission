package Physics::Photoemission::Constants;

use strict;
use warnings;

use parent 'Exporter';
our @EXPORT_OK = ( qw/
  pi
  mass
  charge
  hbar
/ );
our %EXPORT_TAGS = (
  all => \@EXPORT_OK
);

use constant {
  pi => 3.14159,
};
use constant {
  mass => 9.1E-31, # mass of the electron
  charge => 1.6E-19, # charge of the electron
  hbar => (6.63E-34) / (2 * pi),
};

1;

__END__
__POD__

=head1 NAME

Physics::Photoemission::Contants

=head1 SYNOPSIS

 use Physics::Photoemission::Constants ':all';
 say pi;

=head1 DESCRIPTION

A simple package which holds physical contants for use in C<Physics::Photoemission> and its related packages. Export any constants needed as described by C<Exporter> or export all with the tag ':all'.

=head1 CONSTANTS

=over
 
=item * pi		circle constant

=item * mass		mass of a free electron

=item * charge		charge of an electron

=item * hbar		planck's constant over 2*pi

=back

=cut
