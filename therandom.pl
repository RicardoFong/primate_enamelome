#!/usr/bin/perl

=This program 
generates a random number from 0 to 10
using Crypt::Random package.
The output is an integer string
=cut

use strict;
use warnings;
use 5.010;

use Crypt::Random qw( makerandom_itv ); 

my $N = 10;
my $r = makerandom_itv ( Size => $N, Strength => 1, Uniform => 1, Lower => 0, Upper => 10 );
say $r;
