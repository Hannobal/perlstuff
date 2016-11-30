package testmodule;

use 5.014002;
use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw($testvar module_testsub);

our $VERSION = '0.01';

our $testvar=15;
sub module_testsub {
  print $_[0],"\n";
}

# Preloaded methods go here.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

testmodule - Perl extension for blah blah blah

=head1 SYNOPSIS

  use testmodule;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for testmodule, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Hanno Dietrich, E<lt>dietrich@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Hanno Dietrich

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
