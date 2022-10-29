
use FileHandle;
use strict;


my $line;

my @pharm_map = 
  (
   ['n_runs', 5],
   ['n_iterations', 25000],
   ['do_rigid_fit',  'no'],
   ['do_rigid_fit',  'no'],
   ['guess_ga_parameters', 'no'],
  );

generate_conf ('pharm_search.conf', \@pharm_map);

my @rigid_map = 
  (
   ['log_level', 1],
   ['add_hydrogens', 'no'],
   ['solvate', 'no'],
);

generate_conf ('rigid_search.conf', \@rigid_map);

sub generate_conf {
  my ($file, $map) = @_;
  
  my $in_fh = new FileHandle 'superposition.conf', 'r';
  my $out_fh = new FileHandle $file, 'w';

  die unless $in_fh && $out_fh;

  while ($line = <$in_fh>) {
    $line =~ s/[\r\n]+$//;
    
    foreach my $kv (@$map) {
      my $field = $kv->[0];
      my $value = $kv->[1];
      my $find = "^$field";
      $line = "$field = $value" if $line =~ /^($field)/;
    }
    
    print $out_fh "$line\n";
  }
  
  close $in_fh;
  close $out_fh;
}
