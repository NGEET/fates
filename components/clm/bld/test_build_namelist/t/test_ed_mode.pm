package test_ed_mode;

# Unit tests for function: ed_mode

use Data::Dumper;
use Test::More;
use Test::Exception;

use parent qw(Test::Class);

#-------------------------------------------------------------------------------
#
# Common test fixture for all tests:
#
#-------------------------------------------------------------------------------
sub startup : Test(startup => 4) {
  my $self = shift;
  # provide common fixture for all tests, only created once at the
  # start of the tests.

  $self->{cfg} = Build::Config->new("t/input/config_cache_clm4_5_test.xml");
  isnt($self->{cfg}, undef, (caller(0))[3] . " : config object created.");

  $self->{definition} = Build::NamelistDefinition->new("t/input/namelist_definition_clm4_5_test.xml");
  isnt($self->{definition}, undef, (caller(0))[3] . " : namelist_definition object created.");

  $self->{defaults} = Build::NamelistDefaults->new("t/input/namelist_defaults_clm4_5_test.xml");
  isnt($self->{defaults}, undef,  (caller(0))[3] . " : namelist_defaults object created.");

  $self->{physv} = config_files::clm_phys_vers->new( $self->{cfg}->get('phys') );
  isnt($self->{physv}, undef,  (caller(0))[3] . " : phys_vers object created.");

}

sub shutdown : Test(shutdown) {
  # cleanup the single instance test fixtures
}

sub setup : Test(setup => 1) {
  my $self = shift;
  # provide common fixture for all tests, create fresh for each test

  $self->{nl} = Build::Namelist->new();
  isnt($self->{nl}, undef, (caller(0))[3] . " : empty namelist object created.");
}

sub teardown : Test(teardown) {
  # clean up after test
}

#-------------------------------------------------------------------------------
#
# tests
# 
# 1: test that use_century_decom is the default with ed_mode
# 7: test that clm4.0 fails with ed_mode
#
#-------------------------------------------------------------------------------

# Test 1: use_century_decomp is default with ed_mode
sub test_ed_mode__use_century_decomp : Tests {
  my $self = shift;

  my $msg = "Test that the ed_mode set with use_century_decomp results success.\n";

  use CLMBuildNamelist qw(setup_cmdl_ed_mode);

  my $opts = { ed_mode => 1, };
  my $nl_flags = { crop => 1, };  # eq "on"

#  my $nl_flags = { phys => "clm4_5", };
#  my $physv = config_files::clm_phys_vers->new( "clm4_5" );

  CLMBuildNamelist::setup_cmdl_ed_mode($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}, $self->{physv});
  my $group = $self->{definition}->get_group_name("use_century_decomp");
  my $result = $self->{nl}->get_variable_value($group, "use_century_decomp");
  is($result, '.true.' ) || diag($msg);
  #isnt
}


1;
