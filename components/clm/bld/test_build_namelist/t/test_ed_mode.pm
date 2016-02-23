package test_ed_mode;

# Unit tests for function: ed_mode

use Data::Dumper;
use Test::More;
use Test::Exception;

use parent qw(Test::Class);

use CLMBuildNamelist qw(setup_logic_lnd_frac);

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

  $self->{test_files} = 0;

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
# 2: test that use_vertsoilc is the default with ed_mode
# 3: test that use_ed_spit_fire is the default with ed_mode
# 4: test that use_ed = .false. will fail for ed_mode
# 5: test that crop = "on" will fail for ed_mode
# 6: test that clm4_0 will fail for ed_mode (NOT AVAILABLE)
# 7: test that use_ed_spit_fire = false generates a false
# 8: test that use_versoilc = false generates a false
# 9: test that use_century_decomp = false generates a false
#
#-------------------------------------------------------------------------------

# Test 1: use_century_decomp is default with ed_mode
sub test_ed_mode__use_century_decomp : Tests {
  my $self = shift;

  my $msg = "Test that use_century_decomp is the default on ed_mode.\n";

  use CLMBuildNamelist qw(setup_cmdl_ed_mode);
  use CLMBuildNamelist qw(setup_logic_ed);

  my $opts = { ed_mode => 1, bgc => "default" };
  my $nl_flags = { crop => "off", }; # the setup_cmdl_ed_mode code requires this logic

  # include bgc because they have mutually exclusive functionality that needs to be tested
  CLMBuildNamelist::setup_cmdl_bgc($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}, $self->{cfg},$self->{physv});
  CLMBuildNamelist::setup_cmdl_ed_mode($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}, $self->{physv});
  $nl_flags->{"use_ed"} = $self->{nl}->get_value("use_ed");

  CLMBuildNamelist::setup_logic_ed($self->{test_files}, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}, $self->{physv});
  my $result = $self->{nl}->get_variable_value($group, "use_century_decomp");
  is($result, '.true.' ) || diag($msg);
  
}

# Test 2: use_vertsoilc is default with ed_mode
sub test_ed_mode__use_vertsoilc : Tests {
    my $self = shift;

    my $msg = "Tests that use_vertsoilc is default on ed_mode.\n";

    use CLMBuildNamelist qw(setup_cmdl_ed_mode);

    my $opts = { ed_mode => 1, };
    my $nl_flags = { crop => "off", };

    CLMBuildNamelist::setup_cmdl_ed_mode($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}, $self->{physv});
    my $group = $self->{definition}->get_group_name("use_vertsoilc");
    my $result = $self->{nl}->get_variable_value($group, "use_vertsoilc");
    is($result, '.true.' ) || diag($msg);
}

# Test 3: use_ed_spit_fire is default with ed_mode
sub test_ed_mode__use_ed_spit_fire : Tests {
    my $self = shift;

    my $msg = "Tests that use_ed_spit_fire is default on ed_mode.\n";

    use CLMBuildNamelist qw(setup_cmdl_ed_mode);

    my $opts = { ed_mode => 1, };
    my $nl_flags = { crop => "off", };

    CLMBuildNamelist::setup_cmdl_ed_mode($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}, $self->{physv});
    my $group = $self->{definition}->get_group_name("use_ed_spit_fire");
    my $result = $self->{nl}->get_variable_value($group, "use_ed_spit_fire");
    is($result, '.true.' ) || diag($msg);
}

# Test 4: use_ed=false will crash on ed_mode
sub test_ed_mode__use_ed_nl_contradicts_cmdl : Tests {
    my $self = shift;

    my $msg = "Tests that use_ed=false will fail on ed_mode.\n";

    use CLMBuildNamelist qw(setup_cmdl_ed_mode);

    my $opts = { ed_mode => 1, };
    my $nl_flags = { crop => "off", use_ed => ".false."};

    my $group = $self->{definition}->get_group_name("use_ed");
    $self->{nl}->set_variable_value($group, "use_ed", '.false.' );

    dies_ok(sub { CLMBuildNamelist::setup_cmdl_ed_mode($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl},
							$self->{physv}) }) || diag($msg);
}

# Test 5: crop=on will crash on ed_mode
sub test_ed_mode__crop_on_nl_contradicts_cmdl : Tests {
    my $self = shift;

    my $msg = "Tests that crop=on will fail on ed_mode.\n";

    use CLMBuildNamelist qw(setup_cmdl_ed_mode);

    my $opts = { ed_mode => 1, };
    my $nl_flags = { crop => "on", };

    my $group = $self->{definition}->get_group_name("crop");
    $self->{nl}->set_variable_value($group, "crop", 'on' );

    dies_ok(sub { CLMBuildNamelist::setup_cmdl_ed_mode($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl},
						   $self->{physv}) }) || diag($msg);
}

# 7: test that use_ed_spit_fire = false generates a false
sub test_ed_mode__use_ed_spit_fire_false : Tests {
    my $self = shift;

    my $msg = "Tests that use_ed_spit_fire can be turned to false.\n";

    use CLMBuildNamelist qw(setup_cmdl_ed_mode);

    my $opts = { ed_mode => 1, };
    my $nl_flags = { crop => "off", };

    my $group = $self->{definition}->get_group_name("use_ed_spit_fire");
    $self->{nl}->set_variable_value($group, "use_ed_spit_fire", '.false.' );

    CLMBuildNamelist::setup_cmdl_ed_mode($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl},
						       $self->{physv});

    my $result = $self->{nl}->get_variable_value($group, "use_ed_spit_fire");
    is($result, '.false.' ) || diag($msg);

}

# 8: test that use_versoilc = false generates a false
sub test_ed_mode__use_vertsoilc_false : Tests {
    my $self = shift;

    my $msg = "Tests that use_vertsoilc can be turned to false.\n";

    use CLMBuildNamelist qw(setup_cmdl_ed_mode);

    my $opts = { ed_mode => 1, };
    my $nl_flags = { crop => "off", };

    my $group = $self->{definition}->get_group_name("use_vertsoilc");
    $self->{nl}->set_variable_value($group, "use_vertsoilc", '.false.' );

    CLMBuildNamelist::setup_cmdl_ed_mode($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl},
					 $self->{physv});

    my $result = $self->{nl}->get_variable_value($group, "use_vertsoilc");
    is($result, '.false.' ) || diag($msg);

}

# 9: test that use_ed_spit_fire = false generates a false
sub test_ed_mode__use_century_decomp_false : Tests {
    my $self = shift;

    my $msg = "Tests that use_century_decomp can be turned to false.\n";

    use CLMBuildNamelist qw(setup_cmdl_ed_mode);

    my $opts = { ed_mode => 1, };
    my $nl_flags = { crop => "off", };

    my $group = $self->{definition}->get_group_name("use_century_decomp");
    $self->{nl}->set_variable_value($group, "use_century_decomp", '.false.' );

    CLMBuildNamelist::setup_cmdl_ed_mode($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl},
                                         $self->{physv});

    my $result = $self->{nl}->get_variable_value($group, "use_century_decomp");
    is($result, '.false.' ) || diag($msg);

}



1;
