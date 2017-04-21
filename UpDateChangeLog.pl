#!/usr/bin/env perl
#=======================================================================
#
#  This is a script to update the ChangeLog
#
# Usage:
#
# perl ChangeLog tag-name One-line summary
#
#
#=======================================================================

use strict;
use Getopt::Long;
use IO::File;
#use warnings;
#use diagnostics;

use English;

my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!; # name of program
my $ProgDir = $1;                         # name of directory where program lives

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options] <tag-name> <one-line-summary>

OPTIONS
     -compbrnch version     Enter clm branch  version to compare to (under branch_tags in repo).
      [or -cb]
     -comptrunk version     Enter clm trunk version to compare to (under trunk_tags in repo).
      [or -ct]
     -help   [or -h]        Help on this script.
     -update [or -u]        Just update the date/time for the latest tag
                            In this case no other arguments should be given.
ARGUMENTS
     <tag-name>             Tag name of tag to document
     <one-line-summary>     Short summary description of this tag
EXAMPLES:
     To just update the date/time for the latest tag

     $ProgName -update

     To document a new tag

     $ProgName clm4_5_2_r097 "Description of this tag"

     To document a new tag and compare expected fails to previous tag.

     $ProgName clm4_5_2_r097 "Description of this tag" -ct clm4_5_2_r096
EOF
}

my %opts = { 
               help      => 0,
               update    => 0,
               comptrunk => undef,
               compbrnch => undef,
           };
GetOptions( 
    "h|help"         => \$opts{'help'},
    "u|update"       => \$opts{'update'},
    "ct|comptrunk=s" => \$opts{'comptrunk'},
    "cb|compbrnch=s" => \$opts{'compbrnch'},
   );
if ( $opts{'help'} ) {
  usage();
}
my $tag; my $sum;

if ( ! $opts{'update'} ) {
   if ( $#ARGV != 1 ) {
     print "ERROR: wrong number of arguments: $ARGV\n";
     usage();
   }

   $tag = $ARGV[0];
   $sum = $ARGV[1];

   if ( $tag !~ /clm[0-9]+_([0-9]+)_[0-9]+_r[0-9]+/ ) {
     print "ERROR: bad tagname: $tag\n";
     usage();
   }
} else {
   if ( $#ARGV != -1 ) {
     print "ERROR: wrong number of arguments when update option picked: $ARGV\n";
     usage();
   }
}
my $EDITOR = $ENV{EDITOR};
if ( $EDITOR !~ /.+/ ) {
  print "ERROR: editor NOT set -- set the env variable EDITOR to the text editor you would like to use\n";
  usage();
}


my $template      = ".ChangeLog_template";
my $changelog     = "ChangeLog";
my $changesum     = "ChangeSum";
my $changelog_tmp = "ChangeLog.tmp";
my $changesum_tmp = "ChangeSum.tmp";

my $user = $ENV{USER};
if ( $user !~ /.+/ )  {
  die "ERROR: Could not get user name: $user";
}
my @list = getpwnam( $user );
my $fullname = $list[6];
my $date = `date`;
chomp( $date );

if ( $date !~ /.+/ )  {
  die "ERROR: Could not get date: $date\n";
}

#
# Deal with ChangeLog file
#
my $fh = IO::File->new($changelog_tmp, '>') or die "** $ProgName - can't open file: $changelog_tmp\n";

#
# If adding a new tag -- read in template and add information in
#
if ( ! $opts{'update'} ) {
   open( TL, "<$template"     )  || die "ERROR:: trouble opening file: $template";
   while( $_ = <TL> ) {
     if (      $_ =~ /Tag name:/ ) {
        chomp( $_ );
        print $fh "$_ $tag\n";
     } elsif ( $_ =~ /Originator/ ) {
        chomp( $_ );
        print $fh "$_ $user ($fullname)\n";
     } elsif ( $_ =~ /Date:/ ) {
        chomp( $_ );
        print $fh "$_ $date\n";
     } elsif ( $_ =~ /One-line Summary:/ ) {
        chomp( $_ );
        print $fh "$_ $sum\n";
     } elsif ( $_ =~ /CLM tag used for the baseline comparison tests if applicable:/ ) {
        chomp( $_ );
        if (      defined($opts{'comptrunk'}) ) {
           print $fh "$_ $opts{'comptrunk'}\n";
           &AddExpectedFailDiff( $fh, "trunk_tags/$opts{'comptrunk'}" );
        } elsif ( defined($opts{'compbrnch'}) ) {
           print $fh "$_ $opts{'compbrnch'}\n";
           &AddExpectedFailDiff( $fh, "branch_tags/$opts{'compbrnch'}" );
        } else {
           print $fh "$_\n";
        }
     } else {
        print $fh $_;
     }
   }
   close( TL );
}
open( CL, "<$changelog"     ) || die "ERROR:: trouble opening file: $changelog";
my $update = $opts{'update'};
my $oldTag = "";
while( $_ = <CL> ) {
  # If adding a new tag check that new tag name does NOT match any old tag
  if (  $_ =~ /Tag name:[   ]*(clm.+)/ ) {
     $oldTag = $1;
     if ( (! $opts{'update'}) && ($tag eq $oldTag) ) {
        close( CL );
        close( $fh );
        system( "/bin/rm -f $changelog_tmp" );
        print "ERROR:: New tag $tag matches a old tag name\n";
        usage();
     }
  # If updating the date -- find first occurance of data and change it 
  # Then turn the update option to off
  } elsif ( ($update) && ($_ =~ /(Date:)/) ) {
     print $fh "Date: $date\n";
     print "Update $oldTag with new date: $date\n";
     $update = undef;
     $_ = <CL>;
  }
  print $fh $_;
}
# Close files and move to final name
close( CL );
$fh->close( );
system( "/bin/mv    $changelog_tmp $changelog" );
#
# Deal with ChangeSum file
#

open( FH, ">$changesum_tmp" ) || die "ERROR:: trouble opening file: $changesum_tmp";

open( CS, "<$changesum"     ) || die "ERROR:: trouble opening file: $changesum";

my $update = $opts{'update'};

$date = `date "+%m/%d/%Y"`;
chomp( $date );

while( $_ = <CS> ) {
  # Find header line
  if ( $_ =~ /=====================/ ) {
     print FH $_;
     my $format = "%16.16s %8.8s %10.10s %s\n";
     if ( $update ) {
       $_ = <CS>;
       if ( /^(.{16}) (.{8}) (.{10}) (.+)$/ ) {
          $tag  = $1;
          $user = $2;
          $sum  = $4;
       } else {
          die "ERROR: bad format for ChangeSum file\n";
       }
     }
     printf FH $format, $tag, $user, $date, $sum;
     $_ = <CS>;
  }
  print FH $_;
}
# Close files and move to final name
close( CS );
close( FH );
system( "/bin/mv    $changesum_tmp $changesum" );

#
# Edit the files
#
if ( ! $opts{'update'} ) {
  system( "$EDITOR $changelog" );
  system( "$EDITOR $changesum" );
}
system( "/bin/cp -fp $changelog components/clm/doc/." );
system( "/bin/cp -fp $changesum components/clm/doc/." );
system( "/bin/chmod 0444 components/clm/doc/$changelog" );
system( "/bin/chmod 0444 components/clm/doc/$changesum" );

sub AddExpectedFailDiff {
#
# Add information about the expected fail difference
#
  my $fh      = shift;
  my $version = shift;

  my $SVN_MOD_URL  = "https://svn-ccsm-models.cgd.ucar.edu/clm2/";
  my $expectedFail = `find . -name 'expected*Fail*.xml' -print`;
  if ( $expectedFail eq "" ) {
     die "ERROR:: expectedFails file NOT found here\n";
  }

  `svn ls $SVN_MOD_URL/$version`               || die "ERROR:: Bad version to compare to: $version\n";
  `svn ls $SVN_MOD_URL/$version/$expectedFail` || die "ERROR:: expectedFails file NOT found in: $version\n";
  print $fh "\nDifference in expected fails from testing:\n\n";
  my $diff = `svn diff --old $SVN_MOD_URL/$version/$expectedFail \ \n --new $expectedFail`;
  if ( $diff eq "" ) {
     print $fh "    No change in expected failures in testing\n";
  } else {
     print $fh $diff;
  }
}
