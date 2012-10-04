#! /usr/bin/perl
# Author: Eric Lowe

use warnings;
use Linux::Inotify2;

my $file = '/home/elowe/phylosift_uploads/';
my $follower = '/home/koadman/bin/follower';
my $dir = q();
my $ps = '/home/elowe/PhyloSift/bin/phylosift';
# calls inotify constructor
my $notifier = new Linux::Inotify2 or die "Unable to create new inotify object: $!"; 
 # tells inotify where to watch and what to watch for
$notifier->watch('/home/elowe/phylosift_uploads', IN_CREATE, sub {
	my $e = shift;
	$file .= $e->name;
	$dir = $e->{name};
	print "File is named $file\n";
});
while () 
{
    my @events = $notifier->read; # array of events from watched directory
    unless (@events > 0) 
    {
	print "read error: $!";
	last;
    }
    my $out_path = "/home/elowe/public_html/$dir";
    my $ps_command = "$follower $file | $ps all --chunks=1 --stdin --output=$out_path blurg";

    #printf "mask\t%d\n", $_->mask foreach @events;
    
    system("$ps_command");
}

$notifier->close(); # inotify destructor
exit;
