#!/usr/bin/perl
use strict;
use warnings;

my $now_string;
$now_string = localtime;
my $comment;
$comment=$ARGV[0];



my $filename = 'MEND_namelist.nml';
my $folder;
open(FH, '<', $filename) or die $!;
while(<FH>){
   my $line=$_;
   chomp($line);
   if ($line=~ /Dir_Output  = /){
      $folder=$';
	  $folder=~s/'//g;
	  $folder=~s/\//\\/g;
   

}
}
close(FH);

system("mkdir $folder ");
my $logfile = 'log.txt';
open(L, '<', $logfile) or die $!;

open(FH, '>', $folder.$logfile) or die $!;

while(<L>){
   my $line=$_;
   chomp($line);
   print FH "$line\n";
     

}
print FH "$folder\n";
print FH "$now_string\n";

close(FH);
close(L);
system("cp $folder\\$logfile  log.txt");


system("echo $comment > $folder\\comment.txt");
system("type MEND_namelist.nml > $folder$filename");
system(".\\dist\\Debug\\Cygwin-Windows\\mendokw.exe");
