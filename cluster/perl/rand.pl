#!/usr/bin/perl
# selecting N random lines from a text file

die "Usage: $0 <N> file, where N is the number of lines to pick\n"
if @ARGV < 1;
$N = shift @ARGV;

@pick = ();
while (<>) {
if (@pick < $N) {
push @pick, $_;
($r1, $r2) = (rand(@pick), rand(@pick));
($pick[$r1], $pick[$r2]) = ($pick[$r2], $pick[$r1]);
} else {
rand($.) <= $N and $pick[rand(@pick)] = $_;
}
}

print @pick;
