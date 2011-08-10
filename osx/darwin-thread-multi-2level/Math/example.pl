#!/usr/bin/perl -w
use Math::Random qw(:all);
use Carp;
use strict;

my $max_choice = 17;
my $ans = "";
my $input = "";
my @args = ();
my @result = ();

 TEST: while (1) {
     print " Enter number corresponding to choice:\n",
     "      (0) Exit this program\n",
     "      (1) Generate Chi-Square deviates\n",
     "      (2) Generate noncentral Chi-Square deviates\n",
     "      (3) Generate F deviates\n",
     "      (4) Generate noncentral F  deviates\n",
     "      (5) Generate random permutation\n",
     "      (6) Generate uniform integers\n",
     "      (7) Generate uniform reals\n",
     "      (8) Generate beta deviates\n",
     "      (9) Generate binomial outcomes\n",
     "     (10) Generate Poisson outcomes\n",
     "     (11) Generate exponential deviates\n",
     "     (12) Generate gamma deviates\n",
     "     (13) Generate multinomial outcomes\n",
     "     (14) Generate normal deviates\n",
     "     (15) Generate negative binomial outcomes\n",
     "     (16) Generate multivariate normal deviates\n",
     "     (17) Generate random permuted index\n";
     $ans = <stdin>;
     chomp $ans;
     $ans = int($ans);
     last TEST if $ans == 0;
     unless ($ans > 0 and $ans <= $max_choice) {
	 print "Try one of (1, ... $max_choice)\n";
	 next TEST;
     }
     print "Enter phrase to initialize seeds:\n";
     my $phrase = <stdin>;
     chomp $phrase;
     random_set_seed_from_phrase($phrase);

     if ($ans == 1) {
	 print "Enter (space-separated) N, DF for chi-square deviates:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_chi_square(@args);
	 shift(@args);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('chis',@args));
	 next TEST;
     }

     if ($ans == 2) {
	 print "Enter (space-separated) N, DF, NONC for ",
	 "noncentral chi-square deviates:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_noncentral_chi_square(@args);
	 shift(@args);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('ncch',@args));
	 next TEST;
     }

     if ($ans == 3) {
	 print "Enter (space-separated) N, DFN, DFD for F deviates:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_f(@args);
	 shift(@args);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('f',@args));
	 next TEST;
     }

     if ($ans == 4) {
	 print "Enter (space-separated) N, DFN, DFD, NONC for ",
	 "non-central F deviates:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_noncentral_f(@args);
	 shift(@args);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('ncf',@args));
	 next TEST;
     }

     if ($ans == 5) {
	 print "Enter (space-separated) list for random permutation:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 print "Result (' : ' separated):\n",
	 join(" : ", random_permutation(@args)),"\n";
	 next TEST;
     }
	 
     if ($ans == 6) {
	 print "Enter (space-separated) Maximum integer, Replications ",
	 "per integer:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_uniform_integer($args[0] * $args[1], 1, $args[0]);
	 my @totals = (0) x ($args[0] + 1);
	 my $val = 0;
	 foreach $val (@result) { $totals[$val]++; }
	 shift @totals;
	 print "Result (' : ' separated):\n",join(" : ", @totals),"\n";
	 next TEST;
     }
	 
     if ($ans == 7) {
	 print "Enter (space-separated) N, LOWER, UPPER for ",
	 "uniform real deviates:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_uniform(@args);
	 shift(@args);
	 $args[0] = 0 unless defined($args[0]);
	 $args[1] = 1 unless defined($args[1]);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('unif',@args));
	 next TEST;
     }

     if ($ans == 8) {
	 print "Enter (space-separated) N, A, B for ",
	 "beta deviates:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_beta(@args);
	 shift(@args);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('beta',@args));
	 next TEST;
     }

     if ($ans == 9) {
	 print "Enter (space-separated) N, NTrials, P for ",
	 "binomial outcomes:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_binomial(@args);
	 shift(@args);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('bin',@args));
	 next TEST;
     }

     if ($ans == 10) {
	 print "Enter (space-separated) N, MU for ",
	 "poisson outcomes:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_poisson(@args);
	 shift(@args);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('pois',@args));
	 next TEST;
     }

     if ($ans == 11) {
	 print "Enter (space-separated) N, AV for ",
	 "exponential deviates:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_exponential(@args);
	 shift(@args);
	 $args[0] = 1 unless defined($args[0]);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('expo',@args));
	 next TEST;
     }

     if ($ans == 12) {
	 print "Enter (space-separated) N, A, R for ",
	 "gamma deviates:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_gamma(@args);
	 shift(@args);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('gamm',@args));
	 next TEST;
     }
     
     if ($ans == 13) {
	 print "Enter (space-separated) list of prob.s for categories ",
	 "for multinomial outcomes:\n";
	 my $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 print "Please enter number of events to be classified:\n";
	 my $n = <stdin>;
	 chomp $n;
	 @result = random_multinomial($n,@args);
	 print "Result:\n",join(" ", @result),"\n";
	 my $sum = 0;
	 my $val;
	 foreach $val (@result) {$sum += $val; }
	 foreach $val (@result) {$val /= $sum; }
	 print "Observed proportions:\n",join(" ", @result),"\n";
	 pop @args;
	 $sum = 0;
	 foreach $val (@args) {$sum += $val; }
	 push @args, (1 - $sum);
	 print "Expected proportions:\n",join(" ", @args),"\n";
	 next TEST;
     }
	 
     if ($ans == 14) {
	 print "Enter (space-separated) N, AV, SD for ",
	 "normal deviates:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_normal(@args);
	 shift(@args);
	 $args[0] = 0 unless defined($args[0]);
	 $args[1] = 1 unless defined($args[1]);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('norm',@args));
	 next TEST;
     }
	 
     if ($ans == 15) {
	 print "Enter (space-separated) N, NEvents, P for ",
	 "negative binomial outcomes:\n";
	 $input = <stdin>;
	 chomp $input;
	 @args = split(/\s+/,$input);
	 @result = random_negative_binomial(@args);
	 shift(@args);
	 print_results(sampstat(@result),scalar(@result),
		       true_stats('nbin',@args));
	 next TEST;
     }

     if ($ans == 16) {
	 print "Enter dimension of multivariate deviate:\n";
	 my $p = <stdin>;
	 chomp $p;
	 print "Enter mean vector of length $p (space separated):\n";
	 my $temp = <stdin>;
	 chomp $temp;
	 $temp =~ s/^\s*//;
	 my @mean;
	 @mean = split(/\s+/,$temp);
	 print "Enter (symmetric, $p by $p) covariance matrix\n",
	 "One space-separated row per line:\n";
	 my $val;
	 my @covariance = (0) x $p;
	 foreach $val (@covariance) {
	     $temp = <stdin>;
	     chomp $temp;
	     $temp =~ s/^\s*//;
	     $val = [ split(/\s+/,$temp) ];
	 }
	 print "Enter number of observations:\n";
	 my $n = <stdin>;
	 chomp $n;
	 my @ans = random_multivariate_normal($n,@mean,@covariance);
	 my $template = join(" ",('%15.7g') x $p);
	 print "\nResults:\n";
	 foreach $val (@ans) {
	     printf "$template\n",@{$val};
	 }
     }

     if ($ans == 17) {
	 print "Enter N (size of array) for random permuted index:\n";
	 $input = <stdin>;
	 chomp $input;
	 print "Result (' : ' separated):\n",
	 join(" : ", random_permuted_index($input)),"\n";
	 next TEST;
     }
 }
print "Normal termination of tester.\n";

sub sampstat { # gets sample statistics for array - returns array of stats
    my $n = scalar(@_);
    return () unless $n > 0;
    return ($_[0], 0, $_[0], $_[0]) unless $n > 1;
    my($min) = my($max) = $_[0];
    my $val;
    my $avg = 0;
    foreach $val (@_) {
	$avg += $val;
	$min = $val if $val < $min;
	$max = $val if $val > $max;
    }
    $avg /= $n;
    my $var = 0;
    foreach $val (@_) {
	$var += ($val - $avg)**2;
    }
    $var /= ($n - 1);
    return ($avg, $var, $min, $max);
}

sub print_results { # Arguments: $avg, $var, $min, $max, $nobs, $travg, $trvar
    my($avg, $var, $min, $max, $nobs, $travg, $trvar) = @_;
    print "Results:\n";
    printf "Number of observations: %d\n", $nobs;
    printf "Mean    : %15.7g   True    : %15.7g\n", $avg, $travg;
    printf "Variance: %15.7g   True    : %15.7g\n", $var, $trvar;
    printf "Minimum : %15.7g   Maximum : %15.7g\n", $min, $max;
}

sub true_stats { # Arguments: $type, @parin; Returns: $av, $var

########################################################################
#     Returns mean and variance for a number of statistical distribution
#     as a function of their parameters.
#
#
#                              Arguments
#
#
#    $type --> Character string indicating type of distribution
#             'chis' chisquare
#             'ncch' noncentral chisquare
#             'f'    F (variance ratio)
#             'ncf'  noncentral f
#             'unif' uniform
#             'beta' beta distribution
#             'bin'  binomial
#             'pois' poisson
#             'expo' exponential
#             'gamm' gamma
#             'norm' normal
#             'nbin' negative binomial
#
#    @parin --> Array containing parameters of distribution
#              chisquare
#               $parin[0] is df
#              noncentral chisquare
#               $parin[0] is df
#               $parin[1] is noncentrality parameter
#              F (variance ratio)
#               $parin[0] is df numerator
#               $parin[1] is df denominator
#              noncentral F
#               $parin[0] is df numerator
#               $parin[1] is df denominator
#               $parin[2] is noncentrality parameter
#              uniform
#               $parin[0] is LOW bound
#               $parin[1] is HIGH bound
#              beta
#               $parin[0] is A
#               $parin[1] is B
#              binomial
#               $parin[0] is Number of trials
#               $parin[1] is Prob Event at Each Trial
#              poisson
#               $parin[0] is Mean
#              exponential
#               $parin[0] is Mean
#              gamma
#               $parin[0] is A
#               $parin[1] is R
#              normal
#               $parin[0] is Mean
#               $parin[1] is Standard Deviation
#              negative binomial
#               $parin[0] is required Number of events
#               $parin[1] is Probability of event
#
#     $av <-- Mean of specified distribution with specified parameters
#
#     $var <-- Variance of specified distribution with specified paramete
#
#
#                              Note
#
#
#     $av and $var will be returned -1 if mean or variance is infinite
#
#**********************************************************************
    my $type = shift(@_);
    my @parin = @_;
    my($av, $var, $a, $b, $range) = (-1,-1,0,0,0);

  TYPE: {
      if (('chis') eq ($type)){
	  $av = $parin[0];
	  $var = 2.0*$parin[0];
	  last TYPE;}
      
      if (('ncch') eq ($type)) {
	  $a = $parin[0] + $parin[1];
	  $b = $parin[1]/$a;
	  $av = $a;
	  $var = 2.0*$a* (1.0+$b);
	  last TYPE;}

      if (('f') eq ($type)) {
	  unless ($parin[1] <= 2.0001) {
	      $av = $parin[1]/ ($parin[1]-2.0);
	  }
	  unless ($parin[1] <= 4.0001) {
	      $var = (2.0*$parin[1]**2* ($parin[0]+$parin[1]-2.0))/
		  ($parin[0]* ($parin[1]-2.0)**2* ($parin[1]-4.0));
	  }
	  last TYPE;}
      
      if (('ncf') eq ($type)) {
	  unless ($parin[1] <= 2.0001){
	      $av = ($parin[1]* ($parin[0]+$parin[2]))/
		  (($parin[1]-2.0)*$parin[0]);
	  }
	  unless ($parin[1] <= 4.0001) {
	      $a = ($parin[0]+$parin[2])**2 + ($parin[0]+2.0*$parin[2])*
		  ($parin[1]-2.0);
	      $b = ($parin[1]-2.0)**2* ($parin[1]-4.0);
	      $var = 2.0* ($parin[1]/$parin[0])**2* ($a/$b);
	  }
	  last TYPE;}

      if (('unif') eq ($type)) {
	  $range = $parin[1] - $parin[0];
	  $av = $parin[0] + $range/2.0;
	  $var = $range**2/12.0;
	  last TYPE;}

      if (('beta') eq ($type)) { 
	  $av = $parin[0]/ ($parin[0]+$parin[1]);
	  $var = ($av*$parin[1])/ (($parin[0]+$parin[1])*
				   ($parin[0]+$parin[1]+1.0));
	  last TYPE;}

      if (('bin') eq ($type)) { 
	  $av = $parin[0]*$parin[1];
	  $var = $av* (1.0-$parin[1]);
	  last TYPE;}

      if (('pois') eq ($type)) { 
	  $av = $parin[0];
	  $var = $parin[0];
	  last TYPE;}

      if (('expo') eq ($type)) { 
	  $av = $parin[0];
	  $var = $parin[0]**2;
	  last TYPE;}

      if (('gamm') eq ($type)) { 
	  $av = $parin[1] / $parin[0];
	  $var = $av / $parin[0];
	  last TYPE;}

      if (('norm') eq ($type)) { 
	  $av = $parin[0];
	  $var = $parin[1]**2;
	  last TYPE;}

      if (('nbin') eq ($type)) { 
	  $av = $parin[0] * (1.0 - $parin[1]) / $parin[1];
	  $var = $av / $parin[1];
	  last TYPE;}

      croak "Unimplemented \$type: $type in true_stats";
  }
    return ($av,$var);
}
