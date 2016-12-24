#!/usr/bin/env perl

if ($ARGV[1] eq "orf") {
  $hmmfile  = "ptr.std.e100.mer.orf";
  $consfile = "tbn.w3.e100.cons.mer.orf";
  $fpwfile = "tbn.w3.e100.fpw.mer.orf";
} else {
  $hmmfile  = "ptr.std.e100.mer";
  $consfile = "tbn.w3.e100.cons.mer";
  $fpwfile = "tbn.w3.e100.fpw.mer";
}

$eval_min_threshold = 1e-170;

open F, "<$hmmfile";
while ($l = <F>) {
   next if $l =~ /^#/;
   next if $l !~ /^=/;
#   printf("hmm line %10s\n", $l);
   @vals = split /\s+/, $l;
   ($fam, $famseq, $eval, $plusminus) = @vals[2,3,4,5];
   next if ($fam eq "family");
   next if ($plusminus eq "-");
   # skip evalues > 100 until phmmert -domE filter fixed
   next  if ($eval > 100);
   # skip evalues < a threshold 
#   printf("phmmert file %10s %10s\n",$fam,$eval);
   
   $hmm_hits{$fam.$famseq} = 0.0+$eval;   
}

open F, "<$consfile";
while ($l = <F>) {
   next if $l =~ /^#/;
   next if $l !~ /^=/;
#   printf("cons line %10s\n", $l);
   @vals = split /\s+/, $l;
   ($fam, $famseq, $eval, $plusminus) = @vals[2,3,4,5];
   next if ($fam eq "family");
   next if ($plusminus eq "-");
   # skip evalues < a threshold 
#   printf("cons file %10s %10s\n",$fam,$eval);
   next if ($fam eq "family");
#  
   $cons_hits{$fam.$famseq} = 0.0+$eval;   
}

open F, "<$fpwfile";
while ($l = <F>) {
   next if $l =~ /^#/;
   next if $l !~ /^=/;
#   printf("fpw line %10s\n", $l);
   @vals = split /\s+/, $l;
   ($fam, $famseq, $eval, $plusminus) = @vals[2,3,4,5];
   next if ($fam eq "family");
   next if ($plusminus eq "-");
   # skip evalues < a threshold 
#   printf("fpw file %10s %10s\n",$fam,$eval);
   next if ($fam eq "family");
#
#   printf("fpw family:%10s sequence:%10s E-value:%10s\n",$fam, $famseq,0.0+$eval);  
   $fpw_hits{$fam.$famseq} = 0.0+$eval;   
}


if ($ARGV[0] eq "fpw") {
   @gis = sort {$fpw_hits{$a} <=> $fpw_hits{$b}} keys %fpw_hits
} elsif ($ARGV[0] eq "cons") {
   @gis = sort {$cons_hits{$a} <=> $cons_hits{$b}} keys %cons_hits
} elsif ($ARGV[0] eq "hmm") {
#   printf("organizing by hmmi\n");
   @gis = sort {$hmm_hits{$a} <=> $hmm_hits{$b}} keys %hmm_hits
}

printf ( "%50s\t%12s\t%12s\t%12s\n", "gi",  "tblastnfpw", "tblastncons", "phmmert");
   
foreach $k (@gis) {
   #if there there was no hit for a particular
   #tool then give it an E-Value of 100 so 
   #a placeholder for it; this helps when we plot the hists
   if (exists $fpw_hits{$k}){
      $f = $fpw_hits{$k} 
   }
   else {
      $f = 100;
   }

   if (exists $cons_hits{$k}) {
      $c = $cons_hits{$k}
   }
   else {
      $c = 100;
   }

   if (exists $hmm_hits{$k}) {
      $h = $hmm_hits{$k}
   }
   else {
      $h = 100;
   }

   #skip hit if all search tools found the  hit and any search tool recorded an E-Value less
   #than the threshold; we will not be interested in those
   next if (($f < $eval_min_threshold) or
           ($c < $eval_min_threshold) or
           ($h < $eval_min_threshold));


   printf ( "%50s\t%12g\t%12g\t%12g\n", $k, $f, $c, $h);
}



#foreach $k (sort {$b <=> $a} keys %cons_hits) {
#   print "$k [$cons_hits{$k}] [$hmm_hits{$k}]\n";
#}

