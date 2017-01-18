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

#set this to the minimum E-Value to be included in the table
#any tool with an evalue equal to or less than this will
#not be included in the table
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


#if ($ARGV[0] eq "fpw") {
#   @fpwgis = sort {$fpw_hits{$a} <=> $fpw_hits{$b}} keys %fpw_hits
#} elsif ($ARGV[0] eq "cons") {
#   @consgis = sort {$cons_hits{$a} <=> $cons_hits{$b}} keys %cons_hits
#} elsif ($ARGV[0] eq "hmm") {
#   printf("organizing by hmmi\n");
#   @hmmgis = sort {$hmm_hits{$a} <=> $hmm_hits{$b}} keys %hmm_hits
#}

#merge the hits found by each tool so that there is one
#list with all the hits even though one tool may have found
#the positive sequence (hit) and no other tools found it
@hmmgis = keys %hmm_hits;
@fpwgis = keys %fpw_hits;
@consgis = keys %cons_hits;

@hmmseen{@hmmgis} = ();
@hmm_cons_hits = (@hmmgis, grep{!exists $hmmseen{$_}} @consgis);
@hmmconsseen{@hmm_cons_hits} = ();
@all_hits = (@hmm_cons_hits, grep{!exists $hmmconsseen{$_}} @fpwgis);


foreach $k (@all_hits) {
   #if there there was no hit for a particular
   #tool then give it an E-Value of -1 as 
   #a placeholder for it; this helps when we plot the histograms
   if (exists $fpw_hits{$k})  {
       ;
   }
   else {
      $fpw_hits{$k} = -1;
   }

   if (exists $cons_hits{$k}) {
       ;
   }
   else {
      $cons_hits{$k} = -1;
   }

   if (exists $hmm_hits{$k}) {
      ;
   }
   else {
      $hmm_hits{$k} = -1;
   }
}

# sort the E-Values by tool based on the input argument
if ($ARGV[0] eq "fpw") {
   @sorted_keys = sort {$fpw_hits{$a} <=> $fpw_hits{$b}} keys %fpw_hits
} elsif ($ARGV[0] eq "cons") {
   @sorted_keys = sort {$cons_hits{$a} <=> $cons_hits{$b}} keys %cons_hits
} elsif ($ARGV[0] eq "hmm") {
#   printf("organizing by hmmi\n");
   @sorted_keys = sort {$hmm_hits{$a} <=> $hmm_hits{$b}} keys %hmm_hits
}

printf ( "%50s\t%12s\t%12s\t%12s\n", "gi",  "tblastnfpw", "tblastncons", "phmmert");

foreach $k (@sorted_keys) {
   #skip hit if all search tools found the  hit and any search tool recorded an E-Value less
   #than the threshold; we will not be interested in those
   next if ((($fpw_hits{$k} < $eval_min_threshold) and ($fpw_hits{$k} > -1)) or
           (($cons_hits{$k} < $eval_min_threshold) and ($cons_hits{$k} > -1)) or
           (($hmm_hits{$k} < $eval_min_threshold) and ($hmm_hits{$k} > -1)) );


   printf ( "%50s\t%12g\t%12g\t%12g\n", $k, $fpw_hits{$k}, $cons_hits{$k}, $hmm_hits{$k});
}



