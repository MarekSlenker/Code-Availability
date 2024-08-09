package polyconfig;

# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018-20

use strict;
require Exporter;

our @ISA = qw( Exporter );

our @EXPORT = qw( 
  get_label_from_rules
  @diploids @polyploids @subgenomes 
  %outgroups $ROOT %sister_clades @CODES $NODEORDER
  $MINBLOCKLENGTH $MAXGAPSPERBLOCK $MINBLOCKOVERLAP
);

# Abbreviated names of diploid species as found in FASTA and tree files.
# See %sister_clades below for how to define sister species.
our @diploids = qw( matthioli rivularis anatolica SharGramos Dinaric vardousPindicola lazica EBalkan ); 


# note these are diploids as well; hash instead of list
our %outgroups = ( 
  'lazica',1
);

# diploid used to root trees 
our $ROOT = 'lazica';

# Abbreviated names of polyploid species as found in FASTA and tree files.
our @polyploids = ('acraC018.101.h1_acrisPP','acraC018.101.h2_acrisPP','acraC018.101.h3_acrisPP','acraC018.101.h4_acrisPP',
                   'acraC095.109.h1_acrisPP','acraC095.109.h2_acrisPP','acraC095.109.h3_acrisPP',
                   'acraC192.3.h1_acrisPP','acraC192.3.h2_acrisPP','acraC192.3.h3_acrisPP','acraC192.3.h4_acrisPP');

# Abbreviated names of polyploid subgenomes, defined by user.
# These are used to annotate alleles from the subgenomes as ancestral.
# Leave empty if not required
our @subgenomes = ();


# Abbreviated names of labelled polyploid sequences as found in FASTA and tree files.
# Note that the capital letters correpond to @CODES below
our @polyploids_labelled = ( 'acraC018.101.h1_acrisPP','acraC018.101.h2_acrisPP','acraC018.101.h3_acrisPP','acraC018.101.h4_acrisPP','acraC095.109.h1_acrisPP','acraC095.109.h2_acrisPP','acraC095.109.h3_acrisPP','acraC192.3.h1_acrisPP','acraC192.3.h2_acrisPP','acraC192.3.h3_acrisPP','acraC192.3.h4_acrisPP' );


# Optional custom definition of clades that contain >1 diploids, if any.
# Leave empty or comment out the examples otherwise.
# MRCA nodes are internal node names, can be used in rules below; there are two keys:
# # i) MRCA for all species in all sister clades, also called a bifurcation in the code
# # ii) MRCA for each explicitely defined clade (there should be two)
# # finally the arrays contain lists of species in each clade, should be diploid
our %sister_clades = ();



# see rules defined below
our @CODES = qw( matthioli rivularis anatolica SharGramos Dinaric vardousPindicola lazica EBalkan all ); 

# use to ladderize trees
our $NODEORDER = 0; # 1:increasing, 0:decreasing

# Takes 4 parameters: 
# 1) string with abbreviated name of diploid/MRCA taxon (ancestor)
# 2) string with abbreviated name of diploid/MRCA taxon (descendant)
# 3) boolean scalar to indicate ancestor is sister
# 4) boolean scalar to indicate descendat is sister
# Returns a label, which is either a code from @CODES or '-' otherwise
sub get_label_from_rules {

  my ($anc_dip_taxon, $desc_dip_taxon, $anc_is_sister, $desc_is_sister) = @_;
	
  my  $lineage_code = '-';

  # check input params
  my $ancOK = 0;
  if(grep(/^$anc_dip_taxon/,@diploids)){ $ancOK = 1 }
  if(defined($sister_clades{$anc_dip_taxon})){ $ancOK = 1 }
  else {
    foreach my $MRCA (keys(%sister_clades)){
      if(defined($sister_clades{$MRCA}{$anc_dip_taxon})){
        $ancOK = 1; 
        last;
      }
    }
  }	
  if($ancOK == 0){
    print "# ERROR (get_label_from_rules): unrecognized ancestor $anc_dip_taxon\n";
    return $lineage_code;
  }	

  my $descOK = 0;
  if(grep(/^$desc_dip_taxon/,@diploids)){ $descOK = 1 }
  if(defined($sister_clades{$desc_dip_taxon})){ $descOK = 1 }
  else {
    foreach my $MRCA (keys(%sister_clades)){
      if(defined($sister_clades{$MRCA}{$desc_dip_taxon})){
        $descOK = 1;
        last;
      }
    }
  }
  if($descOK == 0){
    print "# ERROR (get_label_from_rules): unrecognized descendant $desc_dip_taxon\n";
    return $lineage_code;
  }

  # start applying rules
	
  ## ancestor is sister or descendant is empty, 
  ## only ancestor diploid is looked up
  if($anc_is_sister == 1 || $desc_dip_taxon eq ''){
    if($anc_dip_taxon eq 'Bsta'){ $lineage_code = 'B' }

    elsif($anc_dip_taxon eq 'matthioli'){ $lineage_code = 'matthioli' }
    elsif($anc_dip_taxon eq 'rivularis'){ $lineage_code = 'rivularis' }
    elsif($anc_dip_taxon eq 'anatolica'){ $lineage_code = 'anatolica' }
    elsif($anc_dip_taxon eq 'SharGramos'){ $lineage_code = 'SharGramos' }
    elsif($anc_dip_taxon eq 'Dinaric'){ $lineage_code = 'Dinaric' }
    elsif($anc_dip_taxon eq 'vardousPindicola'){ $lineage_code = 'vardousPindicola' }
    elsif($anc_dip_taxon eq 'lazica'){ $lineage_code = 'lazica' }
    elsif($anc_dip_taxon eq 'EBalkan'){ $lineage_code = 'EBalkan' }


  }
  else { ## both ancestor/descendant diploids/clades are considered

    if($anc_dip_taxon eq 'Hvul' && $desc_dip_taxon eq 'Bsta'){ 
      $lineage_code = 'A' 
    }
    elsif($anc_dip_taxon eq 'matthioli' && $desc_dip_taxon eq 'rivularis'){ 
      $lineage_code = 'rivularis' 
    }
    elsif($anc_dip_taxon eq 'rivularis' && $desc_dip_taxon eq 'matthioli'){ 
      $lineage_code = 'matthioli' 
    }

    elsif($anc_dip_taxon eq 'lazica' && $desc_dip_taxon eq 'rivularis'){ 
      $lineage_code = 'rivularis' 
    }
    elsif($anc_dip_taxon eq 'lazica' && $desc_dip_taxon eq 'matthioli'){ 
      $lineage_code = 'matthioli' 
    }

    elsif($anc_dip_taxon eq 'rivularis' && $desc_dip_taxon eq 'SharGramos'){ 
      $lineage_code = 'SharGramos' 
    }
    elsif($anc_dip_taxon eq 'matthioli' && $desc_dip_taxon eq 'SharGramos'){ 
      $lineage_code = 'SharGramos' 
    }

    elsif($anc_dip_taxon eq 'SharGramos' && $desc_dip_taxon eq 'anatolica'){ 
      $lineage_code = 'anatolica' 
    }
    elsif($anc_dip_taxon eq 'anatolica' && $desc_dip_taxon eq 'vardousPindicola'){ 
      $lineage_code = 'vardousPindicola' 
    }
    
    elsif($anc_dip_taxon eq 'vardousPindicola' && $desc_dip_taxon eq 'Dinaric'){ 
      $lineage_code = 'Dinaric' 
    }
    elsif($anc_dip_taxon eq 'vardousPindicola' && $desc_dip_taxon eq 'EBalkan'){ 
      $lineage_code = 'EBalkan' 
    }
    
    elsif($anc_dip_taxon eq 'Dinaric' && $desc_dip_taxon eq 'EBalkan'){ 
      $lineage_code = 'EBalkan' 
    }
    elsif($anc_dip_taxon eq 'EBalkan' && $desc_dip_taxon eq 'Dinaric'){ 
      $lineage_code = 'Dinaric' 
    }
  }

  return $lineage_code;
}

1;
