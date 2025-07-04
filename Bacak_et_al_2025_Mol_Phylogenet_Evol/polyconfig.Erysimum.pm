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
our @diploids = qw( carniolicum crassistylum cuspidatum kuemmerlei linariifolium malcolmia pectinatum saxosum sylvestre vitekii witmannii microstylum ); 


# note these are diploids as well; hash instead of list
our %outgroups = ( 
  'malcolmia',1
);

# diploid used to root trees 
our $ROOT = 'malcolmia';

# Abbreviated names of polyploid species as found in FASTA and tree files.
our @polyploids = ('BREITs1.h1_odoratum4x', 'BREITs1.h2_odoratum4x', 'BREITs1.h3_odoratum4x', 'BREITs1.h4_odoratum4x', 'CIUC6R.h1_odoratum22chrom', 'CIUC6R.h2_odoratum22chrom', 'CIUC6R.h3_odoratum22chrom', 'CIUC6R.h4_odoratum22chrom', 'GLAVs1.h1_odoratum6x', 'GLAVs1.h2_odoratum6x', 'GLAVs1.h3_odoratum6x', 'GLAVs1.h4_odoratum6x', 'GLAVs1.h5_odoratum6x', 'GLAVs1.h6_odoratum6x', 'Kavlar2.h1_odoratum4x', 'Kavlar2.h2_odoratum4x', 'Kavlar2.h3_odoratum4x', 'Kavlar2.h4_odoratum4x', 'Maliscak9.h1_odoratum22chrom', 'Maliscak9.h2_odoratum22chrom', 'Maliscak9.h3_odoratum22chrom', 'Maliscak9.h4_odoratum22chrom', 'RET1s.h1_odoratumRetezat', 'RET1s.h2_odoratumRetezat', 'RET1s.h3_odoratumRetezat', 'RET1s.h4_odoratumRetezat', 'SCO3.h1_odoratumRetezat', 'SCO3.h2_odoratumRetezat', 'SCO3.h3_odoratumRetezat', 'SCO3.h4_odoratumRetezat', 'VGrad1.h1_odoratum6x', 'VGrad1.h2_odoratum6x', 'VGrad1.h3_odoratum6x', 'VGrad1.h4_odoratum6x', 'VGrad1.h5_odoratum6x', 'VGrad1.h6_odoratum6x');

# Abbreviated names of polyploid subgenomes, defined by user.
# These are used to annotate alleles from the subgenomes as ancestral.
# Leave empty if not required
our @subgenomes = ();


# Abbreviated names of labelled polyploid sequences as found in FASTA and tree files.
# Note that the capital letters correpond to @CODES below
our @polyploids_labelled = ('BREITs1.h1_odoratum4x', 'BREITs1.h2_odoratum4x', 'BREITs1.h3_odoratum4x', 'BREITs1.h4_odoratum4x', 'CIUC6R.h1_odoratum22chrom', 'CIUC6R.h2_odoratum22chrom', 'CIUC6R.h3_odoratum22chrom', 'CIUC6R.h4_odoratum22chrom', 'GLAVs1.h1_odoratum6x', 'GLAVs1.h2_odoratum6x', 'GLAVs1.h3_odoratum6x', 'GLAVs1.h4_odoratum6x', 'GLAVs1.h5_odoratum6x', 'GLAVs1.h6_odoratum6x', 'Kavlar2.h1_odoratum4x', 'Kavlar2.h2_odoratum4x', 'Kavlar2.h3_odoratum4x', 'Kavlar2.h4_odoratum4x', 'Maliscak9.h1_odoratum22chrom', 'Maliscak9.h2_odoratum22chrom', 'Maliscak9.h3_odoratum22chrom', 'Maliscak9.h4_odoratum22chrom', 'RET1s.h1_odoratumRetezat', 'RET1s.h2_odoratumRetezat', 'RET1s.h3_odoratumRetezat', 'RET1s.h4_odoratumRetezat', 'SCO3.h1_odoratumRetezat', 'SCO3.h2_odoratumRetezat', 'SCO3.h3_odoratumRetezat', 'SCO3.h4_odoratumRetezat', 'VGrad1.h1_odoratum6x', 'VGrad1.h2_odoratum6x', 'VGrad1.h3_odoratum6x', 'VGrad1.h4_odoratum6x', 'VGrad1.h5_odoratum6x', 'VGrad1.h6_odoratum6x');
# our @polyploids_labelled = (
#    'Bmex_A','Bmex_B','Bmex_C','Bmex_D','Bmex_E','Bmex_F','Bmex_G','Bmex_H','Bmex_I',
#    'Bboi_A','Bboi_B','Bboi_C','Bboi_D','Bboi_E','Bboi_F','Bboi_G','Bboi_H','Bboi_I',
#    'Bret_A','Bret_B','Bret_C','Bret_D','Bret_E','Bret_F','Bret_G','Bret_H','Bret_I',
#    'Bhyb_A','Bhyb_B','Bhyb_C','Bhyb_D','Bhyb_E','Bhyb_F','Bhyb_G','Bhyb_H','Bhyb_I',
#    'Brup_A','Brup_B','Brup_C','Brup_D','Brup_E','Brup_F','Brup_G','Brup_H','Brup_I',
#    'Bpho_A','Bpho_B','Bpho_C','Bpho_D','Bpho_E','Bpho_F','Bpho_G','Bpho_H','Bpho_I',
#    'B422_A','B422_B','B422_C','B422_D','B422_E','B422_F','B422_G','B422_H','B422_I'
# );

# Optional custom definition of clades that contain >1 diploids, if any.
# Leave empty or comment out the examples otherwise.
# MRCA nodes are internal node names, can be used in rules below; there are two keys:
# # i) MRCA for all species in all sister clades, also called a bifurcation in the code
# # ii) MRCA for each explicitely defined clade (there should be two)
# # finally the arrays contain lists of species in each clade, should be diploid
our %sister_clades = ();
# $sister_clades{'MRCA_BiHSrbMn_subspAcris'}{'MRCA1'} = ['BiHSrbMn','subspAcris']; # clade 1
# $sister_clades{'MRCA_BiHSrbMn_subspAcris'}{'MRCA2'} = ['BiHSrbMn','subspAcris']; # clade1 again with a different name, hack



# Default values for paremeters controlling how anatolicaid block covered by outgroups & polyploids

# Optional user-defined contribution of diploid species to block width calculations
# in script _trim_MSA_block.pl (outgroup species are not used)
# Can be used to indicate that only a member of a clade is required, see Bsyl example 
#my %diploids4width = (
# Example wheat values
#   'Tura'=>1,
#	'Tmon'=>1,
#	'Asha'=>1,
#	'Atau'=>1,
#	'Aspe'=>1,
# Brachypodium values
#'Bsta'=>1,
#'Bdis'=>1,
#'Barb'=>1,
#'Bpin'=>1,
#'Bsyl'=>0.2,
#);


# see rules defined below
our @CODES = qw( carniolicum crassistylum cuspidatum kuemmerlei linariifolium malcolmia pectinatum saxosum sylvestre vitekii witmannii microstylum all ); 

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
    elsif($anc_dip_taxon eq 'Bdis'){ $lineage_code = 'D' }

    elsif($anc_dip_taxon eq 'carniolicum'){ $lineage_code = 'carniolicum' }
    elsif($anc_dip_taxon eq 'crassistylum'){ $lineage_code = 'crassistylum' }
    elsif($anc_dip_taxon eq 'cuspidatum'){ $lineage_code = 'cuspidatum' }
    elsif($anc_dip_taxon eq 'kuemmerlei'){ $lineage_code = 'kuemmerlei' }
    elsif($anc_dip_taxon eq 'linariifolium'){ $lineage_code = 'linariifolium' }
    elsif($anc_dip_taxon eq 'malcolmia'){ $lineage_code = 'malcolmia' }
    elsif($anc_dip_taxon eq 'pectinatum'){ $lineage_code = 'pectinatum' }
    elsif($anc_dip_taxon eq 'saxosum'){ $lineage_code = 'saxosum' }
    elsif($anc_dip_taxon eq 'sylvestre'){ $lineage_code = 'sylvestre' }
    elsif($anc_dip_taxon eq 'vitekii'){ $lineage_code = 'vitekii' }
    elsif($anc_dip_taxon eq 'witmannii'){ $lineage_code = 'witmannii' }
    elsif($anc_dip_taxon eq 'microstylum'){ $lineage_code = 'microstylum' }


    print "anc_is_sister == 1\n";
  }
  else { ## both ancestor/descendant diploids/clades are considered

    if($anc_dip_taxon eq 'Hvul' && $desc_dip_taxon eq 'Bsta'){ 
      $lineage_code = 'A' 
    }




    elsif($anc_dip_taxon eq 'malcolmia' && $desc_dip_taxon eq 'cuspidatum'){ 
      $lineage_code = 'cuspidatum' 
    }
    elsif($anc_dip_taxon eq 'cuspidatum' && $desc_dip_taxon eq 'witmannii'){ 
      $lineage_code = 'witmannii' 
    }
    elsif($anc_dip_taxon eq 'witmannii' && $desc_dip_taxon eq 'pectinatum'){ 
      $lineage_code = 'pectinatum' 
    }
    elsif($anc_dip_taxon eq 'witmannii' && $desc_dip_taxon eq 'carniolicum'){ 
      $lineage_code = 'carniolicum' 
    }
    elsif($anc_dip_taxon eq 'witmannii' && $desc_dip_taxon eq 'saxosum'){ 
      $lineage_code = 'saxosum' 
    }
    elsif($anc_dip_taxon eq 'witmannii' && $desc_dip_taxon eq 'sylvestre'){ 
      $lineage_code = 'sylvestre' 
    }

    elsif($anc_dip_taxon eq 'pectinatum' && $desc_dip_taxon eq 'crassistylum'){ 
      $lineage_code = 'crassistylum' 
    }
    elsif($anc_dip_taxon eq 'crassistylum' && $desc_dip_taxon eq 'microstylum'){ 
      $lineage_code = 'microstylum' 
    }
    elsif($anc_dip_taxon eq 'crassistylum' && $desc_dip_taxon eq 'linariifolium'){ 
      $lineage_code = 'linariifolium' 
    }
    elsif($anc_dip_taxon eq 'linariifolium' && $desc_dip_taxon eq 'microstylum'){ 
      $lineage_code = 'microstylum' 
    }
    elsif($anc_dip_taxon eq 'microstylum' && $desc_dip_taxon eq 'linariifolium'){ 
      $lineage_code = 'linariifolium' 
    }

    elsif($anc_dip_taxon eq 'carniolicum' && $desc_dip_taxon eq 'vitekii'){ 
      $lineage_code = 'vitekii' 
    }
    elsif($anc_dip_taxon eq 'carniolicum' && $desc_dip_taxon eq 'kuemmerlei'){ 
      $lineage_code = 'kuemmerlei' 
    }
    elsif($anc_dip_taxon eq 'vitekii' && $desc_dip_taxon eq 'kuemmerlei'){ 
      $lineage_code = 'kuemmerlei' 
    }
    elsif($anc_dip_taxon eq 'kuemmerlei' && $desc_dip_taxon eq 'vitekii'){ 
      $lineage_code = 'vitekii' 
    }

    elsif($anc_dip_taxon eq 'saxosum' && $desc_dip_taxon eq 'sylvestre'){ 
      $lineage_code = 'sylvestre' 
    }
    elsif($anc_dip_taxon eq 'sylvestre' && $desc_dip_taxon eq 'saxosum'){ 
      $lineage_code = 'saxosum' 
    }

    
    # elsif($anc_dip_taxon eq 'vardousPindicola' && $desc_dip_taxon eq 'MRCA_BiHSrbMn_subspAcris'){ 
    #   $lineage_code = 'MRCA_BiHSrbMn_subspAcris_anc' 
    # }
    # elsif($anc_dip_taxon eq 'MRCA_BiHSrbMn_subspAcris' && $desc_dip_taxon eq 'subspAcris'){ 
    #   $lineage_code = 'subspAcris' 
    # }
    # elsif($anc_dip_taxon eq 'MRCA_BiHSrbMn_subspAcris' && $desc_dip_taxon eq 'BiHSrbMn'){ 
    #   $lineage_code = 'BiHSrbMn' 
    # }

    elsif($anc_dip_taxon eq 'BiHSrbMn' && $desc_dip_taxon eq 'subspAcris'){ 
      $lineage_code = 'subspAcris' 
    }
    elsif($anc_dip_taxon eq 'subspAcris' && $desc_dip_taxon eq 'BiHSrbMn'){ 
      $lineage_code = 'BiHSrbMn' 
    }

    elsif($anc_dip_taxon eq 'Bsta' && $desc_dip_taxon eq 'Bdis'){ 
      $lineage_code = 'C' 
    }
    elsif($anc_dip_taxon eq 'Bdis' && $desc_dip_taxon eq 'Barb'){ 
      $lineage_code = 'E' 
    }
    elsif($anc_dip_taxon eq 'Barb' && $desc_dip_taxon eq 'MRCABsylBpin'){
      $lineage_code = 'G'
    }
    elsif($anc_dip_taxon eq 'MRCABsylBpin' && $desc_dip_taxon eq 'Bsyl'){
      $lineage_code = 'H'
    }
    elsif($anc_dip_taxon eq 'MRCABsylBpin' && $desc_dip_taxon eq 'Bpin'){
      $lineage_code = 'I'
    }
    elsif($anc_dip_taxon eq 'Bsyl' && $desc_dip_taxon eq 'Bpin'){
      $lineage_code = 'I'
    }
    elsif($anc_dip_taxon eq 'Bpin' && $desc_dip_taxon eq 'Bsyl'){
      $lineage_code = 'H'
    }
    # swapped Barb and Bsyl/Bpin
    elsif($anc_dip_taxon eq 'MRCABsylBpin' || $anc_dip_taxon eq 'Bpin'
          || $anc_dip_taxon eq 'Bsyl'){
      if($desc_dip_taxon eq 'Barb'){
        $lineage_code = 'F'
      }
    }
    elsif($anc_dip_taxon eq 'Bdis'){
      if($desc_dip_taxon eq 'MRCABsylBpin' || $desc_dip_taxon eq 'Bpin' 
         || $desc_dip_taxon eq 'Bsyl'){
        $lineage_code = 'E'
      }
    }
  }

  return $lineage_code;
}

1;
