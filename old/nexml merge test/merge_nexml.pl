#!/usr/bin/perl -w

use Bio::Phylo::IO 'parse';

# load trees, get taxa from first tree
my $trees_proj = parse('-format' => 'newick', '-file' => $treefile, '-as_project' => 1 );
my $trees_taxa = $trees_proj->get_taxa->[0];

# load chars, get taxa from first char matrix
my $chars_proj = parse( '-format' => 'nexus', '-file' => $charfile, '-as_project' => 1 );
my $chars_taxa = $chars_proj->get_taxa->[0];

# merge these two sets (collapses matching names)
my $merged_taxa = $chars_taxa->merge_by_name($ trees_taxa);

# link the trees from $treefile to this merged set of OTUs
$trees_proj->get_ forests->[0]->set_taxa($ merged_taxa);

# link the first chars block from $charsfile to merged set
$chars_proj->get_ matrices->[0]->set_taxa($ merged_taxa);

# collect new merged taxa block and the trees into one project
$chars_proj->insert($ merged_taxa);
$chars_proj->insert($ trees_proj->get_forests->[0]);

# delete the old taxa block read from the nexus file
$chars_proj->delete($ chars_taxa);

# write to nexml
print $chars_proj->to_xml;