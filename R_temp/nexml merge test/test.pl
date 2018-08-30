#!/usr/bin/perl -w

use Bio::Phylo::IO 'parse';
use Bio::Phylo::NeXML::Meta;

my $data = do { local $/; <DATA> };

my $project = parse(
	'-format'  =>   'NeXML',
	'-string'  =>   $data,
	'-as_project'  => 1,
);

my $otus = $project->get_forests->[0]->make_taxa;

print $otus

