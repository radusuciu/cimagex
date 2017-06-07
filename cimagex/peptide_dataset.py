"""Defines a PeptideDataset, which is a collection of Peptides.

Clases and functions here are for dealing with datasets that are grouped by peptide. Not the most elegant tbh but
I can't afford to spend the time so.. :/
"""

from copy import deepcopy
from .peptide import Peptide as MudpitPeptide
from .protein import Protein as MudpitProtein
import statistics
import math
import csv



class Peptide(MudpitPeptide):
    """An individual peptide."""

    def __init__(self, sequence, uniprot, description, mass, charge, segment, ratio, rsquared):
        """Init peptide."""
        self.sequence = sequence
        self.uniprot = uniprot
        self.description = description
        self.mass = mass
        self.charge = charge
        self.segment = segment
        self.ratio = ratio
        self.rsquared = rsquared

    def __repr__(self):
        return 'Peptide(sequence={}, uniprot={}, description={}, mass={}, charge={}, segment={}, ratio={}, rsquared={})'.format(
            self.sequence,
            self.uniprot,
            self.description,
            self.mass,
            self.charge,
            self.segment,
            self.ratio,
            self.rsquared
        )


class PeptideContainer(MudpitProtein):
    """A grouping of other peptides."""

    def __init__(self, sequence, peptides=[], mean=None, median=None, stdev=None):
        """Init peptide."""
        self.sequence = sequence
        self.peptides = peptides
        self.mean = mean
        self.median = median
        self.stdev = stdev

    def get_num_unique_peptides(self):
        """Get number unique peptides by sequence."""
        return len(set([x.clean_sequence for x in self.peptides]))

    def __add__(self, peptide):
        """Add one peptide to another by concatenating peptides with the condition that they have the same sequence."""
        if not type(self) == type(peptide) or not peptide.sequence == self.sequence:
            return

        new_sequence = deepcopy(self)
        new_sequence.peptides = self.peptides + peptide.peptides

        return new_sequence

    def __repr__(self):
        return 'PeptideContainer(sequence={}, peptides={}, mean={}, median={}, stdev={})'.format(
            self.sequence,
            self.peptides,
            self.mean,
            self.median,
            self.stdev
        )


class PeptideDataset():
    """Holds a list of peptides and defines methods for their manipulation."""

    def __init__(self, sequences=[]):
        """Initialize with a list of sequences (peptide containers)."""
        self.sequences = sequences

    def add(self, sequence):
        """Add a protein to the list."""
        if not type(sequence) == PeptideContainer:
            return

        current = self.get(sequence.sequence)

        # if we already have a protein in the dataset with the same uniprot, then just add them together
        # otherwise add it to the list
        if current:
            self.sequences[self.sequences.index(current)] = current + sequence
        else:
            self.sequences.append(sequence)

    def get(self, sequence):
        """Get from the dataset by sequence."""
        for item in self.sequences:
            if item.sequence == sequence:
                return item

        return None

    def dedupe(self):
        """Ensure that we have only unique protein entries."""
        unique_sequences = self.get_unique_sequences()
        # if all of our proteins are unique then don't do any work
        if len(unique_sequences) == len(self.sequences):
            return
        else:
            # assume that our datasets already contain unique uniprots
            # otherwise I need to write some code here. Throwing exception to be safe.
            raise NotImplementedError

    def index_by_sequence(self):
        """Return a dict of all the proteins in dataset with the uniprot id as key."""
        return {p.sequence: p for p in self.sequences}

    def get_unique_sequences(self):
        """Get set of unique uniprots contained in dataset."""
        unique = set()
        for protein in self.sequences:
            unique.add(protein.sequence)
        return unique

    def filter(self, filter_callback):
        """Given any callback, filter the list of proteins contained in a dataset."""
        self.sequences = filter(filter_callback, self.sequences)

    def apply_rsquared_cutoff(self, cutoff):
        """Filter dataset by ratio cutoff."""
        filtered = []

        for sequence in self.sequences:
            sequence.apply_rsquared_cutoff(cutoff)
            if sequence.peptides:
                filtered.append(sequence)

        self.sequences = filtered

    def apply_unique_filter(self, cutoff):
        """Only keep sequences that have a minimum number of unique peptides."""
        self.sequences = [p for p in self.sequences if p.get_num_unique_peptides() >= cutoff]

    def apply_whitelist_filter(self, whitelist):
        """Only keep sequences that are in the passed whitelist."""
        self.sequences = [p for p in self.sequences if p.sequence in whitelist]

    def generate_stats(self, ratio_filter=None):
        """Generate stats for each protein in dataset."""
        for sequence in self.sequences:
            sequence.generate_stats(ratio_filter)

    def to_csv(self, filename, headers=None):
        """Output dataset to .csv file."""
        with open(filename, 'w') as f:
            writer = csv.writer(f, lineterminator='\n')

            if not headers:
                headers = [
                    'uniprot', 'symbol', 'mean',
                    'median', 'stdev', 'n', 'stderr'
                ]

            writer.writerow(headers)

            for sequence in self.sequences:
                writer.writerow([
                    sequence.sequence,
                    sequence.mean,
                    sequence.median,
                    sequence.stdev,
                    sequence.n,
                    sequence.stderr
                ])

    def __add__(self, dataset):
        """Two datasets can be added together by adding all their consituent sequences."""
        new_dataset = deepcopy(self)

        for sequence in dataset.sequences:
            new_dataset.add(sequence)
        return new_dataset

    def __radd__(self, *args, **kwargs):
        return self.__add__(*args, **kwargs)

    def __repr__(self):
        return 'PeptideDataset(sequences={})'.format(self.sequences)
