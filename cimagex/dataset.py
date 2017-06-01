"""Defines a Dataset, which is a collection of Proteins."""

from copy import deepcopy
from .protein import Protein
import csv


class Dataset():
    """Holds a list of proteins and defines methods for their manipulation."""

    def __init__(self, proteins=[], species=None):
        """Initialize with a list of proteins."""
        self.proteins = proteins
        self.species = species

    def add(self, protein):
        """Add a protein to the list."""
        if not type(protein) == Protein:
            return

        current = self.get(protein.uniprot)

        # if we already have a protein in the dataset with the same uniprot, then just add them together
        # otherwise add it to the list
        if current:
            self.proteins[self.proteins.index(current)] = current + protein
        else:
            self.proteins.append(protein)

    def get(self, uniprot):
        """Get a protein from the dataset by uniprot id."""
        for protein in self.proteins:
            if protein.uniprot == uniprot:
                return protein

        return None

    def dedupe(self):
        """Ensure that we have only unique protein entries."""
        unique_uniprots = self.get_unique_uniprots()
        # if all of our proteins are unique then don't do any work
        if len(unique_uniprots) == len(self.proteins):
            return
        else:
            # assume that our datasets already contain unique uniprots
            # otherwise I need to write some code here. Throwing exception to be safe.
            raise NotImplementedError

    def index_by_uniprot(self):
        """Return a dict of all the proteins in dataset with the uniprot id as key."""
        return {p.uniprot: p for p in self.proteins}

    def get_unique_uniprots(self):
        """Get set of unique uniprots contained in dataset."""
        unique = set()
        for protein in self.proteins:
            unique.add(protein.uniprot)
        return unique

    def filter(self, filter_callback):
        """Given any callback, filter the list of proteins contained in a dataset."""
        self.proteins = filter(filter_callback, self.proteins)

    def apply_rsquared_cutoff(self, cutoff):
        """Filter dataset by ratio cutoff."""
        filtered = []

        for protein in self.proteins:
            protein.apply_rsquared_cutoff(cutoff)
            if protein.peptides:
                filtered.append(protein)

        self.proteins = filtered

    def apply_unique_filter(self, cutoff):
        """Only keep proteins that have a minimum number of unique peptides."""
        self.proteins = [p for p in self.proteins if p.get_num_unique_peptides() >= cutoff]

    def apply_whitelist_filter(self, whitelist):
        """Only keep proteins that are in the passed whitelist."""
        self.proteins = [p for p in self.proteins if p.uniprot in whitelist]

    def generate_stats(self, ratio_filter=None):
        """Generate stats for each protein in dataset."""
        for protein in self.proteins:
            protein.generate_stats(ratio_filter)

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

            for protein in self.proteins:
                writer.writerow([
                    protein.uniprot,
                    protein.symbol,
                    protein.mean,
                    protein.median,
                    protein.stdev,
                    protein.n,
                    protein.stderr
                ])

    def __add__(self, dataset):
        """Two datasets can be added together by adding all their consituent proteins."""
        new_dataset = deepcopy(self)

        for protein in dataset.proteins:
            new_dataset.add(protein)
        return new_dataset

    def __radd__(self, *args, **kwargs):
        return self.__add__(*args, **kwargs)

    def __repr__(self):
        return 'Dataset(proteins={})'.format(self.proteins)
