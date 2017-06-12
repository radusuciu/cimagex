"""Defines a PeptideDataset, which is a collection of Peptides.

Clases and functions here are for dealing with datasets that are grouped by peptide. Not the most elegant tbh but
I can't afford to spend the time so.. :/
"""

from copy import deepcopy
from .peptide import Peptide as MudpitPeptide
from .protein import Protein as MudpitProtein
from .parse_combined import ParseCombined
import csv
import itertools
import operator


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

    @property
    def clean_sequence(self):
        """Get totally clean peptide sequence."""
        return self.sequence.split('.')[1].translate(Peptide.delchars)

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

    def __init__(self, sequence, uniprot, description, peptides=[], mean=None, median=None, stdev=None):
        """Init peptide."""
        self._id = '{}_{}'.format(uniprot, sequence)
        self.uniprot = uniprot
        self.description = description
        self.sequence = sequence
        self.peptides = peptides
        self.mean = mean
        self.median = median
        self.stdev = stdev

        self._clean_id = '{}_{}'.format(uniprot, self.clean_sequence)

    @property
    def clean_sequence(self):
        """Get totally clean peptide sequence."""
        return self.sequence.translate(Peptide.delchars)

    def get_num_unique_peptides(self):
        """Get number unique peptides by sequence."""
        return len(set([x.clean_sequence for x in self.peptides]))

    def remove_half_tryptic(self):
        """Removes half tryptic peptides from container."""
        for p in self.peptides:
            if p.sequence[0] not in ['K','R','-']: # nterminus
                self.peptides.remove(p)
            if  p.sequence[-3] not in ['K','R'] and p.sequence[-1] != '-': # cterminus    
                self.peptides.remove(p)

    def make_clean(self):
        """Set sequence and id to clean variants."""
        self._id = self._clean_id
        self.sequence = self.clean_sequence
        return self

    def strip_diff_mods(self):
        """Remove all diff mods from container and child peptides."""
        self.make_clean()

        for peptide in self.peptides:
            peptide.sequence = peptide.clean_sequence

        return self

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

        current = self.get(sequence.sequence, sequence.uniprot)

        # if we already have a protein in the dataset with the same uniprot, then just add them together
        # otherwise add it to the list
        if current:
            self.sequences[self.sequences.index(current)] = current + sequence
        else:
            self.sequences.append(sequence)

    def get(self, sequence, uniprot):
        """Get from the dataset by sequence, uniprot pair."""
        for item in self.sequences:
            if item.sequence == sequence and item.uniprot == uniprot:
                return item
        return None

    def get_by_id(self, sequence_id):
        """Get from the dataset by id"""
        for item in self.sequences:
            if item._id == sequence_id:
                return item
        return None

    def get_by_sequence(self, sequence):
        """Get all items that match a sequence."""
        matches = [s for s in self.sequences if s.sequence == sequence]
        return matches

    def dedupe(self):
        """Ensure that we have only unique protein entries."""
        unique_sequences = self.get_unique_ids()
        # if all of our proteins are unique then don't do any work
        if len(unique_sequences) == len(self.sequences):
            return
        else:
            # assume that our datasets already contain unique uniprots
            # otherwise I need to write some code here. Throwing exception to be safe.
            raise NotImplementedError

    def get_index(self):
        """Return a dict of all the proteins in dataset with the uniprot id as key."""
        return {p._id: p for p in self.sequences}

    def get_unique_ids(self):
        """Get set of unique items contained in dataset."""
        unique = set()
        for s in self.sequences:
            unique.add(s._id)

        return unique

    def get_unique_clean_ids(self):
        """Get set of unique items contained in dataset if we ignore diff mods."""
        unique = set()
        for s in self.sequences:
            unique.add(s._clean_id)
        return unique

    def get_unique_sequences(self):
        """Get set of unique sequences contained in dataset."""
        unique = set()
        for protein in self.sequences:
            unique.add(protein.sequence)

        return unique

    def get_unique_clean_sequences(self):
        """Get set of unique clean sequences contained in dataset."""
        unique = set()
        for sequence in self.sequences:
            unique.add(sequence.clean_sequence)
        return unique

    def filter(self, filter_callback):
        """Given any callback, filter the list of proteins contained in a dataset."""
        self.sequences = list(filter(filter_callback, self.sequences))

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

    def apply_blacklist_filter(self, blacklist):
        """Throw away sequences that are in blacklist."""
        self.sequences = [s for s in self.sequences if s.sequence not in blacklist]

    def remove(self, el):
        """Remove a single element by sequence or by passing the whole PeptideContainer."""
        if type(el) == PeptideContainer:
            self.sequences.remove(el)
        elif type(el) == str:
            self.apply_blacklist_filter(el)

    def remove_oxidized_only(self, oxidized_symbol='+'):
        """Removes sequences that have no non-oxidized variants."""
        oxidized_sequences = [s for s in self.sequences if oxidized_symbol in s.sequence]

        for s in oxidized_sequences:
            non_oxidized_variants = self.get(s.clean_sequence, s.uniprot)
            if not non_oxidized_variants:
                self.remove(s)

    def remove_double_labeled(self, label_symbol='*'):
        """Removes peptides that are doubly labeled."""
        self.filter(lambda s: s.sequence.count(label_symbol) <= 1)

    def remove_half_tryptic(self):
        """Removes half tryptic peptides."""
        for s in self.sequences:
            s.remove_half_tryptic()

    def remove_empty(self):
        """Removes empty peptide containers."""
        for s in self.sequences:
            if not s.peptides:
                self.remove(s)

    def merge_to_clean_sequences(self):
        """Merge peptides with the same base sequence, irrespective of diff mods."""
        for _id in self.get_unique_clean_ids():
            sequences_matching = [s.make_clean() for s in self.sequences if s._clean_id == _id]
            # merge into first match if we have multiple matches
            if len(sequences_matching) > 1:
                sequences_matching[0] = sum(sequences_matching[1:], sequences_matching[0])
                for s in sequences_matching[1:]:
                    self.remove(s)

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
                    'uniprot', 'description', 'sequence', 'mean',
                    'median', 'stdev', 'n', 'stderr'
                ]

            writer.writerow(headers)

            for sequence in self.sequences:
                writer.writerow([
                    sequence.uniprot,
                    sequence.description,
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


def make_dataset(combined_dta_path):
    parser = ParseCombined()
    raw = parser.parse_file(str(combined_dta_path), type='peptide')
    sequences = []

    for item in raw:
        peptides = [make_peptide(peptide) for peptide in item['peptides']]

        # grouping by uniprot in case there are sequences that are assigned to multiple
        # uniprot ids. for assembling our dataset, it would be preferrable to treat these as
        # separate entities even if they are kept together by cimage
        grouped = itertools.groupby(peptides, operator.attrgetter('uniprot'))

        for uniprot, group in grouped:
            sequences.append(make_sequence(uniprot, item, list(group)))

    return PeptideDataset(sequences=sequences)


def make_peptide(raw_peptide):
    return Peptide(
        sequence=raw_peptide['sequence'],
        uniprot=raw_peptide['uniprot_id'],
        description=raw_peptide['description'],
        mass=float(raw_peptide['mass']),
        charge=int(raw_peptide['charge']),
        segment=int(raw_peptide['segment']),
        ratio=float(raw_peptide['mr']),
        rsquared=0
    )


def make_sequence(uniprot, item, peptides):
    return PeptideContainer(
        sequence=item['sequence'],
        uniprot=uniprot,
        description=peptides[0].description,
        peptides=peptides,
        mean=item['mean_ratio'], median=None, stdev=None
    )
