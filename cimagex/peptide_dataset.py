"""Defines a PeptideDataset, which is a collection of Peptides.

Clases and functions here are for dealing with datasets that are grouped by peptide. Not the most elegant tbh but
I can't afford to spend the time so.. :/
"""

from copy import deepcopy
from collections import defaultdict
from Bio import SeqIO
from .peptide import Peptide as MudpitPeptide
from .protein import Protein as MudpitProtein
from .parse_combined import ParseCombined
import csv
import itertools
import operator
import uuid



class Peptide(MudpitPeptide):
    """An individual peptide."""

    def __init__(self, sequence, uniprot, symbol, description, mass, charge, segment, ratio, rsquared, num_ms2=None, residue=None):
        """Init peptide."""
        self.sequence = sequence
        self.uniprot = uniprot
        self.symbol = symbol
        self.description = description
        self.mass = mass
        self.charge = charge
        self.segment = segment
        self.ratio = ratio
        self.rsquared = rsquared
        self.residue = residue

        if num_ms2:
            self.num_ms2 = num_ms2

    @property
    def clean_sequence(self):
        """Get totally clean peptide sequence."""
        return self.sequence.split('.')[1].translate(Peptide.delchars)

    def __repr__(self):
        return 'Peptide(sequence={}, uniprot={}, description={}, mass={}, charge={}, segment={}, ratio={}, rsquared={}, num_ms2={})'.format(
            self.sequence,
            self.uniprot,
            self.description,
            self.mass,
            self.charge,
            self.segment,
            self.ratio,
            self.rsquared,
            self.num_ms2
        )


class PeptideContainer(MudpitProtein):
    """A grouping of other peptides."""

    def __init__(self, sequence, uniprot, symbol, description, peptides=[], mean=None, median=None, stdev=None, uuid=None):
        """Init peptide."""
        self._id = '{}_{}'.format(uniprot, sequence)
        self.uniprot = uniprot
        self.symbol = symbol
        self.description = description
        self.sequence = sequence
        self.peptides = peptides
        self.mean = mean
        self.median = median
        self.stdev = stdev

        self.uuid = uuid

        self._clean_id = '{}_{}'.format(uniprot, self.clean_sequence)

    @property
    def uuid(self):
        return self._uuid

    @uuid.setter
    def uuid(self, uuid):
        for p in self.peptides:
            p.uuid = uuid

        self._uuid = uuid

    @property
    def clean_sequence(self):
        """Get totally clean peptide sequence."""
        return self.sequence.translate(Peptide.delchars)

    def get_num_unique_peptides(self):
        """Get number unique peptides by sequence."""
        return len(set([x.clean_sequence for x in self.peptides]))

    def filter_by_ms2(self, threshold=2):
        """Filter out peptides that don't have a certain number of ms2s."""
        for p in self.peptides:
            if p.num_ms2 < threshold:
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

    def __init__(self, sequences=[], url=None, name=None):
        """Initialize with a list of sequences (peptide containers)."""
        self.uuid = uuid.uuid4()

        for s in sequences:
            if not s.uuid:
                s.uuid = self.uuid

        self.sequences = sequences
        self.url = url
        self.name = name


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

    def standard_processing(self):
        """Convenience method for running common processing steps."""
        self.remove_reverse_matches()
        self.remove_half_tryptic()
        self.remove_double_labeled()
        self.remove_oxidized_only()
        self.remove_empty()
        self.merge_to_clean_sequences()
        self.remove_empty()
        self.generate_stats()

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

    def apply_ms2_filter(self, cutoff=1):
        """Filter dataset to only keep peptides that have at least a certain number of ms2s."""
        for s in self.sequences:
            # apply to peptides that only have one sequence
            if len(s.peptides) == 1:
                s.filter_by_ms2()
                continue

        return self

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
            self.apply_blacklist_filter([el])

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

    def remove_reverse_matches(self):
        """Removes Reverse Uniprot sequences found by IP2."""
        self.sequences = [s for s in self.sequences if 'Reverse' not in s.uniprot]

    def remove_only_zeros(self):
        """Remove peptide containers which contain only ratios of zero."""
        self.sequences = [s for s in self.sequences if any(p.ratio for p in s.peptides)]

    def remove_empty(self):
        """Removes empty peptide containers."""
        self.sequences = [s for s in self.sequences if s.peptides]

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

    def filter_by_stdev(self, stdev_cutoff=0.6, ratio_cutoff=4):
        """Filter peptides that have high standard deviations."""
        for sequence in self.sequences:
            sequence.filter_by_stdev(stdev_cutoff, ratio_cutoff)
        return self

    def filter_20s(self, ratio_cutoff=4):
        """Filter erroneous 20s from data."""
        for sequence in self.sequences:
            sequence.filter_20s(ratio_cutoff)
        return self

    def filter_20s_by_ms2(self):
        """Filter 20s that only have one ms2."""
        for sequence in self.sequences:
            sequence.filter_20s_by_ms2()

    def annotate_residues(self):
        """Annotate residues."""
        db = list(SeqIO.parse('data/db.txt', 'fasta'))
        uniprot_db = {item.id.split('|')[1]: item for item in db if 'Reverse' not in item.id}

        delchars = {ord(c): None for c in map(chr, range(256)) if not c.isalpha() and not c == '*'}

        for sequence in self.sequences:
            protein_sequence = uniprot_db[sequence.uniprot].seq

            for peptide in sequence.peptides:
                sequence = peptide.sequence.split('.')[1].translate(delchars)
                position_in_sequence = sequence.index('*')
                sequence = sequence.replace('*', '')

                sequence_position_in_protein = protein_sequence.find(sequence)
                peptide.residue = sequence_position_in_protein + position_in_sequence

    def to_csv(self, filename, headers=None):
        """Output dataset to .csv file."""
        with open(str(filename), 'w') as f:
            writer = csv.writer(f, lineterminator='\n')

            if not headers:
                headers = [
                    'id', 'uniprot', 'description', 'sequence', 'mean_of_medians',
                    'stdev', 'n'
                ]

            writer.writerow(headers)

            for sequence in self.sequences:
                writer.writerow([
                    sequence._id,
                    sequence.uniprot,
                    sequence.description,
                    sequence.sequence,
                    sequence.mean_of_medians,
                    sequence.stdev_mean_of_medians,
                    sequence.n
                ])

    def __add__(self, dataset):
        """Two datasets can be added together by adding all their constituent sequences."""
        new_dataset = deepcopy(self)

        for sequence in dataset.sequences:
            new_dataset.add(sequence)
        return new_dataset

    def __radd__(self, *args, **kwargs):
        return self.__add__(*args, **kwargs)

    def __repr__(self):
        return 'PeptideDataset(sequences={})'.format(self.sequences)


def make_dataset(combined_dta_path, dtaselect_path=None, name=None):
    parser = ParseCombined()
    raw = parser.parse_file(str(combined_dta_path), type='peptide')
    sequences = []

    if dtaselect_path:
        ms2_counts = parse_dtaselect(dtaselect_path)

    for item in raw:
        peptides = []

        for raw_peptide in item['peptides']:
            peptide = make_peptide(raw_peptide)

            if dtaselect_path:
                peptide.num_ms2 = ms2_counts.get(peptide.sequence)

            peptides.append(peptide)

        # grouping by uniprot in case there are sequences that are assigned to multiple
        # uniprot ids. for assembling our dataset, it would be preferrable to treat these as
        # separate entities even if they are kept together by cimage
        grouped = itertools.groupby(peptides, operator.attrgetter('uniprot'))

        dataset = PeptideDataset(name=name)
        uuid = dataset.uuid

        for uniprot, group in grouped:
            sequences.append(make_sequence(uniprot, item, list(group), uuid=uuid))

        dataset.sequences = sequences

    return dataset


def make_peptide(raw_peptide):
    return Peptide(
        sequence=raw_peptide['sequence'],
        uniprot=raw_peptide['uniprot_id'],
        description=raw_peptide['description'],
        symbol=raw_peptide['symbol'],
        mass=float(raw_peptide['mass']),
        charge=int(raw_peptide['charge']),
        segment=int(raw_peptide['segment']),
        ratio=float(raw_peptide['mr']),
        rsquared=0
    )


def make_sequence(uniprot, item, peptides, uuid=None):
    return PeptideContainer(
        sequence=item['sequence'],
        uniprot=uniprot,
        description=peptides[0].description,
        symbol=peptides[0].symbol,
        peptides=peptides,
        mean=item['mean_ratio'], median=None, stdev=None,
        uuid=uuid
    )


def parse_dtaselect(dtaselect_path):
    with open(str(dtaselect_path)) as f:
        raw = f.readlines()

    data = raw[29:]
    counts = defaultdict(int)

    for row in data:
        if row[0] in ('*', '\t'):
            sequence = row.strip().split('\t')[-1]
            counts[sequence] = counts[sequence] + 1

    return counts
