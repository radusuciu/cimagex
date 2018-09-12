"""Defines a Dataset, which is a collection of Proteins."""
from copy import deepcopy
from collections import OrderedDict
from .protein import Protein
from .peptide import Peptide
from .parse_combined import ParseCombined
import csv
import uuid

class Dataset():
    """Holds a list of proteins and defines methods for their manipulation."""

    def __init__(self, proteins=[], species=None, name=None, inhibitor=None, concentration=None):
        """Initialize with a list of proteins."""
        self.uuid = uuid.uuid4()

        for p in proteins:
            if not p.uuid:
                p.uuid = self.uuid

        self.proteins = proteins
        self.species = species
        self.name = name
        self.inhibitor = inhibitor
        self.concentration = concentration

    def add(self, protein):
        """Add a protein to the list."""
        if not type(protein) == Protein:
            return

        current = self.get(protein.uniprot)

        if not protein.uuid:
            protein.uuid = self.uuid

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

    def get_unique_uuids(self):
        """Get set of all unique uuids which contributed to this dataset."""
        # return(set(p.uuid for p in s.peptides for s in self.sequences))
        return(set(p.uuid for s in self.proteins for p in s.peptides))

    def filter(self, filter_callback):
        """Given any callback, filter the list of proteins contained in a dataset."""
        self.proteins = filter(filter_callback, self.proteins)

    def apply_rsquared_cutoff(self, cutoff):
        """Filter dataset by rsquared cutoff."""
        filtered = []

        for protein in self.proteins:
            protein.apply_rsquared_cutoff(cutoff)
            if protein.peptides:
                filtered.append(protein)

        self.proteins = filtered

    def apply_unique_filter(self, cutoff):
        """Only keep proteins that have a minimum number of unique peptides."""
        self.proteins = [p for p in self.proteins if p.get_num_unique_peptides() >= cutoff]

    def apply_unique_quantified_filter(self, cutoff):
        """Only keep proteins that have a minimum number of unique peptides."""
        self.proteins = [p for p in self.proteins if p.get_num_unique_quantified_peptides() >= cutoff]

    def apply_datasets_quantified_filter(self, cutoff):
        self.proteins = [p for p in self.proteins if p.get_num_datasets_quantified() >= cutoff]

    def apply_whitelist_filter(self, whitelist):
        """Only keep proteins that are in the passed whitelist."""
        self.proteins = [p for p in self.proteins if p.uniprot in whitelist]

    def apply_symbol_blacklist_filter(self, whitelist):
        """Only keep proteins that are in the passed whitelist."""
        self.proteins = [p for p in self.proteins if p.symbol not in whitelist]

    def apply_keratin_filter(self):
        self.proteins = [p for p in self.proteins if 'KRT' != p.symbol[:3]]

    def apply_blacklist_filter(self, blacklist):
        """Throw away sequences that are in blacklist."""
        self.proteins = [p for p in self.proteins if p.uniprot not in blacklist]

    def remove(self, el):
        """Remove a single element by uniprot or by passing the whole Protein."""
        if type(el) == Protein:
            self.proteins.remove(el)
        elif type(el) == str:
            self.apply_blacklist_filter([el])

    def remove_reverse_matches(self):
        """Removes Reverse Uniprot sequences found by IP2."""
        self.proteins = [p for p in self.proteins if 'Reverse' not in p.uniprot]

    def remove_half_tryptic(self):
        """Removes half tryptic peptides."""
        for p in self.proteins:
            p.remove_half_tryptic()

    def remove_oxidized_only(self, oxidized_symbol='+'):
        """Removes sequences that have no non-oxidized variants."""
        for p in self.proteins:
            p.remove_oxidized_only(oxidized_symbol)

    def remove_oxidized_proteins(self, oxidized_symbol='+'):
        for p in self.proteins:
            if all(oxidized_symbol in x.sequence for x in p.peptides):
                self.remove(p)

    def remove_oxidized_methionines(self, oxidized_symbol='+'):
        for p in self.proteins:
            p.remove_oxidized_methionines(oxidized_symbol)

    def remove_empty(self):
        """Removes proteins with no peptides."""
        self.proteins = [p for p in self.proteins if p.peptides]

    def generate_stats(self, ratio_filter=None, inverse=False):
        """Generate stats for each protein in dataset."""
        for protein in self.proteins:
            protein.generate_stats(ratio_filter, inverse, replicate_medians_dict=OrderedDict.fromkeys(self.get_unique_uuids()))
        return self

    def filter_20s(self, ratio_cutoff=4):
        """Filter erroneous 20s from data."""
        for protein in self.proteins:
            # protein.filter_20s_by_ms2()
            protein.filter_20s(ratio_cutoff)
        return self

    def filter_by_stdev(self, stdev_cutoff=0.6, ratio_cutoff=4):
        """Filter peptides that have high standard deviations."""
        for protein in self.proteins:
            protein.filter_by_stdev(stdev_cutoff, ratio_cutoff)
        return self

    def to_csv(self, filename, headers=None):
        """Output dataset to .csv file."""
        with open(str(filename), 'w') as f:
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


def make_dataset(combined_dta_path, name=None, parse_as_file=True):
    parser = ParseCombined(parser_type='protein')

    if parse_as_file:
        raw = parser.parse_file(str(combined_dta_path))
    else:
        raw = parser.parse_io(combined_dta_path)

    dataset = Dataset(name=name)
    _uuid = dataset.uuid
    proteins = []

    for protein in raw:
        peptides = [make_peptide(peptide) for peptide in protein['peptides']]
        proteins.append(make_protein(protein, peptides, uuid=_uuid))

    dataset.proteins = proteins

    return dataset


def make_peptide(raw_peptide):
    """Given a list of attributes and appropriate headers, make a Peptide."""
    return Peptide(
        sequence=raw_peptide['sequence'],
        mass=float(raw_peptide['mass']),
        charge=int(raw_peptide['charge']),
        segment=int(raw_peptide['segment']),
        ratio=float(raw_peptide['mr']),
        rsquared=0
    )


def make_protein(raw_protein, peptides, uuid=None):
    """Given a list of attributes and appropriate headers, make a Protein."""
    return Protein(
        uniprot=raw_protein['uniprot_id'],
        symbol=raw_protein['symbol'],
        description=raw_protein['description'],
        peptides=peptides,
        mean=raw_protein['mean_ratio'], median=None, stdev=None,
        uuid=uuid
    )
