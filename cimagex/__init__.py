from operator import itemgetter
from .dataset import Dataset
from .peptide_dataset import PeptideDataset
from .peptide import Peptide
from .protein import Protein
from .parse_combined import ParseCombined
from .parse_combined import make_dataset as make_dataset_from_combined
import itertools


def make_dataset(headers, data):
    """Given raw data, create a Dataset."""
    dataset = Dataset()

    # sort and group by uniprot
    data.sort(key=itemgetter(headers.index('uniprot')))
    grouped = itertools.groupby(data, key=itemgetter(headers.index('uniprot')))

    proteins = []

    # assemble into a Dataset made up of Proteins and Peptides
    for uniprot, g in grouped:
        group = list(g)
        peptides = [make_peptide(headers, item) for item in group]
        proteins.append(make_protein(headers, group[0], peptides))

    dataset.proteins = proteins
    dataset.dedupe()

    return dataset

