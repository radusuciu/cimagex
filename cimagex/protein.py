"""Defines protein class."""

from copy import deepcopy
import statistics
import math


class Protein():
    """A protein is a collectin of peptides."""

    def __init__(self, uniprot, symbol, description, peptides=[], mean=None, median=None, stdev=None):
        """Init protein."""
        self.uniprot = uniprot
        self.symbol = symbol
        self.description = description
        self.peptides = peptides
        self.mean = mean
        self.median = median
        self.stdev = stdev

    def apply_rsquared_cutoff(self, cutoff):
        """Only keep peptides that meet a specific rsquared cutoff."""
        self.peptides = [peptide for peptide in self.peptides if peptide.rsquared >= cutoff]

    def get_num_unique_peptides(self):
        """Get number unique peptides by sequence."""
        return len(set([x.clean_sequence for x in self.peptides]))

    def generate_stats(self, ratio_filter=None, inverse=False):
        """Set mean, median and stdev."""
        raw_ratios = [p.ratio for p in self.peptides]

        if ratio_filter:
            ratios = ratio_filter(raw_ratios)
        else:
            ratios = [r for r in raw_ratios if r > 0]

        if inverse:
            ratios = [1/x for x in ratios if x > 0]

        self.mean = self.special_mean(ratios)
        self.median = self.special_median(ratios)
        self.stdev = self.special_stdev(ratios)
        self.stderr = self.special_stderr(ratios)
        self.n = len(ratios)

    @staticmethod
    def special_mean(ratios):
        """Return mean ratio given a set of ratios."""
        try:
            return statistics.mean(ratios)
        except statistics.StatisticsError:
            return '-'

    @staticmethod
    def special_median(ratios):
        """Return median ratio given a set of ratios."""
        try:
            return statistics.median(ratios)
        except statistics.StatisticsError:
            return '-'

    @staticmethod
    def special_stdev(ratios):
        """Return the standard deviation of a set of ratios."""
        try:
            return statistics.stdev(ratios)
        except statistics.StatisticsError:
            return '-'

    @staticmethod
    def special_stderr(ratios):
        """Return the standard error of a set of ratios."""
        try:
            return statistics.stdev(ratios) / math.sqrt(len(ratios))
        except:
            return '-'

    def __add__(self, protein):
        """Add one protein to another by concatenating peptides with the condition that they have the same uniprot."""
        if not type(self) == type(protein) or not protein.uniprot == self.uniprot:
            return

        new_protein = deepcopy(self)
        new_protein.peptides = self.peptides + protein.peptides

        return new_protein

    def __radd__(self, *args, **kwargs):
        print('radd called')
        return self.__add__(*args, **kwargs)

    def __repr__(self):
        return 'Protein(uniprot={}, symbol={}, description={}, peptides={}, mean={}, median={}, stdev={})'.format(
            self.uniprot,
            self.symbol,
            self.description,
            self.peptides,
            self.mean,
            self.median,
            self.stdev
        )
