"""Defines protein class."""

from copy import deepcopy
from .utils import get_modified_z_scores, get_iqr_bounds
import statistics
import itertools
import numbers
import operator
import math


class Protein():
    """A protein is a collectin of peptides."""

    def __init__(self, uniprot, symbol, description, peptides=[], mean=None, median=None, stdev=None, uuid=None):
        """Init protein."""
        self.uniprot = uniprot
        self.symbol = symbol
        self.description = description
        self.mean = mean
        self.median = median
        self.stdev = stdev
        self.peptides = peptides
        self.uuid = uuid

    @property
    def uuid(self):
        return self._uuid

    @uuid.setter
    def uuid(self, uuid):
        for p in self.peptides:
            p.uuid = uuid

        self._uuid = uuid

    def copy_with_new_peptides(self, peptides=[]):
        """Makes a copy of self, and removes peptides."""
        me = deepcopy(self)
        me.peptides = peptides
        return me

    def apply_rsquared_cutoff(self, cutoff):
        """Only keep peptides that meet a specific rsquared cutoff."""
        self.peptides = [peptide for peptide in self.peptides if peptide.rsquared >= cutoff]

    def remove_half_tryptic(self):
        """Removes half tryptic peptides from protein."""
        self.peptides = [
            p for p in self.peptides 
            if p.sequence[0] in ['K','R','-']
            and (p.sequence[-3] in ['K','R'] or p.sequence[-1] == '-')
        ]

    def remove_oxidized_only(self, oxidized_symbol='+'):
        """Removes sequences that have no non-oxidized variants."""
        oxidized_peptides = [p for p in self.peptides if oxidized_symbol in p.sequence]

        for ox in oxidized_peptides:
            non_oxidized_variant = next((p for p in self.peptides if p.sequence == ox.clean_sequence), None)
            if not non_oxidized_variant:
                self.peptides.remove(ox)

    def remove_oxidized_methionines(self, oxidized_symbol='+'):
        self.peptides = [p for p in self.peptides if oxidized_symbol not in p.sequence]

    def get_num_datasets_quantified(self):
        return len(set(p.uuid for p in self.peptides if p.ratio > 0))

    def get_num_unique_peptides(self):
        """Get number unique peptides by sequence."""
        return len(set([x.clean_sequence for x in self.peptides]))

    def get_num_quantified_peptides(self):
        """Get number of quantified peptides."""
        return len([x.clean_sequence for x in self.peptides if x.ratio > 0])

    def get_num_unique_quantified_peptides(self):
        """Get number of unique quantified peptides."""
        return len(set([x.clean_sequence for x in self.peptides if x.ratio > 0]))

    def generate_stats(self, ratio_filter=None, inverse=False, replicate_medians_dict={}):
        """Set mean, median and stdev."""
        raw_ratios = [p.ratio for p in self.peptides]

        if ratio_filter:
            ratios = ratio_filter(raw_ratios)
        else:
            ratios = [r for r in raw_ratios if r > 0]

        if inverse:
            ratios = [1/x for x in ratios if x > 0]

        self.ratios = ratios

        self.mean = self.special_mean(ratios)
        self.median = self.special_median(ratios)

        replicate_medians = []
        grouped = itertools.groupby(self.peptides, operator.attrgetter('uuid'))

        for uuid, peptides in grouped:
            median = self.special_median(p.ratio for p in peptides if p.ratio > 0)
            # median = self.special_median(min(p.ratio, 10) for p in peptides if p.ratio > 0)

            if inverse and isinstance(median, numbers.Number) and median > 0:
                median = 1 / median

            replicate_medians_dict[uuid] = median
            replicate_medians.append(median)

        self.replicate_medians_dict = replicate_medians_dict

        # filtering out "-" returned from special median function
        # at least at time of writing, should probably not return non-number values
        replicate_medians  = [r for r in replicate_medians if isinstance(r, numbers.Number)]
        self.replicate_medians = replicate_medians

        self.mean_of_medians = self.special_mean(replicate_medians)
        self.median_of_medians = self.special_median(replicate_medians)

        self.stdev = self.special_stdev(ratios)
        self.stdev_mean_of_medians = self.special_stdev(replicate_medians)

        self.stderr = self.special_stderr(ratios, default=0)
        self.n = len(ratios)
        self.num_unique_peptides = len(set(p.clean_sequence for p in self.peptides))
        self.num_unique_quantified_peptides = len(set(p.clean_sequence for p in self.peptides if p.ratio > 0))
        self.num_datasets = len(set(p.uuid for p in self.peptides))
        self.num_datasets_quantified = self.get_num_datasets_quantified()


    def filter_by_stdev(self, stdev_cutoff=0.6, ratio_cutoff=4):
        """Filter ratios based on a standard deviation cutoff."""
        ratios = [p.ratio for p in self.peptides if p.ratio > 0]

        # no-op if less than two ratios
        if len(ratios) < 2:
            return

        stdev = self.special_stdev(ratios, default=None)
        mean = statistics.mean(ratios)
        minimum = min(ratios)

        # if stdev is tight enough relative to the mean
        # or if all the ratios are above a certain threshold
        # return unchanged
        # else, set everything but the minimum to zero
        if (stdev != None and stdev / mean < stdev_cutoff) or minimum > ratio_cutoff:
            return
        else:
            for p in self.peptides:
                if p.ratio != minimum:
                    p.ratio = 0

    def filter_by_mod_z(self, z_cutoff=3, ratio_cutoff=4):
        """Filter outliers by calculating modified z-score."""
        all_ratios = [p.ratio for p in self.peptides]
        non_zero_ratios = [r for r in all_ratios if r > 0]

        # no-op if less than two ratios or if all non-zero ratios are above cutoff
        if len(non_zero_ratios) < 3 or min(non_zero_ratios) > ratio_cutoff:
            return

        mod_z_scores = get_modified_z_scores(all_ratios)

        for z_score, p in zip(mod_z_scores, self.peptides):
            if z_score > z_cutoff:
                p.ratio = 0


    def filter_by_iqr(self, ratio_cutoff=4):
        """Filter outliers outside of IQR."""
        ratios = [p.ratio for p in self.peptides if p.ratio != 0]

        if len(ratios) < 2 or min(ratios) > ratio_cutoff:
            return

        lower, upper = get_iqr_bounds(ratios)

        for p in self.peptides:
            ratio = p.ratio
            if ratio != 0 and (ratio < lower or ratio > upper):
                p.ratio = 0


    def filter_20s(self, ratio_cutoff=4):
        """Filter out erroneous seeming 20s."""
        ratios = [p.ratio for p in self.peptides]

        if 20 not in ratios:
            return

        non_20s = [x for x in ratios if x != 20]
        non_zero_or_20 = [x for x in non_20s if x > 0]

        # 20s are stripped if the following conditions are met:
        # - the set of ratios is not just composed of 0s and 20s
        # - there is only one 20
        # - the lowest non-zero, non-20 value is below a cutoff
        if not set(ratios).issubset({0, 20}) and ratios.count(20) == 1:
            if not non_zero_or_20 or min(non_zero_or_20) < ratio_cutoff:
                for p in self.peptides:
                    if p.ratio == 20:
                        p.ratio = 0

    def filter_20s_by_stdev(self, stdev_cutoff=0.6, ratio_cutoff=4):
        """Filter out 20s from sets of ratios with high standard deviation."""
        ratios = [p.ratio for p in self.peptides if p.ratio != 0]
        non_20s = [x for x in ratios if x != 20]

        # require that there are at least two ratios
        # that there is a 20 in the set of ratios
        # and that the minimum ratio is below the ratio_cutoff
        # otherwise do nothing
        if len(ratios) < 2 or 20 not in ratios or min(ratios) > 4:
            return

        above_threshold = statistics.stdev(ratios) / statistics.mean(ratios) > stdev_cutoff

        if above_threshold:
            for p in self.peptides:
                if p.ratio == 20:
                    p.ratio = 0

    def filter_20s_by_ms2(self, min_ms2=2):
        """Only keep 20s with a certain number of ms2s."""
        self.peptides = [p for p in self.peptides if not(p.ratio == 20 and p.num_ms2 < min_ms2)]

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
    def special_stdev(ratios, default='-'):
        """Return the standard deviation of a set of ratios."""
        try:
            return statistics.stdev(ratios)
        except statistics.StatisticsError:
            return default

    @staticmethod
    def special_stderr(ratios, default='-'):
        """Return the standard error of a set of ratios."""
        try:
            return statistics.stdev(ratios) / math.sqrt(len(ratios))
        except:
            return default

    def __add__(self, protein):
        """Add one protein to another by concatenating peptides with the condition that they have the same uniprot."""
        if not type(self) == type(protein) or not protein.uniprot == self.uniprot:
            return

        new_protein = deepcopy(self)
        new_protein.peptides = self.peptides + protein.peptides

        return new_protein

    def __radd__(self, *args, **kwargs):
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
