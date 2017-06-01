"""Defines peptide type."""


class Peptide():
    """Container for a peptide."""

    delchars = {ord(c): None for c in map(chr, range(256)) if not c.isalpha()}

    def __init__(self, sequence, mass, charge, segment, ratio, rsquared):
        """Init peptide."""
        self.sequence = sequence
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
        return 'Peptide(sequence={}, mass={}, charge={}, segment={}, ratio={}, rsquared={})'.format(
            self.sequence,
            self.mass,
            self.charge,
            self.segment,
            self.ratio,
            self.rsquared
        )
