# flake8: noqa
from cimagex.dataset import Dataset
from ciamgex.protein import Protein
from cimagex.peptide import Peptide
from cimagex import make_dataset
from copy import deepcopy
import json
import statistics
import unittest


class TestCase(unittest.TestCase):
    def setUp(self):
        with open('test.json') as f:
            self.data = json.load(f)

    def test_rsquared_filter(self):
        dataset = self._make_dataset(self.data['test_rsquared'])

        # test filter at protein level
        self._test_protein_rsquared_cutoff(dataset.proteins[0], 0.8, 0)
        self._test_protein_rsquared_cutoff(dataset.proteins[0], 0.5, 1)

        # test filter at dataset level
        self._test_dataset_rsquared_cutoff(dataset, 0.8, 0)
        self._test_dataset_rsquared_cutoff(dataset, 0.5, 1)

    def _test_protein_rsquared_cutoff(self, protein, cutoff, expect):
        protein = deepcopy(protein)
        protein.apply_rsquared_cutoff(cutoff)
        self.assertEqual(len(protein.peptides), expect)

    def _test_dataset_rsquared_cutoff(self, dataset, cutoff, expect):
        dataset = deepcopy(dataset)
        dataset.apply_rsquared_cutoff(cutoff)
        self.assertEqual(len(dataset.proteins), expect)

    def test_whitelist_filter(self):
        dataset = self._make_dataset(self.data['test_whitelist'])
        dataset.apply_whitelist_filter(['Albatross'])
        self.assertEqual(len(dataset.proteins), 1)
        self.assertEqual(dataset.get_unique_uniprots(), {'Albatross'})

    def test_uniqueness_filter(self):
        dataset = self._make_dataset(self.data['test_unique'])
        self.assertEqual(len(dataset.proteins), 2)

        dataset.apply_unique_filter(1)
        self.assertEqual(len(dataset.proteins), 2)

        dataset.apply_unique_filter(2)
        self.assertEqual(len(dataset.proteins), 1)

        dataset.apply_unique_filter(3)
        self.assertEqual(len(dataset.proteins), 0)

    def test_stats(self):
        dataset = self._make_dataset(self.data['test_stats'])
        dataset.generate_stats(inverse=True)
        self.assertEqual(dataset.proteins[0].mean, 0.4)
        self.assertEqual(dataset.proteins[0].stdev, statistics.stdev([0.2, 0.5, 0.5]))
        self.assertEqual(dataset.proteins[0].median, 0.5)

    def test_add_protein_to_protein(self):
        dataset1 = self._make_dataset(self.data['test_addition_1'])
        dataset2 = self._make_dataset(self.data['test_addition_2'])
        protein1 = dataset1.get('Albatross')
        protein2 = dataset2.get('Albatross')
        combined_protein = protein1 + protein2
        self.assertEqual(len(combined_protein.peptides), 2)

    def test_add_protein_to_dataset(self):
        dataset = self._make_dataset(self.data['test_addition_1'])
        dataset2 = self._make_dataset(self.data['test_addition_2'])

        # check if protein (with extra unique peptide) that already exists gets added correctly
        dataset.add(dataset2.get('Albatross'))
        self.assertEqual(len(dataset.get('Albatross').peptides), 2)

        # check if protein that doesn't already exist gets added correctly
        dataset.add(dataset2.get('Chicken'))
        self.assertEqual(len(dataset.proteins), 3)

    def test_dataset_addition(self):
        dataset1 = self._make_dataset(self.data['test_addition_1'])
        dataset2 = self._make_dataset(self.data['test_addition_2'])
        combined_datasets = dataset1 + dataset2

        # check to see if proteins are getting combined as a superset
        self.assertEqual(len(combined_datasets.proteins), 3)

        # check to see if unique peptides are added together
        protein = combined_datasets.get('Albatross')
        self.assertEqual(len(protein.peptides), 2)

        # check to see if stats are calculated correctly
        combined_datasets.generate_stats()
        self.assertEqual(protein.mean, 0.55)


    def _make_dataset(self, data):
        return make_dataset(data['headers'], data['data'])



if __name__ == '__main__':
    unittest.main()
