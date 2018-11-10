"""Calculate the possible reactions from the parents and daughters in a
nuclear reaction using various models.
"""
# pylint: disable=too-many-instance-attributes, no-self-use, no-member
# pylint: disable=too-few-public-methods, unsubscriptable-object
import gzip
import hashlib
import itertools
import json
import logging
import operator
import os
from os.path import expanduser
import pickle

from .nubase import Nuclides, Electron
from .calculations import (
    IsotopicDecay,
    GamowSuppressionFactor,
    GeigerNuttal,
    Gamow2,
    ReactionEnergy,
)


LENRMC_DIR = os.path.join(expanduser('~'), '.reactions')


class RejectCombination(RuntimeError):
    """Sentinel exception used to flag a specific combination of parent and
    daughter nuclides as being impossible.
    """
    pass


class GammaPhoton:
    """Represent a gamma photon that results during a nuclear reaction."""

    def __init__(self):
        self.mass_number = 0
        self.is_baryon = False
        self.full_label = self.label = 'ɣ'
        self.is_stable = False
        self.spin_and_parity = '1-'
        self.numbers = (0, 0)
        self.notes = {'ɣ'}


class Reaction:
    """Model the parents and possible daughters of a nuclear reaction."""

    @classmethod
    def load(cls, **kwargs):
        """Factory method returning a set of reactions for a given set of
        target parents and target daughters.
        """
        nuclides = Nuclides.data()
        reactants = ((num, nuclides[s]) for num, s in kwargs['reactants'])
        daughters = ((num, nuclides[s]) for num, s in kwargs['daughters'])
        return cls(reactants, daughters, **kwargs)

    _noteworthy = {
        '4He': 'α',
        'n':   'n',
        't':   't',
    }

    def __init__(self, lvalues, rvalues, **kwargs):
        self.model_name = kwargs.get('model')
        self.initial_lvalues = list(lvalues)
        self.rvalues = list(rvalues)
        self.q_value = ReactionEnergy(self).value
        self.is_stable = self._is_stable()
        self.any_excited = self._any_excited()
        self.lvalue_delim = self.rvalue_delim = '+'
        # Photons are not counted as a daughter.
        self.daughter_count = sum(n for (n, d) in self.rvalues)
        if self.is_single_body and len(self.initial_lvalues) > 1:
            self.rvalues.append((1, GammaPhoton()))
        gamow = self.gamow(**kwargs)
        self.gamow_value = gamow.value() if gamow else -1

    @property
    def lvalues(self):
        """The parents of the reaction."""
        if self.model_name == 'stimulated-decay':
            return [(1, Electron())] + self.initial_lvalues
        return self.initial_lvalues

    @property
    def notes(self):
        """Notes that are relevant to this reaction to be included in the
        output.
        """
        notes = set()
        for _, daughter in self.rvalues:
            note = self._noteworthy.get(daughter.label)
            if note:
                notes.add(note)
            for _, parent in self.initial_lvalues:
                if self._neutron_transfer(daughter, parent):
                    notes.add('n-transfer')
        if self.is_stable:
            notes.add('in nature')
        for _, daughter in self.rvalues:
            notes |= daughter.notes
        return notes

    def _neutron_transfer(self, daughter, parent):
        result = map(operator.add, parent.numbers, (1, 0))
        return daughter.numbers == tuple(result)

    def _is_stable(self):
        return all(d.is_stable for num, d in self.rvalues)

    def _any_excited(self):
        combined = self.rvalues + self.initial_lvalues
        return any(n.is_excited for num, n in combined)

    def _decay_components(self):
        values = [p for num, p in self.rvalues if p.is_baryon]
        if len(values) != 2:
            return None
        nuclide0, nuclide1 = values
        smaller = min(values, key=lambda v: v.mass_number)
        larger = nuclide0 if nuclide1.mass_number == smaller.mass_number else nuclide1
        _, parent = max(self.initial_lvalues, key=lambda t: t[1].atomic_number)
        return parent, (smaller, larger)

    @property
    def has_electron_parent(self):
        """Does this reaction involve an electon?"""
        return any(isinstance(p, Electron) for num, p in self.initial_lvalues)

    @property
    def is_single_body(self):
        """Is this reaction a decay?"""
        if len(self.rvalues) > 1:
            return False
        return all(num == 1 for num, c in self.rvalues)

    def geiger_nuttal(self):
        """Do the Geiger-Nuttal computation for the decay components of this
        decay.
        """
        return GeigerNuttal.load(self._decay_components(), self.q_value)

    def gamow(self, **kwargs):
        """Compute the gamow suppression factor and related details for a given
        reaction.
        """
        return GamowSuppressionFactor.load(
            self._decay_components(),
            self.q_value,
            **kwargs,
        )

    def gamow2(self):
        """Compute the gamow suppression factor and related computations for
        a given reaction using a different method from the other calculation
        with the same name.
        """
        return Gamow2.load(self._decay_components(), self.q_value)

    def decay(self, **kwargs):
        """Compute possible decay products for a given pair of parent nuclides."""
        return IsotopicDecay.load(
            self._decay_components(),
            self.q_value,
            **kwargs
        )


def vectors3(integer):
    """Produce 3-tuples that sum up to the input."""
    for i in range(integer):
        j = integer - i
        for k in range(j):
            yield (j - k, k, i)


class RegularCombinations:
    """Compute the possible combinations of neutrons and protons in the
    daughters for a given set of parent nuclides, using a cached result
    if one is available.
    """
    basedir = os.path.expanduser('~/.reactions/objects')

    def __init__(self, totals):
        self.totals = totals
        self.cache_key = self._cache_key(totals)

    def iterator(self):
        """Return an iterator holding the possible combinations of protons
        and neutrons for a given pair of parent nuclides.  These are converted
        to daughter nuclides at a later step.
        """
        iterator = self._cached_results()
        if iterator:
            yield from iterator
            return

        mass_number, atomic_number = self.totals
        results, seen = [], set()
        for masses in vectors3(mass_number):
            for protons in vectors3(atomic_number):
                daughters = []

                try:
                    for i, mass in enumerate(masses):
                        proton_count = protons[i]
                        if mass < proton_count:
                            raise RejectCombination
                        pair = (mass, proton_count)
                        if (0, 0) == pair:
                            continue
                        daughters.append(pair)
                except RejectCombination:
                    continue

                daughters = tuple(sorted(daughters))
                if daughters in seen:
                    continue

                seen.add(daughters)
                results.append(daughters)
                yield daughters
        self._cache_results(results)

    def _cache_key(self, totals):
        string = json.dumps(totals, sort_keys=True).encode('utf-8')
        return hashlib.sha1(string).hexdigest()

    def _cache_results(self, results):
        try:
            os.makedirs(self.basedir)
        except FileExistsError:
            pass
        with gzip.open(self.cache_path, 'wb+') as file:
            pickle.dump(results, file)

    @property
    def cache_path(self):
        """File path for the cached results for this reaction."""
        return os.path.join(self.basedir, self.cache_key)

    def _cached_results(self):
        if not os.path.exists(self.cache_path):
            return None
        logging.info('reading previously computed values from cache')
        with gzip.open(self.cache_path, 'rb') as file:
            combinations = pickle.loads(file.read())
        return combinations


def regular_combinations(totals):
    """Convenience function for invoking the RegularCombinations class."""
    return RegularCombinations(totals).iterator()


def add_numbers(*numbers):
    """Add the mass and atomic numbers for different pairs of input."""
    return tuple(map(operator.add, *numbers))


def normalize(pair):
    """Convert the pair into a normalized list of atomic and mass numbers."""
    o_mass, o_atomic = pair
    if o_mass < 0:
        return [(o_mass, o_atomic)]
    if o_atomic == 0:
        return [(1, 0) for _ in range(o_mass)]
    if o_mass == o_atomic:
        return [(1, 1) for _ in range(o_mass)]
    return [(o_mass, o_atomic)]


class Model:
    """Model is the base class for several different ways of calculating
    a set of possible nuclear reactions from a given set of inputs.
    """
    def sort_key(self, reaction):
        """Sort reactions according to the energy released."""
        kev = reaction.q_value.kev
        return kev > 0, kev

    def parents(self, parents, _):
        """Default no-op implementation of the parents method. In some
        subclasses, the parents are modified here."""
        yield parents

    def _smaller_and_larger(self, reactants):
        "This helper method assumes a two-body reaction."
        num0, smaller = min(reactants, key=lambda t: t[1].mass_number)
        num1, larger = max(reactants, key=lambda t: t[1].mass_number)
        assert 1 == num0 == num1
        return [(num0, smaller), (num1, larger)]


class StandardModel(Model):
    """StandardModel captures a regular nuclear reaction without any bells
    or whistles.
    """
    def __call__(self, reactants):
        numbers = [num * n.numbers for num, n in reactants]
        mass_number, atomic_number = add_numbers(*numbers)
        return regular_combinations((mass_number, atomic_number))


class PionExchangeModel(Model):
    """PionExchangeModel is an abstract base class for two models in whic
    neutrons are exchanged with protons and vice versa.
    """

    # Careful!
    # 7Li -> 7Be
    # 51Ni -> 52Co
    _transformations = [
        [(0, 1), (0, -1)],
        [(0, -1), (0, 1)],
    ]

    def __call__(self, reactants):
        (_, smaller), (_, larger) = self._smaller_and_larger(reactants)
        seen = set()
        for adj_left, adj_right in self._transformations:
            p_left = add_numbers(smaller.numbers, adj_left)
            p_right = add_numbers(larger.numbers, adj_right)
            yield from self._combinations(seen, p_left, p_right)


class PionExchangeAndDecayModel(PionExchangeModel):
    """PionExchangeAndDecayModel is a speculative model in which neutrons are
    allowed to change to protons and vice versa, along with other decays.
    """
    def _combinations(self, seen, p_left, p_right):
        outcomes = regular_combinations(p_left)
        for pairs in outcomes:
            daughters = normalize(p_right)
            for pair in pairs:
                daughters.extend(normalize(pair))
            daughters = tuple(sorted(daughters))
            if daughters in seen:
                continue
            seen.add(daughters)
            yield daughters


class StrictPionExchangeModel(PionExchangeModel):
    """StrictPionExchangeModel is a speculative model in which neutrons are
    allowed to flip to protons and vice versa through the exchange of pions.
    """
    def _combinations(self, seen, p_left, p_right):
        daughters = tuple(sorted(normalize(p_left) + normalize(p_right)))
        if daughters not in seen:
            yield daughters
        seen.add(daughters)


class ElectronMediatedDecayModel(Model):
    """ElectronMediatedDecayModel is a speculative model in which weak-force
    mediated reactions are allowed to be sped up through the presence of
    election screening.
    """
    _transformations = [
        # electron capture
        [(0, -1), [(0, 0), (0, -1)]],
        # β- decay
        [(0, 1), [(0, 0), (0, -1)]],
        # β-β- decay
        [(0, 2), [(0, 0), (0, 0), (0, -1), (0, -1)]],
        # ejection of a proton
        [(-1, -1), [(1, 1)]],
        # ejection of an alpha particle
        [(-4, -2), [(4, 2), (0, -1)]],
    ]

    def parents(self, parents, _):
        yield list(parents) + [(1, Electron())]

    def __call__(self, reactants):
        assert len(reactants) == 1
        num, nuclide = reactants[0]
        assert num == 1
        for transformation, rvalues2 in self._transformations:
            rvalue = add_numbers(nuclide.numbers, transformation)
            yield tuple(normalize(rvalue) + rvalues2)


class SeparatedNuclideModel(Model):
    """Find those parents who, if separated from a specific
    nuclide, would fission in an exothermic reaction, if fission
    were feasible.
    """
    nuclides = Nuclides.data()

    def parents(self, parents, daughters):
        assert len(parents) == 1
        _, daughter = parents[0]
        assert len(daughters) == 2
        parent_numbers = add_numbers(daughters[0].numbers, daughters[1].numbers)
        for parent in self.nuclides.isomers[parent_numbers]:
            if not parent.in_nature:
                continue
            if parent.mass_number == daughter.mass_number:
                continue
            yield [(1, parent)] + [(1, Electron())]

    def __call__(self, reactants):
        assert len(reactants) == 1
        num, nuclide0 = reactants[0]
        assert num == 1
        for nuclide1 in iter(Nuclides.data()):
            yield tuple([nuclide0.numbers, nuclide1.numbers])


class ElectronMediatedFissionModel(Model):
    """ElectronMediatedFissionModel is a speculative model in which electron
    screening is allowed to induce fission in otherwise relatively stable
    isotopes.
    """

    def sort_key(self, reaction):
        gamow = reaction.gamow_value
        kev = reaction.q_value.kev
        return -gamow, kev > 0, kev

    def __call__(self, reactants):
        assert len(reactants) == 1
        num, nuclide0 = reactants[0]
        assert num == 1
        return regular_combinations(nuclide0.numbers)


MODELS = {
    'standard':              StandardModel(),
    'pion-exchange':         PionExchangeAndDecayModel(),
    'strict-pion-exchange':  StrictPionExchangeModel(),
    'induced-decay':         ElectronMediatedDecayModel(),
    'separated-nuclide':     SeparatedNuclideModel(),
    'induced-fission':       ElectronMediatedFissionModel(),
}


class Combinations:
    """For a given set of input isotopes, there will be a number of different
    daughters that could be produced with the same number of neutrons and
    protons. The `Combinations` class encapsulates the calculation of various
    sets of daughters under the preservation of the neutron and proton counts.
    """

    @classmethod
    def load(cls, **kwargs):
        """Compute a set of nuclide combinations from a given set of reactants."""
        parents = kwargs['reactants']
        del kwargs['reactants']
        return cls(parents, **kwargs)

    def __init__(self, parents, **kwargs):
        self.model_name = kwargs.get('model') or 'standard'
        self._model = MODELS[self.model_name]
        self._parents = list(parents)
        self._kwargs = kwargs
        self._lower_bound = float(kwargs.get('lower_bound', 0))
        self._upper_bound = float(kwargs.get('upper_bound', 500000))
        self._excited = kwargs.get('excited')
        self.daughter_count = {int(c) for c in kwargs.get('daughter_count', '').split(',') if c}

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self._parents)

    def sort_key(self, reactions):
        """Sort our reaction output."""
        return self._model.sort_key(reactions)

    def _reactions(self):
        nuclides = Nuclides.data()
        for daughters in self._daughters():
            isomers = [nuclides.isomers[pair] for pair in daughters]
            if not all(isomers):
                continue
            yield from itertools.product(*isomers)

    def _daughters(self):
        return self._model(self._parents)

    def reactions(self):
        """Return an iterator that converts a set of nuclide combinations
        into a set of reactions.
        """
        if 'daughters' in self._kwargs:
            reactants = [(num, (n.label, '0')) for num, n in self._parents]
            reaction = Reaction.load(reactants=reactants, **self._kwargs)
            if self._allowed(reaction):
                yield reaction
            return

        for daughters in self._reactions():
            all_parents = self._model.parents(self._parents, daughters)
            for parents in all_parents:
                rvalues = ((1, d) for d in daughters)
                reaction = Reaction(parents, rvalues, **self._kwargs)
                if not self._allowed(reaction):
                    continue
                yield reaction

    def _allowed(self, reaction):
        conditions = [
            reaction.q_value.kev > self._lower_bound,
            reaction.q_value.kev <= self._upper_bound,
        ]
        if self.daughter_count:
            conditions.append(
                any(c == reaction.daughter_count for c in self.daughter_count)
            )
        if not self._excited:
            conditions.append(not reaction.any_excited)
        return all(conditions)
