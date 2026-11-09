"""
 Copyright (C) 2024 Jessie Fielding

 This file is part of simble.

 simble is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 simble is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with simble.  If not, see <https://www.gnu.org/licenses/>.
 """

import json
import os
import re
import tempfile
import unittest

import numpy as np
import pandas as pd

from simble.cell import Cell
from simble.chain import HeavyChain, LightChain
from simble.location import LocationName
from simble.settings import Settings, s
from simble.simble import do_simulation, process_results
from simble.tree import Node, get_internal_nodes

ENTROPY = 54897022524486695084299880814690718190


def run_clone(internal_nodes, uniform=False, clone_id=1):
    """Runs one short clone with a fixed seed and returns its results."""
    settings = Settings()
    settings.LOCATIONS[0].sample_times = list(range(0, 11, 5))
    settings.LOCATIONS[1].sample_times = list(range(0, 11, 5))
    settings.LOCATIONS[0].migration_rate = 0.5
    settings.INTERNAL_NODES = internal_nodes
    if uniform:
        settings.UNIFORM = True
    seed = np.random.SeedSequence(ENTROPY)
    with tempfile.NamedTemporaryFile(mode="w") as tmpf:
        json.dump(settings, tmpf, default=lambda o: o.encode(), indent=4)
        tmpf.flush()
        return do_simulation(clone_id, seed, None, tmpf.name)


def parse_newick(newick):
    """Parses a simble Newick string without recursion.

    Returns:
        dict: node name -> {"parent": name or None, "length": float,
            "children": [names]}
    """
    newick = newick.strip().rstrip(";")
    # remove the [&...] annotations; names and branch lengths stay
    newick = re.sub(r"\[&[^\]]*\]", "", newick)
    nodes = {}
    stack = [[]]  # each item is the list of child names at that depth
    i = 0
    while i < len(newick):
        char = newick[i]
        if char == "(":
            stack.append([])
            i += 1
        elif char == ",":
            i += 1
        elif char == ")":
            i += 1
        else:
            match = re.match(r"([^(),:]+):([^(),:]+)", newick[i:])
            name, length = match.group(1), float(match.group(2))
            i += match.end()
            children = []
            if newick[i - match.end() - 1] == ")":
                children = stack.pop()
            nodes[name] = {"parent": None, "length": length, "children": children}
            for child in children:
                nodes[child]["parent"] = name
            stack[-1].append(name)
    return nodes


def tips_below(nodes, name):
    """Returns the sorted names of the tips below a node."""
    tips = []
    stack = [name]
    while stack:
        current = stack.pop()
        if len(nodes[current]["children"]) == 0:
            tips.append(current)
        stack.extend(nodes[current]["children"])
    return sorted(tips)


def fasta_records(fasta):
    """Returns a dict of FASTA record ID -> sequence."""
    lines = [x for x in fasta.split("\n") if x != ""]
    return {
        lines[i][1:].split("|")[0]: lines[i + 1]
        for i in range(0, len(lines), 2)
        }


def fasta_headers(fasta):
    """Returns a dict of FASTA record ID -> dict of header fields."""
    headers = {}
    for line in fasta.split("\n"):
        if line.startswith(">"):
            name, *fields = line[1:].split("|")
            headers[name] = dict(x.split("=", 1) for x in fields)
    return headers


def hamming(a, b):
    """Returns the number of positions where two same-length strings differ."""
    assert len(a) == len(b)
    return sum(x != y for x, y in zip(a, b))


class TestInternalNodes(unittest.TestCase):
    """Tests for the internal node (true ancestral sequence) output."""

    @classmethod
    def setUpClass(cls):
        cls.off = run_clone(internal_nodes=False)
        cls.on = run_clone(internal_nodes=True)
        cls.internal = cls.on["internal_airr"]
        cls.tree = parse_newick(cls.on["simplified_tree"])
        cls.time_tree = parse_newick(cls.on["simplified_time_tree"])

    def test_off_by_default(self):
        """The setting is off by default, and gives no internal output."""
        self.assertFalse(Settings().INTERNAL_NODES)
        self.assertIsNone(self.off["internal_airr"])
        self.assertIsNone(self.off["internal_fasta"])

    def test_tip_output_unchanged(self):
        """The flag does not change the sampled-cell output or the RNG stream."""
        drop = ["sequence_id", "cell_id"]
        pd.testing.assert_frame_equal(
            self.off["airr"].drop(columns=drop),
            self.on["airr"].drop(columns=drop))
        pd.testing.assert_frame_equal(self.off["pop_data"], self.on["pop_data"])

    def test_ids_match_both_simplified_trees(self):
        """Each internal node in the Newick trees has rows, and the reverse."""
        self.assertEqual(set(self.tree), set(self.time_tree))
        root = [x for x, y in self.tree.items() if y["parent"] is None][0]
        tree_internal = {
            x for x, y in self.tree.items()
            if len(y["children"]) > 0 or x == root
            }
        self.assertGreater(len(tree_internal), 1)
        self.assertEqual(tree_internal, set(self.internal["cell_id"]))
        root_rows = self.internal[self.internal["node_type"] == "root"]
        self.assertEqual(set(root_rows["cell_id"]), {root})

    def test_two_chains_per_node(self):
        """Paired cells have one heavy row and one light row for each node."""
        counts = self.internal.groupby("cell_id").size()
        self.assertTrue((counts == 2).all())

    def test_root_is_naive(self):
        """The root sequence is the naive (germline) sequence."""
        root_rows = self.internal[self.internal["node_type"] == "root"]
        for _, row in root_rows.iterrows():
            self.assertEqual(row["sequence_alignment"], row["germline_alignment"])
            self.assertEqual(row["branch_mutations"], 0)
            self.assertTrue(pd.isna(row["parent_id"]))

    def test_tree_columns_match_newick(self):
        """parent_id, branch lengths and descendant tips agree with the trees."""
        for cell_id, rows in self.internal.groupby("cell_id"):
            node = self.tree[cell_id]
            row = rows.iloc[0]
            if node["parent"] is None:
                self.assertTrue(pd.isna(row["parent_id"]))
            else:
                self.assertEqual(row["parent_id"], node["parent"])
            self.assertEqual(rows["branch_mutations"].sum(), node["length"])
            self.assertEqual(row["branch_generations"], self.time_tree[cell_id]["length"])
            tips = tips_below(self.tree, cell_id)
            self.assertEqual(row["descendant_tips"].split(";"), tips)
            self.assertEqual(row["n_descendant_tips"], len(tips))

    def test_root_tips_are_sampled_cells(self):
        """The tips below the root are all the sampled cells."""
        root_row = self.internal[self.internal["node_type"] == "root"].iloc[0]
        sampled = sorted(set(self.on["airr"]["cell_id"]))
        self.assertEqual(root_row["descendant_tips"].split(";"), sampled)

    def test_sequences_agree_with_branch_mutations(self):
        """Parent-to-child sequence differences are not more than the mutations
        on that branch. A site can mutate more than once, so they can be fewer."""
        sequences = pd.concat([self.internal, self.on["airr"]])
        sequences = sequences.set_index(["cell_id", "locus"])["sequence_alignment"]
        total_mutations = 0
        for child, node in self.tree.items():
            if node["parent"] is None:
                continue
            differences = sum(
                hamming(sequences[(child, locus)], sequences[(node["parent"], locus)])
                for locus in sequences[child].index
                )
            self.assertLessEqual(differences, node["length"])
            total_mutations += node["length"]
        self.assertGreater(total_mutations, 0)

    def test_fasta(self):
        """The internal node FASTA agrees with the internal node AIRR output."""
        self.assert_fasta_matches_airr(self.on["internal_fasta"], self.internal)
        for header in fasta_headers(self.on["internal_fasta"]).values():
            self.assertNotIn("sample_time", header)
            self.assertIn("node_type", header)

    def test_sampled_fasta(self):
        """The sampled-cell FASTA agrees with the sampled-cell AIRR output, and
        its IDs have the clone prefix, as in the AIRR and Newick output."""
        self.assert_fasta_matches_airr(self.on["fasta"], self.on["airr"])
        for name in fasta_headers(self.on["fasta"]):
            self.assertTrue(name.startswith(f"{self.on['clone_id']}_"))
        tips = {x for x, y in self.tree.items() if len(y["children"]) == 0}
        self.assertEqual(tips, set(self.on["airr"]["cell_id"]))

    def assert_fasta_matches_airr(self, fasta, airr):
        """Checks the FASTA IDs, sequences and loci against the AIRR rows."""
        records = fasta_records(fasta)
        headers = fasta_headers(fasta)
        self.assertEqual(set(records), set(airr["sequence_id"]))
        self.assertEqual(len(records), len(airr))
        for _, row in airr.iterrows():
            self.assertEqual(records[row["sequence_id"]], row["sequence"])
            self.assertEqual(headers[row["sequence_id"]]["locus"], row["locus"])

    def test_process_results_writes_files(self):
        """process_results writes the internal node files only when asked."""
        with tempfile.TemporaryDirectory() as tmpdir:
            old = (s.RESULTS_DIR, s.INTERNAL_NODES, s.FASTA, s.DEV)
            try:
                s.RESULTS_DIR, s.DEV, s.FASTA = tmpdir, False, True
                s.INTERNAL_NODES = False
                process_results([self.off])
                self.assertFalse(os.path.exists(f"{tmpdir}/all_internal_nodes_airr.tsv"))
                self.assertFalse(os.path.exists(f"{tmpdir}/all_internal_nodes.fasta"))
                s.INTERNAL_NODES = True
                process_results([self.on])
                written = pd.read_csv(f"{tmpdir}/all_internal_nodes_airr.tsv", sep="\t")
                self.assertEqual(len(written), len(self.internal))
                self.assertTrue(os.path.exists(f"{tmpdir}/all_internal_nodes.fasta"))
            finally:
                s.RESULTS_DIR, s.INTERNAL_NODES, s.FASTA, s.DEV = old


class TestInternalNodesUniform(unittest.TestCase):
    """Internal node output for single chain (uniform) simulations."""

    def test_one_chain_per_node(self):
        result = run_clone(internal_nodes=True, uniform=True)
        internal = result["internal_airr"]
        self.assertGreater(len(internal), 0)
        self.assertTrue((internal.groupby("cell_id").size() == 1).all())
        self.assertEqual(set(internal["locus"]), {"IGH"})
        for fasta in [result["fasta"], result["internal_fasta"]]:
            headers = fasta_headers(fasta)
            self.assertGreater(len(headers), 0)
            self.assertTrue(all(x.endswith("_heavy") for x in headers))
            self.assertEqual({x["locus"] for x in headers.values()}, {"IGH"})


class TestFastaLocus(unittest.TestCase):
    """The FASTA locus comes from the naive data (for example, kappa chains)."""

    def make_cell(self, light_locus):
        sequence = "CAGGTGCAGCTGGTGCAGTCTGGG"
        heavy = HeavyChain(sequence, gapped_seq=sequence)
        heavy.airr_constants = {"locus": "IGH"}
        light = LightChain(sequence, gapped_seq=sequence)
        light.airr_constants = {"locus": light_locus} if light_locus else {}
        return Cell(heavy, light, created_at=0)

    def test_kappa(self):
        headers = fasta_headers(self.make_cell("IGK").as_fasta(3, name="1_42"))
        self.assertEqual(headers["1_42_heavy"]["locus"], "IGH")
        self.assertEqual(headers["1_42_light"]["locus"], "IGK")
        self.assertEqual(headers["1_42_light"]["sample_time"], "3")

    def test_default_locus(self):
        """Without a locus in the naive data, a light chain is IGL."""
        cell = self.make_cell(None)
        self.assertEqual(cell.light_chain.locus, "IGL")

    def test_no_name_keeps_old_ids(self):
        """Without a name (library use), the ID is the cell id, as before."""
        cell = self.make_cell("IGK")
        headers = fasta_headers(cell.as_fasta(3))
        self.assertIn(f"{id(cell)}_heavy", headers)


class FakeCell: # pylint: disable=too-few-public-methods
    """A minimal cell for tree unit tests."""
    class _Location: # pylint: disable=too-few-public-methods
        name = LocationName.OTHER
    location = _Location()


class TestGetInternalNodes(unittest.TestCase):
    """Unit tests for tree.get_internal_nodes."""

    def test_small_tree(self):
        root = Node(FakeCell(), clone_id=1)
        a = Node(FakeCell(), generation=1, clone_id=1)
        b = Node(FakeCell(), generation=2, clone_id=1)
        c = Node(FakeCell(), generation=2, clone_id=1)
        d = Node(FakeCell(), generation=1, clone_id=1)
        root.add_child(a)
        root.add_child(d)
        a.add_child(b)
        a.add_child(c)
        result = get_internal_nodes(root)
        self.assertEqual([x for x, _ in result], [root, a])
        self.assertEqual(result[0][1], sorted([b.name, c.name, d.name]))
        self.assertEqual(result[1][1], sorted([b.name, c.name]))

    def test_root_only(self):
        root = Node(FakeCell(), clone_id=1)
        self.assertEqual(get_internal_nodes(root), [(root, [root.name])])

    def test_deep_tree_no_recursion_error(self):
        """A tree deeper than the default recursion limit does not fail."""
        root = Node(FakeCell(), clone_id=1)
        current = root
        for i in range(5000):
            child = Node(FakeCell(), generation=i + 1, clone_id=1)
            current.add_child(child)
            current = child
        result = get_internal_nodes(root)
        self.assertEqual(len(result), 5000)
        self.assertEqual(result[0][1], [current.name])


if __name__ == '__main__':
    unittest.main()
