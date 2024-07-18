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

class Node:
    def __init__(self, cell, parent=None, heavy_mutations=0, light_mutations=0, generation=0, clone_id=None):
        self.cell=cell
        self.parent=parent
        self.heavy_mutations=heavy_mutations
        self.light_mutations=light_mutations
        self.generation=0
        self.children=[]
        self.antigen=0
        self.generation=generation
        self.clone_id=clone_id if clone_id else parent.clone_id if parent else -1
        self.sampled_time = None
    
    def add_child(self, child):
        child.parent = self
        self.children.append(child)

    def write_newick(self, time_tree=False):
        name = f"{str(self.clone_id)}_{str(id(self.cell))}_{self.cell.location.value}_{self.generation}"
        if time_tree:
            branch_length = str(1)
        else:
            branch_length = str(self.heavy_mutations+self.light_mutations)
        branch = f':{branch_length}'
        if len(self.children)==0:
            return name + branch
        else:
            children = "(" + ",".join([x.write_newick(time_tree=time_tree) for x in self.children]) + ")"
            return children + name + branch
        
    def copy(self):
        new = Node(
            self.cell,
            None,
            self.heavy_mutations,
            self.light_mutations,
            self.generation,
            clone_id=self.clone_id
        )
        new.antigen = self.antigen
        return new
    
    def prune_subtree(self, to_keep):
        new_tree = _build_tree_to_keep(self, to_keep)
        return new_tree
        
def _build_tree_to_keep(node, to_keep):
    subtrees_to_keep = []
    for child in node.children:
        subtree = _build_tree_to_keep(child, to_keep)
        if subtree is not None:
            subtrees_to_keep.append(subtree)
    if len(subtrees_to_keep) == 0:
        if id(node.cell) in to_keep:
            return node.copy()
        else:
            return None
    else:
        new_node = node.copy()
        for subtree in subtrees_to_keep:
            new_node.add_child(subtree)
        return new_node
