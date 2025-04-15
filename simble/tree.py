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
import logging

logger = logging.getLogger(__package__)

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
        self.occupancy_time = None
        self.time_since_last_split = None
    
    def add_child(self, child):
        child.parent = self
        self.children.append(child)

    def write_newick(self, time_tree=False):
        name = f"cell_id={str(self.clone_id)}_{str(id(self.cell))}|location={self.cell.location.value}|generation={self.generation}"
        if self.occupancy_time is not None:
            name += f"|occupancy_time={self.occupancy_time}"
        if time_tree:
            branch_length = str(self.time_since_last_split) if self.time_since_last_split is not None else str(1)
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
        new.occupancy_time = self.occupancy_time
        new.time_since_last_split = self.time_since_last_split
        return new
    
    def get_occupancy_time(self):
        if self.occupancy_time is not None:
            return self.occupancy_time
        if self.parent is None:
            self.occupancy_time = 0
        elif self.cell.location is None or self.parent.cell.location is None:
            self.occupancy_time = 0
        elif self.cell.location == self.parent.cell.location:
            self.occupancy_time = self.parent.get_occupancy_time()+1
        else:
            self.occupancy_time = 1
        return self.occupancy_time

    def simplify_subtree(self):
        new_tree = _remove_non_splitting_nodes(self)
        return new_tree
    
    def prune_subtree(self, to_keep):
        new_tree = _build_tree_to_keep(self, to_keep)
        return new_tree
    
    def prune_up_tree(self):
        if len(self.children) > 0:
            # we don't want to prune this node or up if it still has children
            return
        if len(self.children) == 0:
            # if there are no children, we will prune this node and move up
            if self.parent is not None:
                self.parent.children.remove(self)
                parent = self.parent
                self.parent = None
                parent.prune_up_tree()           
            
        
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
            logger.warning(f"Pruning is still doing something")
            return None
    else:
        new_node = node.copy()
        for subtree in subtrees_to_keep:
            new_node.add_child(subtree)
        return new_node


def _remove_non_splitting_nodes(node, time_since_last_split=0, heavy_mutations_since_last_split=0, light_mutations_since_last_split=0):
    node.get_occupancy_time()
    new_node = node.copy()
    new_node.time_since_last_split = time_since_last_split+1
    new_node.heavy_mutations += heavy_mutations_since_last_split
    new_node.light_mutations += light_mutations_since_last_split
    if len(node.children) == 1:
        child = node.children[0]
        return _remove_non_splitting_nodes(child, new_node.time_since_last_split, new_node.heavy_mutations, new_node.light_mutations)
    else:
        for child in node.children:
            new_child = _remove_non_splitting_nodes(child)
            if new_child is not None:
                new_node.add_child(new_child)
        return new_node


