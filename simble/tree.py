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
        return _write_newick_iteratively(self, time_tree=time_tree)

    def _write_newick(self, time_tree=False, subtrees=None, recursion=False):
        name = f"{str(self.clone_id)}_{str(id(self.cell))}"
        labels= f"cell_id={str(self.clone_id)}_{str(id(self.cell))},location={self.cell.location.value},generation={self.generation}"
        if self.occupancy_time is not None:
            labels += f",occupancy_time={self.occupancy_time}"
            if self.time_since_last_split is not None:
                labels += f",occupancy={self.occupancy_time/self.time_since_last_split}"
        if time_tree:
            branch_length = str(self.time_since_last_split) if self.time_since_last_split is not None else str(1)
        else:
            branch_length = str(self.heavy_mutations+self.light_mutations)
        branch = f':{branch_length}'
        labels = f"[&{labels}]"
        if len(self.children)==0:
            children = ""
            # return name + labels + branch
        elif subtrees is None and recursion:
            children = "(" + ",".join([x.write_newick(time_tree=time_tree) for x in self.children]) + ")"
        elif subtrees is None or len(subtrees) == 0:
            children = ""
        else:
            children = "(" + ",".join(subtrees) + ")"
        return children + name + labels + branch
        
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
            return None
    else:
        new_node = node.copy()
        for subtree in subtrees_to_keep:
            new_node.add_child(subtree)
        return new_node


def _write_newick_iteratively(tree, time_tree=False):
    stack = [tree]
    children_newick = {}
    newick = ""

    def add_to_newick_dict(node, newick):
        # for memory efficiency, once we're adding this node's newick to its parent
        # we can remove it from the dict
        if node.parent is None:
            return newick
        if node.parent not in children_newick:
            children_newick[node.parent] = []
        children_newick[node.parent].append(newick)
        if node in children_newick:
            children_newick.pop(node)
        return ""
    
    while len(stack) > 0:
        current = stack.pop()
        number_of_children = len(current.children)
        child_newicks = children_newick.get(current, [])
        number_of_child_newicks = len(child_newicks)
        if len(current.children) == 0:
            # leaf node <- but this should be handled by number_of_children == number_of_child_newicks but i want to just test it first
            curr_newick = current._write_newick(time_tree=time_tree)
            newick += add_to_newick_dict(current, curr_newick)
        elif number_of_children == number_of_child_newicks:
            # all children have been processed
            curr_newick = current._write_newick(time_tree=time_tree, subtrees=child_newicks)
            # add the newick string to the parent and if there is no parent
            # i.e. we have the root, then we can just write the newick string
            newick += add_to_newick_dict(current, curr_newick)
        else:
            # not all children have been processed
            # this node can't be processed yet, so push it back onto the stack
            stack.append(current)
            # then push all children onto the stack so the children are above current node
            for child in current.children:
                stack.append(child)
    return newick


def simplify_tree(root):
    root.get_occupancy_time()
    new_root = root.copy()
    subtrees = [(new_root, child, 0, 0, 0) for child in root.children]
    while len(subtrees) > 0:
        parent, current_node, time_since_last_split, heavy_mutations_since_last_split, light_mutations_since_last_split = subtrees.pop(0)
        current_node.get_occupancy_time()
        if len(current_node.children) == 1:
            # we're removing this node
            child = current_node.children[0]
            heavy_mutations_since_last_split += current_node.heavy_mutations
            light_mutations_since_last_split += current_node.light_mutations
            subtrees.append((parent, child, time_since_last_split+1, heavy_mutations_since_last_split, light_mutations_since_last_split))
        else:
            # we're keeping this node
            new_node = current_node.copy()
            new_node.time_since_last_split = time_since_last_split+1
            new_node.heavy_mutations += heavy_mutations_since_last_split
            new_node.light_mutations += light_mutations_since_last_split
            parent.add_child(new_node)
            for child in current_node.children:
                subtrees.append((new_node, child, 0, 0, 0))
    return new_root





