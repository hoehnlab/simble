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
# pylint: disable=expression-not-assigned

import logging
import pandas as pd
import numpy as np
import math

from .settings import s, Encodable
from .location import LocationName

# prob remove these
import matplotlib.pyplot as plt
import matplotlib.patches as plt_ptch
from matplotlib.animation import FuncAnimation

# fig, ax = plt.subplots()

logger = logging.getLogger(__package__)

class SpatialPlane(Encodable):
    """Represents a spatial plane in the simulation.

    Attributes:
    """
    def __init__(
            self,
            boundary,
            rate,
            polygons=None
    ):
        """Initializes a SpatialPlane instance.
        
        Args:
            name (str): The name of the spatial plane.
            settings (SpatialPlaneSettings): The settings for the spatial plane.
     """
        self.boundary = boundary
        self.polygons = []
        if polygons is not None:
            self.polygons = polygons
            # for polygon in polygons:
            #     self.add_polygon(polygon)
        self.rate = rate
        self.antigen_allocation = {}
        self.starting_gc = None

    def get_patches(self):
        all_patches = []
        bndry = plt_ptch.Polygon([(p.x, p.y) for p in self.boundary.vertices], 
                  edgecolor='green', facecolor='lightgreen')
        all_patches.append(bndry)

        for poly in self.polygons:
            # curr_antigen = self.antigen_allocation.get(poly.name)
            plot_poly = plt_ptch.Polygon([(p.x, p.y) for p in poly.vertices], 
                  edgecolor='blue', facecolor='lightblue')
            all_patches.append(plot_poly)

        return all_patches


    # def make_animation(self, positions, step):
    #     fig, ax = plt.subplots()

    #     bndry = plt_ptch.Polygon([(p.x, p.y) for p in self.boundary.vertices], 
    #               edgecolor='green', facecolor='lightgreen')
    #     ax.add_patch(bndry)

    #     for poly in self.polygons:
    #         plot_poly = plt_ptch.Polygon([(p.x, p.y) for p in poly.vertices], 
    #               edgecolor='blue', facecolor='lightblue')
    #         ax.add_patch(plot_poly)

    #     # curr_positions = []
    #     grouped = positions.groupby("gen")
    #     data_positions = ax.scatter([], [], c="r", s=20, edgecolor='k', lw=1.5, label='B cells')
    #     plt.plot()

    #     def update(i):
    #         curr_positions = [(row["x"], row["y"]) for _, row in grouped.get_group(i).iterrows()]
    #         x = [a for a, _ in curr_positions]
    #         y = [a for _, a in curr_positions]
    #         data_positions.set_offsets(np.stack([x, y]).T)
    #         return data_positions,

    #     ani = FuncAnimation(fig, update, frames=range(0,len(grouped), step), interval=500, blit=True)
    #     plt.show()


    def initialize(self):
        self.create_antigen_allocation()
        max_polygon_name = max(self.antigen_allocation, key=self.antigen_allocation.get)
        max_polygon = [x for x in self.polygons if x.name == max_polygon_name][0]
        self.starting_gc = max_polygon

    def pick_random_root_starting_point(self):
        """Picks a random starting point for a root cell in whichever polygon has the highest antigen allocation."""
        self.initialize()
        while True:
            x = s.RNG.uniform(self.starting_gc.min_x, self.starting_gc.max_x)
            y = s.RNG.uniform(self.starting_gc.min_y, self.starting_gc.max_y)
            point = Point(x, y)
            if self.starting_gc.is_point_in_polygon(point):
                return point
            # polygon = self.get_polygon_containing_point(point)
            # if polygon is not None and polygon.name == self.starting_gc.name:
            #     self.starting_gc = polygon
            #     return point


    def create_antigen_allocation(self):
        """Creates an antigen allocation for the spatial plane."""
        concentration = [10, 8, 6, 4, 2]
        concentration.extend(np.ones(len(self.polygons)))
        concentration = concentration[:len(self.polygons)]
        s.RNG.shuffle(concentration)
        dirichlet = s.RNG.dirichlet(concentration, size=1)[0]
        max_population = [x for x in s.LOCATIONS if x.name == LocationName.GC][0].max_population
        self.antigen_allocation = {
            self.polygons[i].name: int(dirichlet[i] * 2 * max_population)
            for i in range(len(self.polygons))
            }
        logger.warning(self.antigen_allocation)


    def get_polygon_containing_point(self, point):
        """Returns the polygon containing a given point."""
        for polygon in self.polygons:
            if polygon.is_point_in_polygon(point):
                return polygon
        return None

    def do_outside_random_walk(self, point):
        """Performs a random walk for a point outside the spatial plane."""
        # Simulate brownian motion with reflection at the boundary
        count = 0
        while True:
            count += 1
            new_point = self.boundary.do_random_walk(point, self.rate)
            if self.get_polygon_containing_point(new_point) is None:
                return new_point


    def do_entry_random_walk(self, point, last_GC=None):
        """Performs a random walk for a point that should enter a polygon."""
        count = 0
        if last_GC is None:
            raise ValueError("last_polygon was none for a cell, which is not allowed")
        while True:
            count += 1
            new_point = self.boundary.do_random_walk(point, self.rate)
            new_poly = self.get_polygon_containing_point(new_point)
            if new_poly is not None and new_poly.name != last_GC.name:
                return new_point


    def do_random_walk(self, point, transition=False, last_GC=None):
        """Performs a random walk for a point on the spatial plane."""
        polygon = self.get_polygon_containing_point(point)
        if not transition:
            if polygon is None:
                return self.do_outside_random_walk(point)
            else:
                return polygon.do_random_walk(point, self.rate)
        else:
            if polygon is None:
                # if polygon is none, it's outside, and if not inside, then that 
                # means transition
                return self.do_entry_random_walk(point, last_GC=last_GC)
            else:
                return self.do_outside_random_walk(point)

    @classmethod
    def from_dict(cls, json_dict):
        """Creates a SpatialPlane instance from a dictionary."""
        boundary = Polygon.from_dict(json_dict["boundary"])
        rate = json_dict["rate"]
        polygons = [Polygon.from_dict(poly) for poly in json_dict["polygons"]]
        return cls(boundary, rate, polygons)

class Polygon(Encodable):
    """Represents a polygon in the spatial plane.

    Attributes:
    """
    def __init__(
            self,
            name,
            vertices
    ):
        """Initializes a Polygon instance.
        
        Args:
            vertices (list): A list of vertices defining the polygon.
     """
        self.name = name
        self.vertices = self.sort_vertices(vertices)
        self.edges = []
        self.build_edges()
        self.min_x = min(vertex.x for vertex in self.vertices)
        self.max_x = max(vertex.x for vertex in self.vertices)
        self.min_y = min(vertex.y for vertex in self.vertices)
        self.max_y = max(vertex.y for vertex in self.vertices)

    def build_edges(self):
        """Builds the edges of the polygon from its vertices."""
        n = len(self.vertices)
        for i in range(n):
            self.edges.append(Line(self.vertices[i], self.vertices[(i + 1) % n]))

    def sort_vertices(self, vertices):
        """Sorts the vertices of the polygon in a consistent order."""
        # Sort vertices by angle from centroid
        centroid_x = sum(vertex.x for vertex in vertices) / len(vertices)
        centroid_y = sum(vertex.y for vertex in vertices) / len(vertices)
        return sorted(vertices, key=lambda v: np.arctan2(v.y - centroid_y, v.x - centroid_x))

    def is_point_in_bounding_box(self, point):
        """Determines if a point is within the bounding box of the polygon."""
        return self.min_x <= point.x <= self.max_x and self.min_y <= point.y <= self.max_y

    def is_point_in_polygon(self, point):
        """Determines if a point is inside the polygon."""
        # Implementing the ray-casting algorithm to determine if the point is inside the polygon
        if not self.is_point_in_bounding_box(point):
            return False

        inside = False

        # make a line from the point in any direction but make sure the other end of
        # the line is definitely outside the polygon -- easy way to do that is
        # make sure it is outside the bounding box
        ray = Line(point, Point(self.max_x + 10, self.max_y + 10))

        for edge in self.edges:
            intersection = edge.get_intersection_point(ray)
            if intersection is not None:
                if (intersection == edge.start and edge.end.y < intersection.y or 
                    intersection == edge.end and edge.start.y < intersection.y or
                    intersection != edge.start and intersection != edge.end):
                    inside = not inside

        return inside
    
    # def is_point_in_polygon_debug(self, point):
    #     """Determines if a point is inside the polygon."""
    #     # Implementing the ray-casting algorithm to determine if the point is inside the polygon
    #     if not self.is_point_in_bounding_box(point):
    #         return False

    #     inside = False

    #     print(f"in is_point_in_polygon_debug with point {point}", flush=True)


    #     # make a line from the point in any direction but make sure the other end of
    #     # the line is definitely outside the polygon -- easy way to do that is
    #     # make sure it is outside the bounding box
    #     ray = Line(point, Point(self.max_x + 10, self.max_y + 10))

    #     plt.plot([ray.start.x, ray.end.x], [ray.start.y, ray.end.y], "x-k")
    #     # plt.show()

    #     for edge in self.edges:
    #         intersection = edge.get_intersection_point_debug(ray)
    #         if intersection is not None:
    #             plt.plot([edge.start.x, edge.end.x], [edge.start.y, edge.end.y], "x-r")
    #             if (intersection == edge.start and edge.end.y < intersection.y or 
    #                 intersection == edge.end and edge.start.y < intersection.y or
    #                 intersection != edge.start and intersection != edge.end):
    #                 inside = not inside
        
    #     plt.show()

    #     return inside

    def do_random_walk(self, point, rate):
        """Performs a random walk for a point within the polygon."""
        # Simulate brownian motion
        count = 0
        while True:
            count += 1
            new_x = s.RNG.normal(point.x, np.sqrt(1*rate))
            new_y = s.RNG.normal(point.y, np.sqrt(1*rate))
            new_point = Point(new_x, new_y)
            if self.is_point_in_polygon(new_point):
                return new_point
            
            # if (count > 5000 and count % 5000 == 0):
            #     suggestions.append(new_point)
            # if (count > 50000 and count % 50000 == 0):
            #     print(f"\tpolygon id: {self.name} \tstarting point: {point} \tsuggested point: {new_point} \tcount: {count}", flush=True)
            #     print(f"start in box?? {self.is_point_in_bounding_box(point)} \t start in poly?? {self.is_point_in_polygon(point)}")
            #     print(f"curr in box?? {self.is_point_in_bounding_box(new_point)} \t curr in poly?? {self.is_point_in_polygon(new_point)}")
            #     x_sugs = [p.x for p in suggestions]
            #     y_sugs = [p.y for p in suggestions]
            #     plt.plot(x_sugs, y_sugs, 'ro')
            #     plt.plot(point.x, point.y, 'k+')
            #     plt.plot(new_point.x, new_point.y, 'ks')
            #     plt.show()
                

    @classmethod
    def from_dict(cls, json_dict):
        """Creates a Polygon instance from a dictionary."""
        name = json_dict["name"]
        vertices = [Point(vertex["x"], vertex["y"]) for vertex in json_dict["vertices"]]
        return cls(name, vertices)


class Point(Encodable):
    """Represents a point in the spatial plane.

    Attributes:
    """
    def __init__(
            self,
            x,
            y
    ):
        """Initializes a Point instance.
        
        Args:
            x (float): The x-coordinate of the point.
            y (float): The y-coordinate of the point.
     """
        self.x = x
        self.y = y

    def __eq__(self, other): 
        if not isinstance(other, Point):
            # don't attempt to compare against unrelated types
            raise NotImplementedError

        return self.x == other.x and self.y == other.y
    
    def __repr__(self):
        return f"{self.x}, {self.y}"

class Line(Encodable):
    """Represents a line in the spatial plane.

    Attributes:
    """
    def __init__(
            self,
            a,
            b
    ):
        """Initializes a Line instance.
        
        Args:
            start (Point): The starting point of the line.
            end (Point): The ending point of the line.
     """
        if a.x < b.x or a.x == b.x and a.y < b.y:
            start, end = a, b
        else:
            start, end = b, a
        self.start = start
        self.end = end
        self.slope = (end.y - start.y) / (end.x - start.x) if (end.x - start.x) != 0 else float('inf')

    def get_slope(self):
        """Returns the slope of the line."""
        return self.slope

    def get_intersection_point(self, other):
        """Calculates the intersection point between two lines, if it exists."""
        a_x = self.start.x
        a_y = self.start.y
        b_x = self.end.x
        b_y = self.end.y
        c_x = other.start.x
        c_y = other.start.y
        d_x = other.end.x
        d_y = other.end.y

        # skip some computation if they have no overlapping x or y values
        if (b_x < c_x or 
            d_x < a_x or 
            max(a_y, b_y) < min(c_y, d_y) or 
            max(c_y, d_y) < min(a_y, b_y)):
            return None

        slope_ab = self.get_slope()
        slope_cd = other.get_slope()

        if slope_ab == slope_cd:
            return None

        if slope_ab == float('inf'):
            x = a_x
            y = slope_cd * (x - c_x) + c_y
        elif slope_cd == float('inf'):
            x = c_x
            y = slope_ab * (x - a_x) + a_y
        else:
            x = ((slope_ab * a_x) - a_y - (slope_cd * c_x) + c_y) / (slope_ab - slope_cd)
            if not self.is_point_in_segment(x=x):
                return None
            y = slope_ab * (x - a_x) + a_y

        if not (self.is_point_in_segment(x=x, y=y) and other.is_point_in_segment(x=x, y=y)):
            return None
        
        return Point(x, y)
    
    def is_point_in_segment(self, x=None, y=None, point=None):
        if point is None and x is None and y is None:
                raise ValueError("nothing specified")
        if point is not None:
            if x is not None or y is not None:
                raise ValueError("x and y cannot be specified if point is provided")
            return self.is_point_in_segment(x=point.x, y=point.y)
        x_in_segment = min(self.start.x, self.end.x) <= x <= max(self.start.x, self.end.x) if x is not None else True
        y_in_segment = min(self.start.y, self.end.y) <= y <= max(self.start.y, self.end.y) if y is not None else True
        return x_in_segment and y_in_segment


def create_spatial_plane_from_csv(file_path, rate):
    """Creates a SpatialPlane instance from a CSV file.

    Args:
        file_path (str): The path to the CSV file containing polygon data.

    Returns:
        SpatialPlane: A SpatialPlane instance containing the polygons defined in the CSV file.
    """
    df = pd.read_csv(file_path)
    boundary = None
    polygons = []
    for polygon, polygon_data in df.groupby("polygon"):
        current_polygon_vertices = []
        if polygon == "bp":
            for _, vertex_row in polygon_data.iterrows():
                current_polygon_vertices.append(Point(vertex_row['x'], vertex_row['y']))
            boundary = Polygon(polygon, current_polygon_vertices)
            continue

        for _, vertex_row in polygon_data.iterrows():
            current_polygon_vertices.append(Point(vertex_row['x'], vertex_row['y']))
        polygons.append(Polygon(polygon, current_polygon_vertices))

    return SpatialPlane(boundary=boundary, polygons=polygons, rate=rate)
