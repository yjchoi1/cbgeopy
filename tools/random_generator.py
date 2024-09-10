import random
from typing import List
import numpy as np


def generate_cubes(
        domain_range, min_length, max_length, num_cubes
):
    n_dims = len(domain_range)
    cubes = []

    def cubes_fully_contained(cube1, cube2):
        # Extract the origins and lengths
        origin1, length1 = cube1
        origin2, length2 = cube2

        contained = True  # Assume fully contained until proven otherwise

        for i in range(3):  # Check each dimension
            # Check if cube1 is fully inside cube2 along dimension i
            if not (origin1[i] >= origin2[i] and
                    origin1[i] + length1[i] <= origin2[i] + length2[i]):
                contained = False
                break

        return contained

    for _ in range(num_cubes):
        while True:
            # Generate random origin and lengths
            length = [random.uniform(min_length[i], max_length[i]) for i in range(n_dims)]
            origin = [random.uniform(domain_range[i][0], domain_range[i][1] - length[i]) for i in range(n_dims)]

            new_cube = (origin, length)

            # Check if this new cube is fully contained in any existing cube
            fully_contained = False
            for cube in cubes:
                if cubes_fully_contained(new_cube, cube) or cubes_fully_contained(cube, new_cube):
                    fully_contained = True
                    break

            # If it's not fully contained, we accept the new cube
            if not fully_contained:
                cubes.append(new_cube)
                break

    return cubes



def generate_soils(
        n_soil_range: List, friction_range: List
):
    """
    Randomly generate materials with specified number of soils and friction range.
    Note that the material id starts from 1. This is to accommodate the bedrock id as 0.
    Args:
        n_soil_range (int): [min, max]
        friction_range (List): [min_friction, max_friction]

    Returns:
        Dict mpm input for materials
    """
    n_soils = random.randint(*n_soil_range)
    soils = []
    for i in range(1, n_soils + 1):
        soil = {
            "id": i,
            "density": 1800,
            "youngs_modulus": 20000000.0,
            "poisson_ratio": 0.3,
            "friction": round(random.uniform(*friction_range), 2),
            "dilation": 0.0,
            "cohesion": 1000,
            "tension_cutoff": 100,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 30.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb2D"
        }
        soils.append(soil)

    return soils


def generate_bedrock_line(
    x_range: List,
    y_range: List,
    n_middle_points: int
):
    """
    Make a list of points that describes a line
    Args:
        x_range (List): [x_min, x_max]
        y_range (List): [y_min, y_max]
        n_middle_points (int): number of intermediate points

    Returns:
        A list of points
    """
    start = np.array([x_range[0], np.random.uniform(*y_range)])
    end = np.array([x_range[1], np.random.uniform(*y_range)])

    x_points = np.random.uniform(*x_range, n_middle_points)
    x_points = np.sort(x_points)
    y_points = np.random.uniform(*y_range, n_middle_points)
    xy_points = np.stack((x_points, y_points), axis=1)

    line_points = np.concatenate(([start], xy_points, [end]))

    return line_points.tolist()




