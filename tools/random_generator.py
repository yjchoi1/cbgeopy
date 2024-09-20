import random
from typing import List
import numpy as np

def generate_cubes(
        domain_range, min_length, max_length, num_cubes, tolerance=0):
    n_dims = len(domain_range)
    cubes = []

    # Validate inputs
    for i in range(n_dims):
        if min_length[i] > max_length[i]:
            raise ValueError(f"min_length[{i}] cannot be greater than max_length[{i}]")
        if domain_range[i][1] - domain_range[i][0] < min_length[i]:
            raise ValueError(f"Domain range in dimension {i} is too small for min_length[{i}]")

    def is_fully_contained(cube1, cube2, tolerance):
        # Returns True if cube1 is fully contained within cube2, with some tolerance
        origin1, length1 = cube1
        origin2, length2 = cube2
        for i in range(n_dims):
            if not (origin1[i] >= origin2[i] - tolerance and
                    origin1[i] + length1[i] <= origin2[i] + length2[i] + tolerance):
                return False
        return True

    attempts = 0
    max_attempts = num_cubes * 1000  # To prevent infinite loops

    while len(cubes) < num_cubes and attempts < max_attempts:
        # Generate random lengths for the cube
        length = [random.uniform(min_length[i], max_length[i]) for i in range(n_dims)]

        # Ensure that we can generate a valid origin
        valid = True
        for i in range(n_dims):
            if domain_range[i][1] - domain_range[i][0] < length[i]:
                valid = False
                break

        if not valid:
            attempts += 1
            continue

        # Generate random origin within the domain
        origin = [random.uniform(domain_range[i][0], domain_range[i][1] - length[i]) for i in range(n_dims)]

        new_cube = (origin, length)

        # Check if the new cube is fully contained within any existing cube, or vice versa
        contained = False
        for cube in cubes:
            if is_fully_contained(new_cube, cube, tolerance) or is_fully_contained(cube, new_cube, tolerance):
                contained = True
                break

        # If not contained, add to list
        if not contained:
            cubes.append(new_cube)

        attempts += 1

    if len(cubes) < num_cubes:
        raise ValueError(f"Warning: Only generated {len(cubes)} cubes after {attempts} attempts.")
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




