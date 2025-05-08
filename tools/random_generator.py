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


def generate_rectangles(
        domain_range, width_range, aspect_ratio_range, num_rectangles, tolerance=0):
    """
    Generates 2D rectangles with random width and aspect ratio.
    Height is calculated as width * aspect_ratio.

    Args:
        domain_range (List[List[float]]): [[min_x, max_x], [min_y, max_y]]
        width_range (List[float]): [min_width, max_width]
        aspect_ratio_range (List[float]): [min_aspect_ratio, max_aspect_ratio]
        num_rectangles (int): Number of rectangles to generate.
        tolerance (float): Tolerance for checking if a rectangle is contained within another.

    Returns:
        List[Tuple[List[float], List[float]]]: A list of rectangles, where each
                                                rectangle is (origin, [width, height]).
    """
    n_dims = 2  # This function is specifically for 2D
    if len(domain_range) != n_dims:
        raise ValueError(f"Domain range must be 2D, got {len(domain_range)}D")
    for i in range(n_dims):
        if len(domain_range[i]) != 2:
            raise ValueError(f"Domain range for dimension {i} must be [min, max]")

    rectangles = []

    # Validate inputs
    if width_range[0] <= 0 or width_range[1] <= 0:
        raise ValueError("Width range values must be positive.")
    if width_range[0] > width_range[1]:
        raise ValueError("min_width cannot be greater than max_width.")
    if aspect_ratio_range[0] <= 0 or aspect_ratio_range[1] <= 0:
        raise ValueError("Aspect ratio range values must be positive.")
    if aspect_ratio_range[0] > aspect_ratio_range[1]:
        raise ValueError("min_aspect_ratio cannot be greater than max_aspect_ratio.")

    min_height_possible = width_range[0] * aspect_ratio_range[0]
    if domain_range[0][1] - domain_range[0][0] < width_range[0]:
        raise ValueError(f"Domain range in x-dimension is too small for min_width {width_range[0]}")
    if domain_range[1][1] - domain_range[1][0] < min_height_possible:
        raise ValueError(f"Domain range in y-dimension is too small for min_height {min_height_possible}")

    def is_fully_contained(rect1, rect2, tolerance):
        # Returns True if rect1 is fully contained within rect2, with some tolerance
        origin1, dims1 = rect1  # dims1 is [width1, height1]
        origin2, dims2 = rect2  # dims2 is [width2, height2]
        for i in range(n_dims):
            if not (origin1[i] >= origin2[i] - tolerance and
                    origin1[i] + dims1[i] <= origin2[i] + dims2[i] + tolerance):
                return False
        return True

    attempts = 0
    max_attempts = num_rectangles * 1000  # To prevent infinite loops

    while len(rectangles) < num_rectangles and attempts < max_attempts:
        # Generate random width and aspect ratio
        width = random.uniform(width_range[0], width_range[1])
        aspect_ratio = random.uniform(aspect_ratio_range[0], aspect_ratio_range[1])
        height = width * aspect_ratio

        rect_dims = [width, height]

        # Ensure that we can generate a valid origin
        valid_origin_possible = True
        for i in range(n_dims):
            if domain_range[i][1] - domain_range[i][0] < rect_dims[i]:
                valid_origin_possible = False
                break
        
        if not valid_origin_possible:
            attempts += 1
            continue

        # Generate random origin within the domain
        origin = [random.uniform(domain_range[i][0], domain_range[i][1] - rect_dims[i]) for i in range(n_dims)]

        new_rectangle = (origin, rect_dims)

        # Check if the new rectangle is fully contained within any existing rectangle, or vice versa
        contained = False
        for rect in rectangles:
            if is_fully_contained(new_rectangle, rect, tolerance) or \
               is_fully_contained(rect, new_rectangle, tolerance):
                contained = True
                break

        # If not contained, add to list
        if not contained:
            rectangles.append(new_rectangle)

        attempts += 1

    if len(rectangles) < num_rectangles:
        # Changed to a warning as per the original generate_cubes, but consider if an error is more appropriate
        print(f"Warning: Only generated {len(rectangles)} rectangles after {attempts} attempts. Requested {num_rectangles}.")
        # raise ValueError(f"Warning: Only generated {len(rectangles)} rectangles after {attempts} attempts.")
    return rectangles


def generate_soils(
        n_soil_range: List,
        friction_range: List=None,
        friction_options: List=None,
        cohesion_range: List=None,
        cohesion_options: List=None
):
    """
    Randomly generate materials with specified number of soils and material property range.
    Note that the material id starts from 1. This is to accommodate the bedrock id as 0.
    Args:
        n_soil_range (int): [min, max]
        friction_range (List): [min_friction, max_friction]
        friction_options (List): [a list of friction angles to choose from]
        cohesion_range (List): [min_cohesion, max_cohesion]
        cohesion_options (List): [a list of cohesions to choose from]

    Returns:
        Dict mpm input for materials
    """
    n_soils = random.randint(*n_soil_range)
    soils = []
    for i in range(1, n_soils + 1):
        # Friction
        friction = None
        if friction_range is not None:
            friction = round(random.uniform(*friction_range), 2)
        elif friction_options is not None:
            friction = random.choice(friction_options)
        else:
            friction = None

        # Cohesion
        cohesion = None
        if cohesion_range is not None:
            cohesion = round(random.uniform(*cohesion_range), 2)
        elif cohesion_options is not None:
            cohesion = random.choice(cohesion_options)
        else:
            cohesion = None
            
        if friction is None and cohesion is None:
            raise ValueError("Friction and cohesion cannot be both None")

        soil = {
            "id": i,
            "density": 1800,
            "youngs_modulus": 40000000.0,
            "poisson_ratio": 0.3,
            "friction": round(friction, 2) if friction is not None else 0,
            "dilation": 0.0,
            "cohesion": round(cohesion) if cohesion is not None else 100,
            "tension_cutoff": 10,
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




