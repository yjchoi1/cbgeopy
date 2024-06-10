import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt


def plot_surfaces(*meshes, save_path, resolution=1, points=None):
    """
    Plots multiple 3D mesh surfaces and multiple groups of points in a single plot.

    Parameters:
    meshes (tuple): A tuple of mesh objects, each containing vertices and faces.
    points (dict, optional): A dictionary of points groups, each group with shape (n_points, 3).
                             Example: {id-1: np.array([...]), id-2: np.array([...])}
    """
    # Create the figure
    fig = go.Figure()

    # Define color scales for different meshes
    color_scales = ['Viridis', 'Magma', 'Thermal', 'Icefire']
    markers = ['circle', 'square', 'diamond', 'cross']

    # Initialize lists to track the coordinate ranges
    all_x, all_y, all_z = [], [], []

    # Plot each mesh
    for index, mesh in enumerate(meshes):
        # Extract vertices and faces for the current mesh
        x, y, z = mesh.vertices.T
        i, j, k = mesh.faces.T

        # Update coordinate ranges
        all_x.extend(x)
        all_y.extend(y)
        all_z.extend(z)

        # Add the mesh to the figure
        fig.add_trace(
            go.Mesh3d(
                x=x, y=y, z=z,
                i=i, j=j, k=k,
                intensity=z,
                colorscale=color_scales[index % len(color_scales)],
                opacity=0.5,
                name=f'Layer {index + 1}'))

        # Plot each group of points if provided
    if points is not None:
        for i, particle_group in points.items():
            # Ensure particle group has the correct shape (n_points, 3)
            if particle_group.shape[1] != 3:
                raise ValueError("Each group of particles must have shape (n_points, 3)")

            px, py, pz = particle_group[::resolution].T

            # Add points group to the figure
            fig.add_trace(
                go.Scatter3d(
                    x=px, y=py, z=pz,
                    mode='markers',
                    marker=dict(
                        size=5,
                        symbol=markers[i],
                        line=dict(width=1)
                    ),
                    name=f"group-{i}"))

    # Update layout with the actual aspect ratio
    fig.update_layout(
        title='Comparison of Multiple Layers',
        scene=dict(
            xaxis=dict(title='X'),
            yaxis=dict(title='Y'),
            zaxis=dict(title='Depth'),
            aspectmode='data'
        ))

    fig.write_html(save_path)
    print(f"Plot saved to {save_path}")


def plot_cross_section(particle_groups, plane, location, tolerance=1e-5):
    """
    Plots the particles in the specified plane at a specified location.

    Args:
        particle_groups (dict): Dictionary containing particle groups.
                                Keys are group names, values are (n_particles, 3) arrays.
        plane (str): The plane to plot ('xy', 'xz', or 'yz').
        location (float): The location on the remaining axis to filter particles.
        tolerance (float): The tolerance range for filtering location (default is 1e-5).
    """
    """
    Plots the particles in the specified plane at a specified location.
    Args:
        particle_groups (dict): Dictionary containing particle groups.
                                Keys are group names, values are (n_particles, 3) arrays.
        plane (str): The plane to plot ('xy', 'xz', or 'yz').
        location (float): The location on the remaining axis to filter particles.
        tolerance (float): The tolerance range for filtering location (default is 1e-5).
    """
    plt.figure(figsize=(5, 3.5))
    markers = ['o', 's', '^', 'd', 'v', '<', '>', 'p', '*', 'h', 'H', 'D']  # Marker styles
    colors = ['blue', 'red', 'green', 'purple', 'orange', 'cyan', 'magenta', 'brown', 'gray', 'olive', 'pink', 'teal']  # Marker edge colors
    for i, (group_name, particles) in enumerate(particle_groups.items()):
        marker = markers[i % len(markers)]  # Assign a unique marker style for each group
        color = colors[i % len(colors)]  # Assign a unique edge color for each group
        if plane == 'xy':
            filtered_particles = particles[
                (particles[:, 2] >= location - tolerance) & (particles[:, 2] <= location + tolerance)]
            if filtered_particles.size > 0:
                plt.scatter(filtered_particles[:, 0], filtered_particles[:, 1], marker=marker, facecolors='none', edgecolors=color, label=group_name)
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.title(f'Particles in X-Y plane at Z = {location}')
        elif plane == 'xz':
            filtered_particles = particles[
                (particles[:, 1] >= location - tolerance) & (particles[:, 1] <= location + tolerance)]
            if filtered_particles.size > 0:
                plt.scatter(filtered_particles[:, 0], filtered_particles[:, 2], marker=marker, facecolors='none', edgecolors=color, label=group_name)
            plt.xlabel('X')
            plt.ylabel('Z')
            plt.title(f'Particles in X-Z plane at Y = {location}')
        elif plane == 'yz':
            filtered_particles = particles[
                (particles[:, 0] >= location - tolerance) & (particles[:, 0] <= location + tolerance)]
            if filtered_particles.size > 0:
                plt.scatter(filtered_particles[:, 1], filtered_particles[:, 2], marker=marker, facecolors='none', edgecolors=color, label=group_name)
            plt.xlabel('Y')
            plt.ylabel('Z')
            plt.title(f'Particles in Y-Z plane at X = {location}')
        else:
            raise ValueError("Invalid plane. Choose from 'xy', 'xz', or 'yz'.")
    plt.legend()
    plt.grid(True)
    plt.gca().set_aspect('equal')
    plt.show()
