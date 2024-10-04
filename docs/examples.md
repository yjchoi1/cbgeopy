# Examples

## Make 3D topography model from mesh files

This example generates material point method model from external `mesh.obj` files which defines layer boundaries of
ground of mountain.

![output example for 3d](img/mp3d-example.png "Material point output example for 3d slope")

```{literalinclude} ../examples/fundao-3d/input_script.py
:language: python
```

## Make MP model based on 2D lines

This example generates material point method model from user-defined points defining the boundary lines that
distinguishes the layers. 

![output example for 3d](img/mp2d-example.png "Material point output example for 2d slope")

```{literalinclude} ../examples/fundao-2d/input_script.py
:language: python
```

## Make MP model for sand cube collision

This example generates material point method model from two colliding cubes. 
This saves 

![output example for 3d](img/sand_cubes2d.png "Material point output example for 2d slope")

```{literalinclude} ../examples/sand_layers-2d/input_script.py
:language: python
```

