import pygmsh
import meshio

with pygmsh.geo.Geometry() as geom:
    geom.add_polygon(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ],
        mesh_size=0.03,
    )
    mesh = geom.generate_mesh(dim=2)




# remove z-coordinate
mesh = meshio.Mesh(mesh.points[:, :2], {"triangle":  mesh.get_cells_type("triangle")})
meshio.svg.write(
    "laser.svg", mesh, float_fmt=".3f", stroke_width="0.1"
)
mesh.write("laser.xml")