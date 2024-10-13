import pygmsh
import meshio
from numpy import min

def save_laser_mesh(name: str, resolution: float):
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

        geom.set_mesh_size_callback(
        lambda dim, tag, x, y, z, lc: min(
            [2e-2 + 6.0e-1 * ((x - 0.75) ** 2 + (y - 0.75) ** 2), 
            2e-2 + 6.0e-1 * ((x - 0.9) ** 2 + (y - 0.1) ** 2),
            0.06])/resolution/1.5
        )

        mesh = geom.generate_mesh()

    # remove z-coordinate
    mesh = meshio.Mesh(mesh.points[:, :2], {"triangle":  mesh.get_cells_type("triangle")})
    meshio.svg.write(
        f"{name}.svg", mesh, float_fmt=".3f", stroke_width="0.1"
    )
    mesh.write(f"{name}.xml")

if __name__ == '__main__':
    save_laser_mesh('laser/laser', 1)