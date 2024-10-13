import pygmsh
import meshio
from numpy import min

def save_rect_mesh(name: str, resolution: float):
    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(
            [
                [-2.0, 0.0],
                [1.0, 0.0],
                [1.0, 1.0],
                [2.0, 1.0],
                [2.0, 0.0],
                [3.0, 0.0],
                [3.0, 2.0],
                [4.0, 2.0],
                [4.0, 0.0],
                [8.0, 0.0],
                [8.0, 2.5],
                [8.0, 3.0],
                [-2.0, 3.0],
                [-2.0, 0.5],
            ],
            mesh_size=0.1/resolution,
        )
        mesh = geom.generate_mesh()
    # remove z-coordinate
    mesh = meshio.Mesh(mesh.points[:, :2], {"triangle":  mesh.get_cells_type("triangle")})
    meshio.svg.write(
        f"{name}.svg", mesh, float_fmt=".3f", stroke_width="0.1"
    )
    mesh.write(f"{name}.xml")

if __name__ == '__main__':
    save_rect_mesh('rect/rect', 1)