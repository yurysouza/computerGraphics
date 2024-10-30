from grid import Grid
from algorithm import *

# Initialize grid
grid = Grid(extent=10, size=500)

# Adds the algorithm to the grid
grid.add_algorithm(name="Render cells", parameters=None, algorithm=lambda cells, rendered, params: my_render_cells_algorithm(cells, rendered, params, grid))
grid.add_algorithm(name='Translate', parameters=['x', 'y'], algorithm=lambda cells, rendered, params: translate(cells, rendered, params, grid))
grid.add_algorithm(name="Bresenham", parameters=None, algorithm=lambda cells, rendered, params: bresenham(cells, rendered, params, grid))
grid.add_algorithm(name="Circle", parameters=['r'], algorithm=lambda cells, rendered, params: circle(cells, rendered, params, grid))
grid.add_algorithm(name="Elipse", parameters=['a','b'], algorithm=lambda cells, rendered, params: elipse(cells, rendered, params, grid))
grid.add_algorithm(name="Bezier Quadratic", parameters=None, algorithm=lambda cells, rendered, params: bezier_quadratic(cells, rendered, params, grid))
grid.add_algorithm(name="Bezier Cubic", parameters=None, algorithm=lambda cells, rendered, params: bezier_cubic(cells, rendered, params, grid))
grid.add_algorithm(name="Polilines", parameters=None, algorithm=lambda cells, rendered, params: polilines(cells, rendered, params, grid))
grid.add_algorithm(name="Scanline", parameters=None, algorithm=lambda cells, rendered, params: scanline(cells, rendered, params, grid))
#grid.add_algorithm(name="Flood Fill", parameters=None, algorithm=lambda cells, rendered, params: flood_fill(cells, rendered, params, grid))



grid.show()
