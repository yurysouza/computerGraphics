# algorithms.py
from grid import Grid

def my_render_cells_algorithm(selected_cells, rendered_cells, parameters, grid):
    for cell in selected_cells:
        grid.render_cell(cell)

def translate(selected_cells, rendered_cells, parameters, grid):
    x_offset = int(parameters['x'])
    y_offset = int(parameters['y'])
    min_x = min(cell[0] for cell in selected_cells)
    max_x = max(cell[0] for cell in selected_cells)
    min_y = min(cell[1] for cell in selected_cells)
    max_y = max(cell[1] for cell in selected_cells)
    for cell in rendered_cells:
        if min_x <= cell[0] <= max_x and min_y <= cell[1] <= max_y:
            grid.clear_cell(cell)
            new_cell = (cell[0] + x_offset, cell[1] + y_offset)
            grid.render_cell(new_cell)

def bresenham(selected_cells, rendered_cells, parameters, grid):
    (x1, y1), (x2, y2) = selected_cells[0], selected_cells[1]

    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    sx = 1 if x1 < x2 else -1
    sy = 1 if y1 < y2 else -1
    er1 = dx - dy

    while True:
        grid.render_cell((x1, y1))
        if x1 == x2 and y1 == y2:
            break
        er2 = er1
        if er2 > -dy:
            er1 -= dy
            x1 += sx
        if er2 < dx:
            er1 += dx
            y1 += sy


def circle(selected_cells, rendered_cells, radius, grid):
    x0, y0 = selected_cells[0]
    x = 0  # Initialize x to 0
    y = int(radius['r'])  # y is set to the radius
    
    e = 1 - int(radius['r'])  # Decision parameter

    points = []

    while x <= y:  # Change the condition to x <= y
        # Add circle points
        points.append((x0 + x, y0 + y))
        points.append((x0 + y, y0 + x))
        points.append((x0 - y, y0 + x))
        points.append((x0 - x, y0 + y))
        points.append((x0 - x, y0 - y))
        points.append((x0 - y, y0 - x))
        points.append((x0 + y, y0 - x))
        points.append((x0 + x, y0 - y))

        x += 1  # Increment x
        if e <= 0:
            e += 2 * x + 1  # Mid-point is inside or on the circle
        else:
            y -= 1  # Decrease y
            e += 2 * x - 2 * y + 1  # Mid-point is outside the circle

    for point in points:
        grid.render_cell(point)  # Render the calculated points

def elipse(selected_cells, render_cells, parameters, grid):
    x0, y0 = selected_cells[0]
    a = int(parameters['a'])
    b = int(parameters['b'])
    x = 0
    y = b
    e = b**2 - (a**2 * b) + (0.25 * a**2)

    points = []

    # First half of the ellipse
    while (2 * b**2 * x) <= (2 * a**2 * y):
        points.append((x0 + x, y0 + y))
        points.append((x0 - x, y0 + y))
        points.append((x0 + x, y0 - y))
        points.append((x0 - x, y0 - y))

        x += 1
        if e < 0:
            e += 2 * b**2 * x + b**2
        else:
            y -= 1
            e += 2 * b**2 * x - 2 * a**2 * y + b**2

    # Second half of the ellipse
    x = a
    y = 0
    e = a**2 - (b**2 * a) + (0.25 * b**2)

    while (2 * a**2 * y) <= (2 * b**2 * x):
        points.append((x0 + x, y0 + y))
        points.append((x0 - x, y0 + y))
        points.append((x0 + x, y0 - y))
        points.append((x0 - x, y0 - y))

        y += 1
        if e < 0:
            e += 2 * a**2 * y + a**2
        else:
            x -= 1
            e += 2 * a**2 * y - 2 * b**2 * x + a**2

    for point in points:
        grid.render_cell(point)  # Render the calculated points

def bezier_quadratic(selected_cells, render_cells, parameters, grid):
    num_points=100
    p0 = selected_cells[0] #start point
    p1 = selected_cells[1] #control point
    p2 = selected_cells[2] #end point

    points = []

    for t in range(num_points + 1):
        t /= num_points
        x = (1 - t) ** 2 * p0[0] + 2 * (1 - t) * t * p1[0] + t ** 2 * p2[0]
        y = (1 - t) ** 2 * p0[1] + 2 * (1 - t) * t * p1[1] + t ** 2 * p2[1]
        points.append((int(x), int(y)))

    for point in points:
        grid.render_cell(point)  # Render the calculated points

def bezier_cubic(selected_cells, render_cells, parameters, grid):
    num_points = 100
    p0 = selected_cells[0]  # Start point
    p1 = selected_cells[1]  # Control point 1
    p2 = selected_cells[2]  # Control point 2
    p3 = selected_cells[3]  # End point

    points = []

    for t in range(num_points + 1):
        t /= num_points
        # Cubic Bézier formula
        x = (1 - t) ** 3 * p0[0] + 3 * (1 - t) ** 2 * t * p1[0] + 3 * (1 - t) * t ** 2 * p2[0] + t ** 3 * p3[0]
        y = (1 - t) ** 3 * p0[1] + 3 * (1 - t) ** 2 * t * p1[1] + 3 * (1 - t) * t ** 2 * p2[1] + t ** 3 * p3[1]
        points.append((int(x), int(y)))

    for point in points:
        grid.render_cell(point)  # Render the calculated points

def polilines(selected_cells, render_cells, parameters, grid):
    n = len(selected_cells)
    while n > 1:  
        selected_pair = (selected_cells[n - 1], selected_cells[n - 2])  
        bresenham(selected_pair, render_cells, parameters, grid)  
        n -= 1  

def scanline(selected_cells, render_cells, parameters, grid):
    new_color = '#0068ca'
    start_point = selected_cells[0]
    x_start = start_point
    
    # Extract the polygon vertices
    polygon = render_cells

    # Extract y-coordinates of the vertices and find the bounding box
    y_min = min(y for _, y in polygon)
    y_max = max(y for _, y in polygon)

    # Loop through each scanline from y_min to y_max
    for y in range(y_min, y_max + 1):
        intersections = []

        # Find intersections with polygon edges
        for i in range(len(polygon)):
            p1 = polygon[i]
            p2 = polygon[(i + 1) % len(polygon)]

            # Check if the edge crosses the scanline
            if (p1[1] <= y < p2[1]) or (p2[1] <= y < p1[1]):
                x = p1[0] + (y - p1[1]) * (p2[0] - p1[0]) / (p2[1] - p1[1])
                intersections.append(x)

        intersections.sort()

        for i in range(0, len(intersections), 2):
            if i + 1 < len(intersections): 
                x_start = int(intersections[i])
                x_end = int(intersections[i + 1])

                for x in range(x_start, x_end + 1):                    
                    grid.color_cell((x,y),new_color)

