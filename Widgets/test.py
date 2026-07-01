import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import unary_union
import ipywidgets as widgets
from IPython.display import display, clear_output

# ==========================================
# Geometry Generation Functions
# ==========================================
def create_irregular_shrub(center_x, center_y, base_radius, num_points=12):
    """Generates an irregular polygon to represent a shrub."""
    angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
    points = []
    for angle in angles:
        # Perturb the radius by +/- 40% to make it irregular
        r = base_radius * np.random.uniform(0.6, 1.4)
        points.append((center_x + r * np.cos(angle), center_y + r * np.sin(angle)))
    return Polygon(points)

def generate_random_shrubs(num_shrubs=10, plot_radius=55):
    """Generates a list of randomly placed shrub polygons within the plot."""
    shrubs = []
    for _ in range(num_shrubs):
        # Random radius between 1m and 6m
        radius = np.random.uniform(1.0, 6.0)
        # Random position within the plot (keeping it mostly inside the boundary)
        r = np.random.uniform(0, plot_radius - radius)
        theta = np.random.uniform(0, 2 * np.pi)
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        shrubs.append(create_irregular_shrub(x, y, radius))
    return shrubs

def get_nri_spokes(hub_radius=5, plot_radius=55):
    """Generates the 3 NRI transect lines at 0, 120, and 240 degrees."""
    angles = [0, 120, 240]
    spokes = []
    for a in angles:
        rad = np.radians(a)
        x1, y1 = hub_radius * np.cos(rad), hub_radius * np.sin(rad)
        x2, y2 = plot_radius * np.cos(rad), plot_radius * np.sin(rad)
        spokes.append(LineString([(x1, y1), (x2, y2)]))
    return spokes

# ==========================================
# Main Simulation Class
# ==========================================
class NRISimulation:
    def __init__(self):
        self.plot_radius = 55
        self.hub_radius = 5
        self.plot_circle = Point(0, 0).buffer(self.plot_radius)
        self.spokes = get_nri_spokes(self.hub_radius, self.plot_radius)
        self.total_transect_length = 3 * (self.plot_radius - self.hub_radius) # 150m
        self.shrubs = []
        self.shrub_union = None
        
        # Output widget for matplotlib to prevent cell output spam
        self.out = widgets.Output()
        
        # Setup UI FIRST so the dropdown attribute exists before plotting
        self.setup_ui()
        
        # Generate initial data (this calls update_plot internally)
        self.generate_new_shapes()

    def generate_new_shapes(self, b=None):
        """Creates entirely new shrub geometries and positions."""
        self.shrubs = generate_random_shrubs(10, self.plot_radius)
        self.shrub_union = unary_union(self.shrubs)
        self.update_plot()

    def randomize_positions(self, b=None):
        """Moves existing shrubs to new locations."""
        new_shrubs = []
        for shrub in self.shrubs:
            # Calculate current centroid
            cx, cy = shrub.centroid.x, shrub.centroid.y
            
            # Generate new random position
            radius = np.sqrt(shrub.area / np.pi) # approximate base radius
            r = np.random.uniform(0, self.plot_radius - radius)
            theta = np.random.uniform(0, 2 * np.pi)
            new_cx, new_cy = r * np.cos(theta), r * np.sin(theta)
            
            # Translate polygon
            from shapely.affinity import translate
            new_shrub = translate(shrub, xoff=(new_cx - cx), yoff=(new_cy - cy))
            new_shrubs.append(new_shrub)
            
        self.shrubs = new_shrubs
        self.shrub_union = unary_union(self.shrubs)
        self.update_plot()

    def setup_ui(self):
        self.protocol_dropdown = widgets.Dropdown(
            options=[
                ('0m (Continuous Line)', 0.0),
                ('0.25m (Point Intercept)', 0.25),
                ('0.50m (Point Intercept)', 0.50),
                ('1.00m (Point Intercept)', 1.00),
                ('2.00m (Point Intercept)', 2.00)
            ],
            value=0.0,
            description='Protocol:'
        )
        self.protocol_dropdown.observe(lambda change: self.update_plot() if change['name'] == 'value' else None)

        self.btn_randomize = widgets.Button(description="Randomize Positions", button_style='info')
        self.btn_randomize.on_click(self.randomize_positions)
        
        self.btn_new_shapes = widgets.Button(description="Generate New Shrubs", button_style='warning')
        self.btn_new_shapes.on_click(self.generate_new_shapes)

        controls = widgets.HBox([self.protocol_dropdown, self.btn_randomize, self.btn_new_shapes])
        display(widgets.VBox([controls, self.out]))

    def update_plot(self):
        with self.out:
            clear_output(wait=True)
            interval = self.protocol_dropdown.value
            
            # 1. Calculate True Bare Ground
            # Clip shrubs to plot boundary just in case they stick out
            shrubs_in_plot = self.shrub_union.intersection(self.plot_circle)
            true_bare_area = self.plot_circle.area - shrubs_in_plot.area
            true_bg_pct = (true_bare_area / self.plot_circle.area) * 100

            # 2. Plotting setup
            fig, ax = plt.subplots(figsize=(5, 5))
            ax.set_aspect('equal')
            ax.set_facecolor('#f0f0f0')
            ax.axis('off')
            
            # Draw plot boundary and hub
            plot_patch = plt.Circle((0, 0), self.plot_radius, color='white', ec='black', lw=2)
            hub_patch = plt.Circle((0, 0), self.hub_radius, color='#e0e0e0', ec='black', lw=1, ls='--')
            ax.add_patch(plot_patch)
            ax.add_patch(hub_patch)

            # Draw Shrubs
            if self.shrub_union.geom_type == 'Polygon':
                geoms = [self.shrub_union]
            else:
                geoms = self.shrub_union.geoms
                
            for geom in geoms:
                x, y = geom.exterior.xy
                ax.fill(x, y, alpha=0.5, fc='forestgreen', ec='darkgreen', lw=1.5)

            # 3. Calculate Sampled Bare Ground & Draw Sampling Tools
            sampled_bg_pct = 0
            
            if interval == 0.0: # Continuous Line Intercept
                total_intercept = 0
                for spoke in self.spokes:
                    # Draw base spoke
                    x, y = spoke.xy
                    ax.plot(x, y, color='gray', lw=2, zorder=3)
                    
                    # Intersect with shrubs
                    intersection = spoke.intersection(self.shrub_union)
                    total_intercept += intersection.length
                    
                    # Highlight intersections
                    if not intersection.is_empty:
                        if intersection.geom_type == 'LineString':
                            ix, iy = intersection.xy
                            ax.plot(ix, iy, color='red', lw=3, zorder=4)
                        elif intersection.geom_type == 'MultiLineString':
                            for line in intersection.geoms:
                                ix, iy = line.xy
                                ax.plot(ix, iy, color='red', lw=3, zorder=4)

                sampled_bg_pct = (1 - (total_intercept / self.total_transect_length)) * 100
                
            else: # Point Intercept
                total_points = 0
                bare_points = 0
                
                # Draw faint lines for context
                for spoke in self.spokes:
                    x, y = spoke.xy
                    ax.plot(x, y, color='lightgray', lw=1, zorder=2)
                
                # Generate and test points
                for spoke in self.spokes:
                    # spoke length is 50m
                    for d in np.arange(0, 50.001, interval):
                        pt = spoke.interpolate(d)
                        total_points += 1
                        
                        if self.shrub_union.contains(pt):
                            ax.plot(pt.x, pt.y, 'ro', markersize=4, zorder=4) # Shrub hit
                        else:
                            ax.plot(pt.x, pt.y, 'ko', markersize=3, zorder=4) # Bare ground hit
                            bare_points += 1
                            
                sampled_bg_pct = (bare_points / total_points) * 100 if total_points > 0 else 0

            # 4. Display Metrics
            error = sampled_bg_pct - true_bg_pct
            
            title_text = (
                f"True Bare Ground: {true_bg_pct:.2f}%\n"
                f"Sampled Bare Ground: {sampled_bg_pct:.2f}%\n"
            )
            
            # Color code the error text based on severity
            color = 'black'
            if abs(error) > 5: color = 'darkorange'
            if abs(error) > 10: color = 'red'
            
            plt.title(title_text, fontsize=14, pad=15)
            plt.figtext(0.5, 0.05, f"Error (Bias): {error:+.2f}%", ha="center", fontsize=16, fontweight='bold', color=color)
            
            plt.xlim(-60, 60)
            plt.ylim(-60, 60)
            plt.tight_layout()
            plt.show()

# Run the simulation
simulation = NRISimulation()
