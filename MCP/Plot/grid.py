import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Circle, Rectangle
from matplotlib.path import Path

# Parameters
radius = 1.0         # circle radius
n = 6                # number of grid cells per side
gap = 0.1            # spacing of grid holes
hole_size = 0.2      # hole size

# Big outer rectangle (covers whole plot)
outer = Rectangle((-2,-2), 4, 4)  # much larger than circle

verts = []
codes = []

# Add outer rectangle first
path = outer.get_path().transformed(outer.get_transform())
verts.extend(path.vertices)
codes.extend(path.codes)

# Subtract circle (keep transparent)
circle = Circle((0,0), radius)
path = circle.get_path().transformed(circle.get_transform())
verts.extend(path.vertices[::-1])   # reverse direction for subtraction
codes.extend(path.codes[::-1])

# Add back circle grid bars (so they stay filled, not cut out)
for i in range(-n//2, n//2):
    for j in range(-n//2, n//2):
        x, y = i*(hole_size+gap), j*(hole_size+gap)
        square = Rectangle((x-hole_size/2, y-hole_size/2), hole_size, hole_size)
        sq_path = square.get_path().transformed(square.get_transform())
        verts.extend(sq_path.vertices)
        codes.extend(sq_path.codes)

compound_path = Path(verts, codes)

# Plot
fig, ax = plt.subplots()
patch = PathPatch(compound_path, facecolor="slateblue", alpha=0.5, lw=0)
ax.add_patch(patch)

ax.set_aspect("equal")
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)
plt.show()
