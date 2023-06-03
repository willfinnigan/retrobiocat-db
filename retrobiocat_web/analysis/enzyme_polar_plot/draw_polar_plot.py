import numpy as np

from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, output_file, show, output_notebook
from bokeh.models import HoverTool, LassoSelectTool, TapTool, ImageURL, Text
from bokeh.core import properties

from rdkit import Chem
import math
from bokeh.palettes import Spectral10 as palette
import itertools

from retrobiocat_web.analysis.enzyme_polar_plot.mol_drawing import create_svg, get_svg_url

def translate(X, Y, angle, distance):
    dY = distance * math.sin(angle)  # change in y
    dX = distance * math.cos(angle)  # change in x
    Xfinal = X + dX
    Yfinal = Y + dY
    return Xfinal, Yfinal

def get_svg_data(polar_dict, wedge_edges):
    data_svg = {'url': [], 'x': [], 'y': []}
    for j, v in enumerate(polar_dict.values()):
        smi = v[1]
        mol = Chem.MolFromSmiles(smi)
        svg = create_svg(mol)
        data_svg['url'].append(get_svg_url(svg))
        data_svg['x'].append(wedge_edges[j][0])
        data_svg['y'].append(wedge_edges[j][1])

    return data_svg

def ensure_min_y(ys, min_y=3):
    new_ys = []
    for y in ys:
        if y < min_y:
            new_ys.append(min_y)
        else:
            new_ys.append(y)
    return new_ys

def create_polar_plot(polar_dict, inner_radius=25, outer_radius=100, gap_size=0, num_axes=0,
                      center=(0, 0), mol_dims=(50, 50), mols_to_outer=False, min_y_value=5,
                      enzyme_name=None):

    from bokeh.plotting import figure

    outer_radius += inner_radius

    ys = [v[0] for v in polar_dict.values()]
    ys = ensure_min_y(ys, min_y=min_y_value)
    ys = [y+inner_radius for y in ys]

    width = outer_radius * 1.7
    height = outer_radius * 1.7

    p = figure(plot_width=500, plot_height=500, title="",
               x_axis_type=None, y_axis_type=None,
               x_range=(-width, width), y_range=(-height, height),
               min_border=0)

    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None

    up = np.pi * 0.5
    segment = ((2 * np.pi) - gap_size) / (len(ys))
    initial_start = up + gap_size / 2
    starts = []
    ends = []
    wedge_edges = []
    start = initial_start
    min_height = 0.7 * outer_radius
    for y in ys:
        end = start + segment
        starts.append(start)
        middle = start + 0.5 * segment
        if mols_to_outer == True:
            mol_distance = outer_radius * 1.2
        else:
            mol_distance = y + outer_radius * 0.2
            if mol_distance < min_height:
                mol_distance = min_height
        eX, eY = translate(center[0], center[1], middle, mol_distance)
        wedge_edges.append([eX, eY])
        ends.append(end)
        start = end

    radii = np.linspace(inner_radius * 1.01, outer_radius * 1.01, num_axes)
    p.circle(0, 0, radius=radii, fill_color=None, line_color="black")

    # create a color iterator
    colors = itertools.cycle(palette)
    fill_colours = []
    for y, colour in zip(ys, colors):
        fill_colours.append(colour)

    source = ColumnDataSource(data={'ys': ys,
                                    'starts': starts,
                                    'ends': ends,
                                    'fill_colour': fill_colours})

    p.annular_wedge(0, 0, inner_radius, 'ys', 'starts', 'ends', source=source, line_color='white',
                    fill_color='fill_colour')
    if gap_size != 0:
        space = 0.1 * segment
        p.annular_wedge(0, 0, inner_radius * 1.01, outer_radius * 1.01, initial_start - gap_size + space,
                        initial_start - space, line_color='black', color='white')
        p.annular_wedge(0, 0, outer_radius, outer_radius * 1.1, initial_start - gap_size + space, initial_start - space,
                        line_color=None, fill_color='white')

    mol_height = mol_dims[0]
    mol_width = mol_dims[1]
    svg_data = get_svg_data(polar_dict, wedge_edges)
    source = ColumnDataSource(svg_data)
    image1 = ImageURL(url='url', x='x', y='y', anchor='center', h=mol_height, w=mol_width)
    p.add_glyph(source, image1)

    if enzyme_name != None:
        text_glyph = Text(x=center[0], y=center[1], text=properties.value(enzyme_name),
                          angle=0, text_font_size='14px',
                          text_align='center', text_baseline='middle',
                          text_font=properties.value("Arial"))  # text_color='colour'
        p.add_glyph(source, text_glyph)

    return p

