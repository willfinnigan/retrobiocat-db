import time

from bokeh.colors import RGB

from retrobiocat_web.analysis.data_query import get_data
from rdkit import DataStructs
from bokeh.models.sources import ColumnDataSource
from bokeh.models import HoverTool, ImageURL, Text, TapTool
from bokeh.core.properties import value
from bokeh import palettes
from bokeh.models.callbacks import CustomJS
from bokeh import plotting

from sklearn.metrics.pairwise import pairwise_distances
from scipy.cluster.hierarchy import linkage, dendrogram

from requests.utils import quote

import pandas as pd
import numpy as np


class HeatMapper(object):

    smi_col = 'product_1_smiles'
    enzyme_col = 'enzyme_name'

    block_x_size = 1
    block_y_size = 1
    dendogram_width = 10  # 20-50 is reasonable
    dendogram_height = 5
    scale = 20  # 20-30 is reasonable
    mol_height = 1.1
    mol_width = 2.75
    molSize = (250, 100)
    start_x = (block_x_size * 0.75) + mol_width
    start_y = (block_y_size * 0.75)
    y_title_height = 4
    y_bottom_space = 5
    empty_alpha = 0.1
    data_alpha = 0.8
    hover_alpha = 0.8
    active_colour = RGB(44,181,0)
    inactive_colour = RGB(135,12,0)
    medium_colour = RGB(232,236,0)
    low_colour = RGB(203, 138, 38)
    empty_colour = RGB(191, 191, 191)
    show_toolbar = False
    min_height = 600
    min_width = 700
    enz_letters = 12

    def __init__(self, data_query, log_level=0, plot_mode='binary',
                 smi_col='product_1_smiles'):

        self.log_level = log_level
        self.plot_mode = plot_mode  # alternative is categorical

        self.log("HeatMapper initialised")

        self.smi_col = smi_col
        self.smi_map = {'product_1_smiles': 'Product',
                        'substrate_1_smiles': 'Substrate 1',
                        'substrate_2_smiles': 'Substrate 2'}

        self.activity_data = data_query.activity_data

    def enzyme_hierarchical_clustering(self, enzymes, embeddings_map):
        t0 = time.time()
        if len(embeddings_map) <= 1:
            return [], enzymes, [], [], [], []

        enz_with_embeddings = list(embeddings_map.keys())
        embeddings = list(embeddings_map.values())

        distances = pairwise_distances(embeddings, metric='euclidean')
        linked = linkage(distances, 'ward')
        results = dendrogram(linked,
                             no_plot=True,
                             distance_sort='descending',
                             labels=enz_with_embeddings)

        colours, icoord, dcoord = results['color_list'], results['icoord'], results['dcoord']
        leaf_colours = results['leaves_color_list']
        ordered_enz = results['ivl']
        other_enzymes = [enz for enz in enzymes if enz not in ordered_enz]

        t1 = time.time()
        self.log(f"Enzyme hierarchical clustering completed in {round(t1-t0, 2)} seconds. "
                 f"{len(ordered_enz)} out of {len(enzymes)} enzymes with embeddings could be clustered", 1)
        return ordered_enz, other_enzymes, icoord, dcoord, colours, leaf_colours

    def substrate_hierarchical_clustering(self, substrates, fps):

        if len(substrates) <= 1:
            return substrates, [], [], []

        size = len(substrates)
        hmap = np.empty(shape=(size, size))
        for i, fp_i in enumerate(fps):
            similarities = DataStructs.BulkTanimotoSimilarity(fp_i, list(fps))
            for j, sim in enumerate(similarities):
                hmap[i, j] = sim

        linked = linkage(hmap, 'single')
        results = dendrogram(linked, no_plot=True, distance_sort='ascending')

        colours, icoord, dcoord = results['color_list'], results['icoord'], results['dcoord']
        sub_numbered = list(map(int, results['ivl']))

        substrates_ordered = []
        for i in sub_numbered:
            substrates_ordered.append(substrates[i])

        return substrates_ordered, icoord, dcoord, colours

    def get_activity_info(self, filtered_activity_dict, num_rows, svg_map, substrate):
        cat_colours = {'High': 'green',
                       'Medium': 'goldenrod',
                       'Low': 'orange',
                       'None': 'red'}

        info = ""

        short_citation = filtered_activity_dict.get('short_citation', ['nan']*num_rows)
        active = filtered_activity_dict.get('binary', ['nan']*num_rows)
        categorical = filtered_activity_dict.get('categorical', ['nan']*num_rows)
        conv = filtered_activity_dict.get('conversion', ['nan']*num_rows)
        conv_time = filtered_activity_dict.get('conversion_time', ['nan']*num_rows)
        sa = filtered_activity_dict.get('specific_activity', ['nan']*num_rows)
        reactions = filtered_activity_dict.get('reaction', ['nan']*num_rows)

        svg = svg_map.get(substrate, substrate)
        enzyme_name = filtered_activity_dict.get('enzyme_name', ['nan'])[0]
        enzyme_type = filtered_activity_dict.get('enzyme_type', ['nan'])[0]
        info += f"{svg}"
        info += f"<p>Reaction {self.smi_map[self.smi_col]}</p>"
        info += f"<h5>{enzyme_name}</h3>"
        info += f"<small>{enzyme_type}</small>"
        info += "<hr/>"

        for i, enzyme_name in enumerate(short_citation):
            activity_info = f"<p> {reactions[i]} - "
            if active[i] == 'True':
                activity_info += """<span style="color: green">Active</span>"""
            else:
                activity_info += """<span style="color: red">Inactive</span>"""

            colour = 'darkgrey'
            if str(categorical[i]) != "nan":
                colour = cat_colours.get(categorical[i], 'darkgrey')
                activity_info += f"""<span style="color: {colour}">, {categorical[i]}</span>"""

            if str(conv[i]) != "nan":
                activity_info += f"""<span style="color: {colour}">, {conv[i]} % conversion over {conv_time[i]} hours</span>"""

            if str(sa[i]) != "nan":
                activity_info += f"""<span style="color: {colour}">, {sa[i]} umols/min/mg</span>"""

            activity_info += "</p>"

            info += f"<p>{short_citation[i]}</p>"
            info += activity_info
            info += "<hr/>"

        return info

    def get_best_categorical(self, filtered_df):
        """Take filtered df containing the data for a specific point on heatmap, get best categorical value and return colour"""
        if len(filtered_df[filtered_df['categorical'] == 'High']) != 0:
            return self.active_colour
        if len(filtered_df[filtered_df['categorical'] == 'Medium']) != 0:
            return self.medium_colour
        if len(filtered_df[filtered_df['categorical'] == 'Low']) != 0:
            return self.low_colour
        if len(filtered_df[filtered_df['categorical'] == 'None']) != 0:
            return self.inactive_colour
        if len(filtered_df[filtered_df['binary'] == 'True']) != 0:
            return self.empty_colour

        return self.inactive_colour

    def get_plot_data(self, substrates, enzymes, svg_map):
        activity_df = pd.DataFrame(list(self.activity_data))
        activity_df = activity_df.fillna('nan')

        data = {'info': [], 'xs': [], 'ys': [], 'colors': [], 'ids': []}

        for i, enz in enumerate(enzymes):
            enzyme_df = activity_df[(activity_df[self.enzyme_col] == enz)]
            substrates_for_enzyme = enzyme_df[self.smi_col].unique()
            for j, sub in enumerate(substrates):
                x = self.start_x + (i * self.block_x_size)
                y = self.start_y + (j * self.block_y_size)

                if sub in substrates_for_enzyme:
                    filtered_df = enzyme_df[(enzyme_df[self.smi_col] == sub)]
                    filtered_activity_dict = filtered_df.to_dict(orient='list')
                    num_rows = len(filtered_df.index)
                    data['ids'].append(filtered_activity_dict.get('_id', ['nan']*num_rows))
                    data['info'].append(self.get_activity_info(filtered_activity_dict, num_rows, svg_map, sub))
                    data['xs'].append(x)
                    data['ys'].append(y)
                    if self.plot_mode == 'categorical':
                        # if categorical, need to loop through all the data for this point, and take the highest categorical
                        colour = self.get_best_categorical(filtered_df)
                        data['colors'].append(colour)
                    else:
                        if len(filtered_df[filtered_df['binary'] == 'True']) != 0:
                            data['colors'].append(self.active_colour)
                        else:
                            data['colors'].append(self.inactive_colour)

        return data

    def get_svg_data(self, substrates, svg_url_map):
        data_svg = {'url': [], 'x': [], 'y': []}
        for j, sub in enumerate(substrates):
            y = self.start_y + (j * self.block_y_size)
            data_svg['url'].append(svg_url_map.get(sub, ''))
            data_svg['x'].append(0)
            data_svg['y'].append(y)

        return data_svg

    def get_plot_enzyme_names_data(self, substrates, enzymes, dendogram_enzymes):
        text, text_xs, text_ys, text_colours = [], [], [], []

        y_text = (len(substrates)-1 * self.block_y_size) + self.start_y + 0.75*self.block_y_size

        for i, enz in enumerate(enzymes):
            x = self.block_x_size + (i * self.block_x_size) + self.mol_width
            if enz in dendogram_enzymes:
                if len(enz) >= self.enz_letters:
                    enz = enz[0:self.enz_letters] + '-'
                elif len(enz) < self.enz_letters:
                    enz += '-'*(self.enz_letters-len(enz)) + '-'

            text.append(enz)
            text_xs.append(x)
            text_ys.append(y_text)

        data_titles = {'text': text,
                       'x': text_xs,
                       'y': text_ys}

        return data_titles

    def get_lines(self, substrates, enzymes):

        h_lines = []
        v_lines = []

        y_line_start = self.start_y - (self.block_y_size/2)
        x_line_start = self.start_x - (self.block_x_size/2)

        y_line_end = (self.start_y + (len(substrates) * self.block_y_size)) - (self.block_y_size/2)
        x_line_end = (self.start_x + (len(enzymes) * self.block_x_size)) - (self.block_x_size / 2)

        for i in range(len(enzymes)+1):
            x_point = (self.start_x + (i * self.block_x_size)) - (self.block_x_size / 2)
            xs = [x_point, x_point]
            ys = [y_line_start, y_line_end]
            v_lines.append((xs, ys))

        for j in range(len(substrates)+1):
            y_point = (self.start_y + (j * self.block_y_size)) - (self.block_y_size / 2)
            ys = [y_point, y_point]
            xs = [x_line_start, x_line_end]
            v_lines.append((xs, ys))

        return h_lines, v_lines

    def process_substrate_dendogram_coords(self, icoord, dcoord, y_height):
        icoord = pd.DataFrame(icoord)
        icoord -= icoord.min().min()
        imax = icoord.max().max()
        if imax == 0:
            imax = 1
        icoord = icoord * (y_height / imax)
        icoord += self.start_y
        icoord = icoord.values

        dcoord = pd.DataFrame(dcoord)
        dmax = dcoord.max().max()
        if dmax == 0:
            dmax = 1
        dcoord = dcoord * (self.dendogram_width / dmax)
        dcoord = dcoord.values

        return icoord, dcoord

    def process_enzyme_dendogram_coords(self, icoord, dcoord, x_width, substrates):
        y_text = self.y_title_height + (len(substrates) - 1 * self.block_y_size) + self.start_y + 0.75 * self.block_y_size

        # icoord is x, dcoord is y

        icoord = pd.DataFrame(icoord)
        icoord -= icoord.min().min()
        imax = icoord.max().max()
        if imax == 0:
            imax = 1
        icoord = icoord * (x_width / imax)
        icoord += self.start_x
        icoord = icoord.values

        dcoord = pd.DataFrame(dcoord)
        dcoord -= dcoord.min().min()
        dmax = dcoord.max().max()
        if dmax == 0:
            dmax = 1
        dcoord = dcoord * (self.dendogram_height / dmax)
        dcoord += y_text
        dcoord = dcoord.values

        return icoord, dcoord

    def create_heatmap(self):

        def create_figure(substrates, enzymes):
            x_start = -self.dendogram_width
            x_end = (len(enzymes) * self.block_x_size) + self.start_x + self.mol_width
            width = int((x_end - x_start) * self.scale)
            if width < self.min_width:
                print(f"Width {width} is less than min width of {self.min_width}")
                width = self.min_width
                x_end = (width / self.scale) + x_start

            y_start = -self.y_bottom_space
            y_end = (len(substrates) * self.block_y_size) + self.start_y + self.y_title_height
            y_end += self.dendogram_height
            height = int((y_end - y_start) * self.scale)
            if height < self.min_height:
                print(f"Height {height} is less than min height of {self.min_height}")
                height = self.min_height
                y_start = -1*((height / self.scale)-y_end)

            hm = plotting.figure(x_range=[x_start, x_end],
                                y_range=[y_start, y_end],
                                height=height,
                                width=width,
                                match_aspect=True,
                                tools=['reset,box_zoom,pan,save'],
                                toolbar_location="left")
            return hm

        def plot_gridlines(plot, substrates, enzymes):
            h_lines, v_lines = self.get_lines(substrates, enzymes)
            for line in h_lines:
                plot.line(x=line[0], y=line[1], line_color='darkgrey', alpha=0.7, line_width=0.5)
            for line in v_lines:
                plot.line(x=line[0], y=line[1], line_color='darkgrey', alpha=0.7, line_width=0.5)

        def plot_enzyme_dendogram(plot, ordered_enzymes, substrates, icoord, dcoord, enz_colours):
            enz_colour_pallet = list(palettes.d3['Category10'][5])
            enz_colour_map = {'C0': 'black'}
            for c in enz_colours:
                if (c != 'C0') and (c not in enz_colour_map):
                    enz_colour_map[c] = enz_colour_pallet.pop(0)
                if len(enz_colour_pallet) == 0:
                    enz_colour_pallet = list(palettes.d3['Category10'][5])

            x_width = len(ordered_enzymes) - 1 * self.block_x_size
            enz_icoord, enz_dcoord = self.process_enzyme_dendogram_coords(icoord, dcoord, x_width, substrates)
            for i, d, c in zip(enz_icoord, enz_dcoord, enz_colours):
                d = list(map(lambda x: x, d))
                line_colour = enz_colour_map[c]
                plot.line(x=i, y=d, line_color=line_colour, line_width=1)

        def plot_substrate_dendogram(plot, substrates, icoord, dcoord):
            colour_map = {'C0': 'black'}
            colour_pallet = list(palettes.d3['Category10'][5])
            for c in colours:
                if (c != 'C0') and (c not in colour_map):
                    colour_map[c] = colour_pallet.pop(0)
                if len(colour_pallet) == 0:
                    colour_pallet = list(palettes.d3['Category10'][5])

            y_height = len(substrates) - 1 * self.block_y_size
            icoord, dcoord = self.process_substrate_dendogram_coords(icoord, dcoord, y_height)
            for i, d, c in zip(icoord, dcoord, colours):
                d = list(map(lambda x: -x, d))
                line_colour = colour_map[c]
                plot.line(x=d, y=i, line_color=line_colour, line_width=1)

        def plot_boxes(plot, source):
            b = plot.rect(x='xs', y='ys',
                    height=self.block_y_size,
                    width=self.block_x_size,
                    fill_color='colors',
                    line_color='black',
                    source=source,
                    line_alpha=0.2,
                    fill_alpha=self.data_alpha,
                    hover_alpha=self.hover_alpha,
                    name='blocks'
                    )

            b.selection_glyph = b.glyph
            b.nonselection_glyph = b.glyph

            HT = HoverTool(names=['blocks'])
            HT.tooltips = [('', '@info{safe}')]
            HT.mode = 'mouse'
            plot.add_tools(HT)

            callback = CustomJS(args=dict(source=source),
                                code=f"""
                                    var id = source.selected.indices[0]
                                    var act_ids = source.data.ids[id]
                                    heatmap_launch_box_modal(act_ids)
                                    """)

            #plot.js_on_event(Tap, callback)

            TT = TapTool(names=['blocks'])
            TT.callback = callback
            plot.add_tools(TT)

        def plot_enzyme_names(plot, source):
            text_glyph = Text(x="x", y="y", text="text",
                              angle=1.57, text_font_size='10px', text_align='left',
                              text_font=value("Courier"))  # text_color='colour'
            plot.add_glyph(source, text_glyph)

        def plot_mol_svgs(plot, source):
            image1 = ImageURL(url='url', x='x', y='y', anchor='center_left', h=self.mol_height, w=self.mol_width)
            plot.add_glyph(source, image1)

        def final_plot_modifications(plot):
            plot.axis.major_tick_line_color = None
            plot.axis.minor_tick_line_color = None
            plot.axis.major_label_text_color = None
            plot.axis.major_label_text_font_size = '0pt'
            plot.axis.axis_line_color = None
            plot.grid.grid_line_color = None
            plot.outline_line_color = None
            plot.toolbar.active_drag = None

            if self.show_toolbar == False:
                plot.toolbar.logo = None
                # hm.toolbar_location = None

        unique_smiles = get_data.get_unique_smiles(self.activity_data, self.smi_col)
        fps = get_data.get_substrates_fps(unique_smiles)
        enzymes = get_data.get_unique_enzymes(self.activity_data, self.enzyme_col)
        embeddings = get_data.get_enzyme_embeddings(enzymes)

        substrates, icoord, dcoord, colours = self.substrate_hierarchical_clustering(unique_smiles, fps)
        ordered_enzymes, other_enzymes, enz_icoord, enz_dcoord, enz_colours, enz_leaf_colours = self.enzyme_hierarchical_clustering(enzymes, embeddings)

        enzymes = ordered_enzymes + other_enzymes

        svg_map, svg_url_map = get_data.get_mol_svgs(substrates)
        data = self.get_plot_data(substrates, enzymes, svg_map)
        data_svg = self.get_svg_data(substrates, svg_url_map)
        data_titles = self.get_plot_enzyme_names_data(substrates, enzymes, ordered_enzymes)

        plot = create_figure(substrates, enzymes)
        plot_gridlines(plot, substrates, enzymes)
        plot_enzyme_dendogram(plot, ordered_enzymes, substrates, enz_icoord, enz_dcoord, enz_colours)
        plot_substrate_dendogram(plot, substrates, icoord, dcoord)
        plot_boxes(plot, ColumnDataSource(data))
        plot_enzyme_names(plot, ColumnDataSource(data_titles))
        plot_mol_svgs(plot, ColumnDataSource(data_svg))
        final_plot_modifications(plot)

        return plot

    def log(self, msg, level=1):
        if self.log_level >= level:
            print(f"HeatMapper ({level}): {msg}")




if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    dq = get_data.DataQuery(enzyme_type='CAR', reaction='Carboxylic acid reduction', log_level=1)
    hm = HeatMapper(dq, log_level=1, plot_mode='categorical')
    heat_map = hm.create_heatmap()

    from bokeh.plotting import figure, output_file, save
    output_file('heatmap_test.html')
    save(heat_map)

    #from bokeh.embed import components
    #script, div = components(heat_map)