import pandas as pd
from sklearn.decomposition import PCA

from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool
from bokeh.models.mappers import CategoricalColorMapper
from bokeh.palettes import Category20

from retrobiocat_web.analysis.enzyme_polar_plot.mol_drawing import get_mol_svgs

def pca(high_d):
    pca = PCA(n_components=2)
    crds = pca.fit_transform(high_d)
    pca_df = pd.DataFrame(crds, columns=["x", "y"])
    return pca_df

def get_fp_pca(substrates, fp_dict):
    fp_df = pca(list(fp_dict.values()))
    xs, ys = list(fp_df['x']), list(fp_df['y'])
    xy_dict = {}
    for i, smi in enumerate(substrates):
        xy_dict[smi] = (xs[i], ys[i])
    return xy_dict, xs, ys

def make_bokeh_plot(xs, ys, smis, labels):
    from bokeh.plotting import figure
    svgs, svg_urls = get_mol_svgs(smis)

    data = ColumnDataSource(data={'x': xs,
                                  'y': ys,
                                  'smis': smis,
                                  'labels': labels,
                                  'svg': svgs})

    fig = figure(plot_width=600, plot_height=500, tools=['reset,box_zoom,wheel_zoom,zoom_in,zoom_out,pan,tap'],
                 title="PCA of molecules")

    mapper = CategoricalColorMapper(factors=list(set(labels)), palette=Category20[20])

    fig.circle('x', 'y', name='buyable', source=data, color={'field': 'labels', 'transform': mapper})

    HT = HoverTool(names=['buyable'])
    HT.tooltips = [("", "<div>@svg{safe}</div>@smis")]
    HT.mode = 'mouse'
    fig.add_tools(HT)

    return fig