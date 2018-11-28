from os.path import join, dirname

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.cluster import KMeans

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from matplotlib import pyplot as plt

from bokeh.plotting import figure
from bokeh.layouts import layout, column, row
from bokeh.models import ColumnDataSource, Div, HoverTool, TapTool, OpenURL
from bokeh.models.widgets import RangeSlider, Select
from bokeh.io import curdoc

plot_df = pd.read_csv(join(dirname(__file__), 'plot_df.csv'))
selected = plot_df[
        (plot_df.mmgbsa <= -60.0) &
        (plot_df.docking_score <= np.max(plot_df.docking_score.values))
    ]

color_map = {
    "Molecular Fingerprint" : "chem_colors",
    "Binding Site" : "site_colors",
}

desc = Div(text=open(join(dirname(__file__), "description.html")).read(), width=1200)

mmgbsa = RangeSlider(title="MM-GBSA dG Ligand Binding", value=(np.min(plot_df.mmgbsa.values), -60.0),
               start=np.min(plot_df.mmgbsa.values), end=np.max(plot_df.mmgbsa.values),
               step=np.floor((np.max(plot_df.mmgbsa.values)-np.min(plot_df.mmgbsa.values))/1000))
docking_score = RangeSlider(title="Glide XP Docking Score", 
               value=(np.min(plot_df.docking_score.values), np.max(plot_df.docking_score.values)),
               start=np.min(plot_df.docking_score.values), end=np.max(plot_df.docking_score.values),
               step=np.floor((np.max(plot_df.docking_score.values)-np.min(plot_df.docking_score.values))/1000))
set_colors = Select(title="Color Main Plot by Clustering", 
                     options=sorted(color_map.keys()), value="Binding Site")

                       
source = ColumnDataSource(data=dict(x0=[], y0=[], x1=[], y1=[], x2=[], y2=[], color0=[], color1=[], 
                                    color2=[], png=[], url=[], label=[]))

TOOLTIPS=[
        ("label", "@label"),
        ("url", "@url"),
        ("png", "@png"),
        ("x0", "@x0"),
        ("y0", "@y0")
]

hover = HoverTool(tooltips="""
        <div>
                <div>
                        <img
                                src="@png" height="150" alt="@png" width="150"
                                style="float: left; margin: 0px 15px 15px 0px;"
                                border="0"
                        ></img>
                </div>
                <div>@label</div>
                <div>MMGBSA = @x0</div>
                <div>XP = @y0</div>
        </div>
        """
)

p0 = figure(plot_height=500, plot_width=500, title="XP Glide Score vs MM-GBSA dG Ligand Binding",
           tools=['reset, box_zoom, wheel_zoom, zoom_out, pan, save, tap', hover])
p0.circle(x="x0", y="y0", source=source, color="color0", size=6, fill_alpha=0.8)
p0.xaxis.axis_label = 'MM-GBSA dG Ligand Binding'
p0.yaxis.axis_label = 'XP Glide Docking Score'

url = "@url"
taptool0 = p0.select(type=TapTool)
taptool0.callback = OpenURL(url=url)

p1 = figure(plot_height=250, plot_width=500, title="Binding Site PCA",
           tools=['reset, box_zoom, wheel_zoom, zoom_out, pan, save, tap', hover])
p1.circle(x="x1", y="y1", source=source, color="color1", size=6, fill_alpha=0.8)
p1.xaxis.axis_label = 'PC1'
p1.yaxis.axis_label = 'PC2'

taptool1 = p1.select(type=TapTool)
taptool1.callback = OpenURL(url=url)

p2 = figure(plot_height=250, plot_width=500, title="Molecular Fingerprints PCA",
           tools=['reset, box_zoom, wheel_zoom, zoom_out, pan, save, tap', hover])
p2.circle(x="x2", y="y2", source=source, color="color2", size=6, fill_alpha=0.8)
p2.xaxis.axis_label = 'PC1'
p2.yaxis.axis_label = 'PC2'

taptool2 = p2.select(type=TapTool)
taptool2.callback = OpenURL(url=url)

def slider_select():
    global selected
    selected = plot_df[
        (plot_df.mmgbsa >= mmgbsa.value[0]) &
        (plot_df.mmgbsa <= mmgbsa.value[1]) &
        (plot_df.docking_score <= docking_score.value[1]) &
        (plot_df.docking_score >= docking_score.value[0])
    ]
    update()

def update():
    source.data = dict(
        x0 = selected['mmgbsa'],
        y0 = selected['docking_score'],        
        x1 = selected['site_pc1'],
        y1 = selected['site_pc2'],
        x2 = selected['chem_pc1'],
        y2 = selected['chem_pc2'],
        color0 = selected[color_map[set_colors.value]],
        color1 = selected['site_colors'],
        color2 = selected['chem_colors'],
        url = selected['url'],
        png = selected['png'],
        label = selected['label']
    )
    
sliders = [mmgbsa, docking_score]
for slider in sliders:
    slider.on_change('value', lambda attr, old, new: slider_select())

selectors = [set_colors]
for selector in selectors:
    selector.on_change('value', lambda attr, old, new: update())

widgets = column(set_colors, mmgbsa, docking_score)
main = row(widgets, p0, column(p1, p2))

l = layout([
    [desc],
    [main]
], sizing_mode='scale_width')

update()

curdoc().add_root(l)
curdoc().title = 'Docking Summary'