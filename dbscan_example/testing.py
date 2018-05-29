import plotly
from plotly.graph_objs import  Layout,Scatter3d

plotly.offline.plot({
    "data": [Scatter3d(x=[1, 2, 3, 4], y=[4, 3, 2, 1],z=[4, 3, 2, 1])],
    "layout": Layout(title="hello world")
})