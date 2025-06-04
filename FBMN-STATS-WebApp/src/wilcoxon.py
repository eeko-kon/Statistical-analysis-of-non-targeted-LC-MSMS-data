import streamlit as st
import pandas as pd
import pingouin as pg
import plotly.express as px
import numpy as np


@st.cache_data
def gen_wilcoxon_data(wilcoxon_attribute, target_groups, alternative, p_correction):
    df = pd.concat([st.session_state.data, st.session_state.md], axis=1)
    wilcoxon = []
    for col in st.session_state.data.columns:
        group1 = df[col][df[wilcoxon_attribute] == target_groups[0]]
        group2 = df[col][df[wilcoxon_attribute] == target_groups[1]]
        if len(group1) != len(group2): #Check if the two groups have equal samples
            raise ValueError("Unequal Sample Sizes")
        result = pg.wilcoxon(group1, group2, alternative)
        result["metabolite"] = col

        wilcoxon.append(result)

    wilcoxon = pd.concat(wilcoxon).set_index("metabolite")
    wilcoxon = wilcoxon.dropna()

    wilcoxon.insert(5, "p-corrected", pg.multicomp(wilcoxon["p-val"].astype(float), method=p_correction)[1])
    # add significance
    wilcoxon.insert(6, "significance", wilcoxon["p-corrected"] < 0.05)
    wilcoxon.insert(7, "wilcoxon_attribute", wilcoxon_attribute)
    wilcoxon.insert(8, "A", target_groups[0])
    wilcoxon.insert(9, "B", target_groups[1])

    return wilcoxon.sort_values("p-corrected")

@st.cache_resource
def get_wilcoxon_plot(wilcoxon):
    # first plot insignificant features
    fig = px.scatter(
        x=wilcoxon[wilcoxon["significance"] == False]["W-val"].apply(np.log),
        y=wilcoxon[wilcoxon["significance"] == False]["p-corrected"].apply(
            lambda x: -np.log(x)),
        template="plotly_white",
        width=600,
        height=600,
    )
    fig.update_traces(marker_color="#696880")

    # plot significant features
    fig.add_scatter(
        x=wilcoxon[wilcoxon["significance"]]["W-val"].apply(np.log),
        y=wilcoxon[wilcoxon["significance"]]["p-corrected"].apply(lambda x: -np.log(x)),
        mode="markers+text",
        text=wilcoxon.index[:6],
        textposition="top left",
        textfont=dict(color="#ef553b", size=14),
        name="significant",
    )

    fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={
            "text": f"wilcoxon - {st.session_state.wilcoxon_attribute.upper()}",
            "font_color": "#3E3D53"
        },
        xaxis_title="log(W-val)",
        yaxis_title="-log(p-val)",
        showlegend=False
    )
    fig.update_yaxes(title_standoff=10)
  
    # fig.update_yaxes(title_font_size=20)
    # fig.update_xaxes(title_font_size=20)

    return fig


@st.cache_resource
def wilcoxon_boxplot(df_wilcoxon, metabolite):
    df = pd.concat([st.session_state.md, st.session_state.data], axis=1)
    df1 = pd.DataFrame(
        {
            metabolite: df[df[st.session_state.wilcoxon_attribute] == st.session_state.wilcoxon_options[0]].loc[:, metabolite],
            "option": st.session_state.wilcoxon_options[0],
        }
    )
    df2 = pd.DataFrame(
        {
            metabolite: df[df[st.session_state.wilcoxon_attribute] == st.session_state.wilcoxon_options[1]].loc[:, metabolite],
            "option": st.session_state.wilcoxon_options[1],
        }
    )
    df = pd.concat([df1, df2])
    fig = px.box(
        df,
        x="option",
        y=metabolite,
        color="option",
        width=350,
        height=400,
        points="all",
    )
    fig.update_layout(
        showlegend=False,
        xaxis_title=st.session_state.wilcoxon_attribute.replace("st.session_state.wilcoxon_attribute_", ""),
        yaxis_title="intensity",
        template="plotly_white",
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={
            "text": metabolite,
            "font_color": "#3E3D53",
        },
    )
    fig.update_yaxes(title_standoff=10)
    pvalue = df_wilcoxon.loc[metabolite, "p-corrected"]
    if pvalue >= 0.05:
        symbol = "ns"
    elif pvalue >= 0.01:
        symbol = "*"
    elif pvalue >= 0.001:
        symbol = "**"
    else:
        symbol = "***"

    top_y = max(df[metabolite]) * 1.2

    if isinstance(df["option"][0], str) and isinstance(df["option"][1], str):
        x0, x1 = 0, 1
    else:
        x0, x1 = st.session_state.wilcoxon_options[0], st.session_state.wilcoxon_options[1]
    
    # horizontal line
    fig.add_shape(
        type="line",
        x0=x0,
        y0=top_y,
        x1=x1,
        y1=top_y,
        line=dict(width=1, color="#000000"),
    )
    if symbol == "ns":
        y_margin = max(df[metabolite]) * 0.05
    else:
        y_margin = max(df[metabolite]) * 0.1
    fig.add_annotation(
        x=(x1-x0)/2,
        y=top_y + y_margin,
        text=f"<b>{symbol}</b>",
        showarrow=False,
        font_color="#555555",
    )
    return fig
