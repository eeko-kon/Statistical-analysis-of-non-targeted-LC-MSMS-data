import streamlit as st
import pandas as pd
import pingouin as pg
import plotly.express as px
import numpy as np


@st.cache_data
def gen_mwu_data(mwu_attribute, target_groups, alternative, p_correction):
    df = pd.concat([st.session_state.data, st.session_state.md], axis=1)
    mwu = []
    for col in st.session_state.data.columns:
        group1 = df[col][df[mwu_attribute] == target_groups[0]]
        group2 = df[col][df[mwu_attribute] == target_groups[1]]
        result = pg.mwu(group1, group2, alternative, method="exact")
        result["metabolite"] = col
        mwu.append(result)

    mwu = pd.concat(mwu).set_index("metabolite")
    mwu = mwu.dropna()

    mwu.insert(5, "p-corrected", pg.multicomp(mwu["p-val"].astype(float), method=p_correction)[1])
    # add significance
    mwu.insert(6, "significance", mwu["p-corrected"] < 0.05)
    mwu.insert(7, "mwu_attribute", mwu_attribute)
    mwu.insert(8, "A", target_groups[0])
    mwu.insert(9, "B", target_groups[1])

    return mwu.sort_values("p-corrected")

@st.cache_resource
def get_mwu_plot(mwu):
    # first plot insignificant features
    fig = px.scatter(
        x=mwu[mwu["significance"] == False]["U-val"].apply(np.log),
        y=mwu[mwu["significance"] == False]["p-corrected"].apply(
            lambda x: -np.log(x)),
        template="plotly_white",
        width=600,
        height=600,
    )
    fig.update_traces(marker_color="#696880")

    # plot significant features
    fig.add_scatter(
        x=mwu[mwu["significance"]]["U-val"].apply(np.log),
        y=mwu[mwu["significance"]]["p-corrected"].apply(lambda x: -np.log(x)),
        mode="markers+text",
        text=mwu.index[:6],
        textposition="top left",
        textfont=dict(color="#ef553b", size=14),
        name="significant",
    )

    fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={
            "text": f"mwu - {st.session_state.mwu_attribute.upper()}",
            "font_color": "#3E3D53"
        },
        xaxis_title="log(U-val)",
        yaxis_title="-log(p-val)",
        showlegend=False
    )
    fig.update_yaxes(title_standoff=10)
  
    # fig.update_yaxes(title_font_size=20)
    # fig.update_xaxes(title_font_size=20)

    return fig

@st.cache_resource
def mwu_boxplot(df_mwu, metabolite):
    df = pd.concat([st.session_state.md, st.session_state.data], axis=1)
    df1 = pd.DataFrame(
        {
            metabolite: df[df[st.session_state.mwu_attribute] == st.session_state.mwu_options[0]].loc[:, metabolite],
            "option": st.session_state.mwu_options[0],
        }
    )
    df2 = pd.DataFrame(
        {
            metabolite: df[df[st.session_state.mwu_attribute] == st.session_state.mwu_options[1]].loc[:, metabolite],
            "option": st.session_state.mwu_options[1],
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
        xaxis_title=st.session_state.mwu_attribute.replace("st.session_state.mwu_attribute_", ""),
        yaxis_title="intensity",
        template="plotly_white",
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={
            "text": metabolite,
            "font_color": "#3E3D53",
        },
    )
    fig.update_yaxes(title_standoff=10)
    pvalue = df_mwu.loc[metabolite, "p-corrected"]
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
        x0, x1 = st.session_state.mwu_options[0], st.session_state.mwu_options[1]
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