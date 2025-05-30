import streamlit as st

from src.common import *
from src.wilcoxon import *

page_setup()

st.markdown("# Wilcoxon Signed-Rank Test")

with st.expander("📖 About"):
    st.markdown(
        "The wilcoxon signed-rank test is a non-parametric statistical test used to compare two ***dependent*** samples (matched pairs) or to compare a single sample against a hypothetical value. It's a non-parametric alternative to the paired t-test and is used when the data is not normally distributed or when the assumptions of the t-test are not met."
    )

if not st.session_state.data.empty:
    c1, c2 = st.columns(2)
    c1.selectbox(
        "Select attribute of interest",
        options=[c for c in st.session_state.md.columns if len(set(st.session_state.md[c])) > 1],
        key="wilcoxon_attribute",
    )
    attribute_options = list(
        set(st.session_state.md[st.session_state.wilcoxon_attribute].dropna())
    )
    attribute_options.sort()
    c2.multiselect(
        "Select **two** options from the attribute for comparison",
        options=attribute_options,
        default=attribute_options[:2],
        key="wilcoxon_options",
        max_selections=2,
        help="Select two options.",
    )
    c1.selectbox("Alternative", options=["two-sided", "greater", "less"], key="wilcoxon_alternative", help="Defines the alternative hypothesis, or tail of the test.")


    if c2.button("Run Wilcoxon", type="primary", disabled=(len(st.session_state.wilcoxon_options) != 2)):
        st.session_state.df_wilcoxon = gen_wilcoxon_data(
            st.session_state.wilcoxon_attribute,
            st.session_state.wilcoxon_options,
            st.session_state.wilcoxon_alternative,
            corrections_map[st.session_state.p_value_correction]
        )
        st.rerun()

    if not st.session_state.df_wilcoxon.empty:
        tabs = st.tabs(
            ["📈 Feature significance", "📊 Single metabolite plots", "📁 Data"]
        )
        with tabs[0]:
            fig = plot_wilcoxon(st.session_state.df_wilcoxon)
            show_fig(fig, "t-test")
        with tabs[1]:
            cols = st.columns(2)
            cols[0].selectbox(
                "metabolite", st.session_state.df_wilcoxon.index, key="wilcoxon_metabolite"
            )
            fig = wilcoxon_boxplot(st.session_state.df_wilcoxon,
                st.session_state.wilcoxon_metabolite
            )
            show_fig(fig, f"wilcoxon-boxplot-{st.session_state.wilcoxon_metabolite}", False)

        with tabs[2]:
            show_table(st.session_state.df_wilcoxon, "wilcoxon-data")
