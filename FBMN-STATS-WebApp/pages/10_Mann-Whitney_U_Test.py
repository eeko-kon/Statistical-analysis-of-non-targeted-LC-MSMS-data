import streamlit as st

from src.common import *
from src.mwu import *

page_setup()

st.markdown("# Mann-Whitney U Test")

with st.expander("📖 About"):
    st.markdown(
        "The Mann-Whitney U test, also known as the Wilcoxon rank sum test, is a non-parametric statistical test used to compare two ***independent*** groups. It's particularly useful when the data is not normally distributed or when the assumptions of parametric tests like the independent t-test are violated. Essentially, it assesses whether two groups are likely drawn from the same population, meaning they have the same distribution."
    )

if not st.session_state.data.empty:
    c1, c2 = st.columns(2)
    c1.selectbox(
        "select attribute of interest",
        options=[c for c in st.session_state.md.columns if len(set(st.session_state.md[c])) > 1],
        key="mwu_attribute",
    )
    attribute_options = list(
        set(st.session_state.md[st.session_state.mwu_attribute].dropna())
    )
    attribute_options.sort()
    c2.multiselect(
        "select **two** options from the attribute for comparison",
        options=attribute_options,
        default=attribute_options[:2],
        key="mwu_options",
        max_selections=2,
        help="Select two options.",
    )
    c1.selectbox("alternative", options=["two-sided", "greater", "less"], key="mwu_alternative", help="Defines the alternative hypothesis, or tail of the test.")


    if c2.button("Run MWU", type="primary", disabled=(len(st.session_state.mwu_options) != 2)):
        st.session_state.df_mwu = gen_mwu_data(
            st.session_state.mwu_attribute,
            st.session_state.mwu_options,
            st.session_state.mwu_alternative,
            corrections_map[st.session_state.p_value_correction]
        )
        st.rerun()

    if not st.session_state.df_mwu.empty:
        tabs = st.tabs(
            ["📈 Feature significance", "📊 Single metabolite plots", "📁 Data"]
        )
        with tabs[0]:
            fig = get_mwu_plot(st.session_state.df_mwu)
            show_fig(fig, "Mann-Whitney U")
        with tabs[1]:
            cols = st.columns(2)
            cols[0].selectbox(
                "metabolite", st.session_state.df_mwu.index, key="mwu_metabolite"
            )
            fig = mwu_boxplot(st.session_state.df_mwu,
                st.session_state.mwu_metabolite
            )
            show_fig(fig, f"mwu-boxplot-{st.session_state.mwu_metabolite}", False)

        with tabs[2]:
            show_table(st.session_state.df_mwu, "mwu-data")
