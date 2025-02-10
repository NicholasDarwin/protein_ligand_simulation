import pandas as pd
import streamlit as st

def export_analysis(data):
    df = pd.DataFrame(data)
    st.write(df)
    df.to_csv("analysis_results.csv")
    st.download_button("Download Results", "analysis_results.csv")
