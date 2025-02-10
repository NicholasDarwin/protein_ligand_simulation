import streamlit as st

def ligand_screening():
    st.title("Ligand Screening")
    st.text("Upload ligand library and start screening.")
    file = st.file_uploader("Upload Ligand Library", type=["sdf", "mol2"])
    if file:
        st.write("Ligand library uploaded.")
