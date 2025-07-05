# Enhanced Molecule Analyzer - Complete Rewrite
# Advanced chemical compound analysis with 3D visualization and reliable search

import streamlit as st
import requests
import pandas as pd
import numpy as np
import json
import io
from datetime import datetime
from typing import Dict, Optional
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import warnings
warnings.filterwarnings('ignore')

# Page Configuration - MUST BE FIRST STREAMLIT COMMAND
st.set_page_config(
    page_title="Enhanced Molecule Analyzer",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Dark Theme CSS
st.markdown("""
<style>
    /* Force dark theme */
    .stApp {
        background-color: #0e1117 !important;
        color: #fafafa !important;
    }
    
    /* Main container */
    .main .block-container {
        background-color: #0e1117 !important;
        color: #fafafa !important;
    }
    
    /* Headers */
    h1, h2, h3, h4, h5, h6 {
        color: #fafafa !important;
    }
    
    /* Text elements */
    p, div, span {
        color: #fafafa !important;
    }
    
    /* Cards and boxes */
    .molecule-card {
        background: linear-gradient(135deg, #1e1e1e 0%, #2d2d2d 100%);
        padding: 1.5rem;
        border-radius: 15px;
        border: 1px solid #404040;
        margin: 1rem 0;
        color: #fafafa !important;
    }
    
    .property-box {
        background: linear-gradient(135deg, #2d2d2d 0%, #404040 100%);
        padding: 1rem;
        border-radius: 10px;
        border-left: 4px solid #667eea;
        margin: 0.5rem 0;
        color: #fafafa !important;
    }
    
    .drug-like {
        background: linear-gradient(135deg, #28a745 0%, #20c997 100%);
        color: white;
        padding: 1rem;
        border-radius: 10px;
        text-align: center;
        margin: 1rem 0;
    }
    
    .not-drug-like {
        background: linear-gradient(135deg, #dc3545 0%, #fd7e14 100%);
        color: white;
        padding: 1rem;
        border-radius: 10px;
        text-align: center;
        margin: 1rem 0;
    }
    
    /* Buttons */
    .stButton > button {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        border-radius: 25px;
        padding: 0.5rem 2rem;
        font-weight: bold;
    }
    
    /* Input fields */
    .stTextInput > div > div > input {
        background-color: #2d2d2d !important;
        color: #fafafa !important;
        border: 1px solid #404040 !important;
    }
    
    /* Select boxes */
    .stSelectbox > div > div > select {
        background-color: #2d2d2d !important;
        color: #fafafa !important;
    }
    
    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        background-color: #1e1e1e !important;
    }
    
    .stTabs [data-baseweb="tab"] {
        color: #fafafa !important;
    }
</style>
""", unsafe_allow_html=True)

# Import RDKit for molecular analysis
RDKIT_AVAILABLE = False
RDKIT_DRAW_AVAILABLE = False
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem, Draw, DataStructs
    
    # Try to import rdMolDraw2D (optional)
    try:
        from rdkit.Chem import rdMolDraw2D
        RDKIT_DRAW_AVAILABLE = True
    except ImportError:
        pass
    
    # Test if RDKit is working
    test_mol = Chem.MolFromSmiles("C")
    if test_mol:
        RDKIT_AVAILABLE = True
except Exception as e:
    st.warning("⚠️ RDKit not available. Some advanced features will be disabled.")

# Import scientific libraries
try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def search_molecule(query: str) -> Optional[Dict]:
    """
    Search for molecule information using multiple APIs
    """
    query = query.strip()
    if not query:
        return None
    
    # Try PubChem first
    try:
        # Get basic info
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/property/InChIKey,Title,MolecularFormula,ExactMass/JSON"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                result = data['PropertyTable']['Properties'][0]
                
                # Get SMILES from CACTUS
                try:
                    smiles_url = f"https://cactus.nci.nih.gov/chemical/structure/{query}/smiles"
                    smiles_response = requests.get(smiles_url, timeout=10)
                    if smiles_response.status_code == 200:
                        # Clean the SMILES string - remove any extra characters
                        smiles_text = smiles_response.text.strip()
                        # Remove any trailing commas, semicolons, or other punctuation
                        smiles_text = smiles_text.rstrip(',;.!?')
                        result['CanonicalSMILES'] = smiles_text
                    else:
                        result['CanonicalSMILES'] = 'N/A'
                except:
                    result['CanonicalSMILES'] = 'N/A'
                
                # Get IUPAC name from CACTUS
                try:
                    iupac_url = f"https://cactus.nci.nih.gov/chemical/structure/{query}/iupac_name"
                    iupac_response = requests.get(iupac_url, timeout=10)
                    if iupac_response.status_code == 200:
                        result['IUPACName'] = iupac_response.text.strip()
                    else:
                        result['IUPACName'] = 'N/A'
                except:
                    result['IUPACName'] = 'N/A'
                
                return result
    except Exception as e:
        st.error(f"PubChem API error: {str(e)}")
    
    # Try CACTUS as fallback
    try:
        cactus_url = f"https://cactus.nci.nih.gov/chemical/structure/{query}/inchikey"
        response = requests.get(cactus_url, timeout=10)
        if response.status_code == 200:
            inchikey = response.text.strip()
            if inchikey:
                return {"InChIKey": inchikey, "Title": query}
    except:
        pass
    
    return None

def get_molecule_details(inchikey: str, smiles: str = None, iupac_name: str = None) -> Dict:
    """
    Get detailed molecule information from PubChem
    """
    # Always use the best available SMILES
    best_smiles = smiles if smiles and smiles != 'N/A' else None
    if best_smiles:
        # Clean the SMILES string - remove any extra characters
        best_smiles = best_smiles.strip().rstrip(',;.!?')
    
    result = {
        "InChIKey": inchikey,
        "Name": "Unknown",
        "Formula": "N/A",
        "Molecular_Weight": 0.0,
        "SMILES": best_smiles if best_smiles else "N/A",
        "IUPAC_Name": iupac_name if iupac_name else "N/A",
        "LogP": 0.0,
        "TPSA": 0.0,
        "NumRotatableBonds": 0,
        "NumHDonors": 0,
        "NumHAcceptors": 0,
        "Drug_Likeness_Score": 0,
        "ToxicityWarnings": []
    }
    
    try:
        # Get basic properties
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/MolecularFormula,ExactMass/JSON"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                props = data["PropertyTable"]["Properties"][0]
                result.update({
                    "Formula": props.get("MolecularFormula", "N/A"),
                    "Molecular_Weight": float(props.get("ExactMass", 0.0))
                })
                
                # Get compound name
                try:
                    name_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/synonyms/JSON"
                    name_response = requests.get(name_url, timeout=10)
                    if name_response.status_code == 200:
                        name_data = name_response.json()
                        if "InformationList" in name_data and "Information" in name_data["InformationList"]:
                            synonyms = name_data["InformationList"]["Information"][0].get("Synonym", [])
                            if synonyms:
                                result["Name"] = synonyms[0]
                except:
                    pass
                
                # Calculate RDKit properties
                if RDKIT_AVAILABLE and result["SMILES"] != "N/A":
                    mol = Chem.MolFromSmiles(result["SMILES"])
                    if mol:
                        result.update({
                            "LogP": round(Descriptors.MolLogP(mol), 2),
                            "TPSA": round(Descriptors.TPSA(mol), 2),
                            "NumRotatableBonds": int(Descriptors.NumRotatableBonds(mol)),
                            "NumHDonors": int(Descriptors.NumHDonors(mol)),
                            "NumHAcceptors": int(Descriptors.NumHAcceptors(mol)),
                            "HeavyAtomCount": int(Descriptors.HeavyAtomCount(mol))
                        })
                        
                        # Calculate aromatic rings and ring count manually
                        try:
                            aromatic_rings = Descriptors.NumAromaticRings(mol)
                            result["NumAromaticRings"] = int(aromatic_rings)
                        except:
                            result["NumAromaticRings"] = "N/A"
                        
                        try:
                            ring_count = Descriptors.RingCount(mol)
                            result["RingCount"] = int(ring_count)
                        except:
                            result["RingCount"] = "N/A"
                        
                        # Calculate drug-likeness score
                        score = calculate_drug_likeness(result)
                        result["Drug_Likeness_Score"] = score
                        
                        # Check for toxicity warnings
                        warnings = check_toxicity_warnings(mol)
                        result["ToxicityWarnings"] = warnings
                    else:
                        pass
                else:
                    pass
    except Exception as e:
        st.error(f"Error fetching molecule details: {str(e)}")
    
    return result

def calculate_drug_likeness(mol_data: Dict) -> int:
    """
    Calculate drug-likeness score (0-100)
    """
    score = 0
    mw = mol_data.get("Molecular_Weight", 0)
    logp = mol_data.get("LogP", 0)
    hbd = mol_data.get("NumHDonors", 0)
    hba = mol_data.get("NumHAcceptors", 0)
    rotatable = mol_data.get("NumRotatableBonds", 0)
    
    # Lipinski's Rule of Five
    if 160 <= mw <= 500: score += 20
    if -0.4 <= logp <= 5.6: score += 20
    if hbd <= 5: score += 20
    if hba <= 10: score += 20
    if rotatable <= 10: score += 20
    
    return score

def check_toxicity_warnings(mol) -> list:
    """
    Check for potential toxicity warnings
    """
    warnings = []
    
    if not mol:
        return warnings
    
    # Check molecular weight
    mw = Descriptors.ExactMolWt(mol)
    if mw > 500:
        warnings.append("High molecular weight (>500 Da) - may have poor bioavailability")
    
    # Check LogP
    logp = Descriptors.MolLogP(mol)
    if logp > 5:
        warnings.append("High lipophilicity (LogP >5) - may cause toxicity")
    elif logp < -0.4:
        warnings.append("Low lipophilicity (LogP <-0.4) - may have poor absorption")
    
    # Check hydrogen bond donors/acceptors
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    
    if hbd > 5:
        warnings.append("Too many hydrogen bond donors (>5)")
    if hba > 10:
        warnings.append("Too many hydrogen bond acceptors (>10)")
    
    # Check for toxic substructures
    toxic_patterns = [
        ("[N+](=O)[O-]", "Contains nitro group"),
        ("[S-]", "Contains thiol group"),
        ("[CH2]=O", "Contains aldehyde group"),
        ("C(=O)Cl", "Contains acid chloride"),
        ("[N-]=[N+]=[N-]", "Contains azide group")
    ]
    
    for pattern, warning in toxic_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            warnings.append(warning)
    
    return warnings

def create_3d_visualization(smiles: str) -> Optional[str]:
    """
    Create 3D molecular visualization using py3Dmol
    """
    if not RDKIT_AVAILABLE or not smiles or smiles == "N/A":
        st.error("No valid SMILES for 3D structure.")
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            st.error("RDKit could not parse SMILES for 3D structure.")
            return None

        mol = Chem.AddHs(mol)
        # Try standard embedding
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception as e1:
            # Try ETKDG as fallback
            try:
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception as e2:
                st.error(f"3D structure generation failed: {e1} / {e2}")
                return None

        pdb = Chem.MolToPDBBlock(mol)
        if not pdb or len(pdb.strip()) == 0:
            st.error("3D structure could not be generated (empty PDB).")
            return None

        viewer_html = f"""
        <script src="https://3dmol.org/build/3Dmol-min.js"></script>
        <style>
            .mol-container {{
                width: 100%;
                height: 500px;
                background: #1e1e1e;
                border-radius: 10px;
                margin: 1rem 0;
            }}
        </style>
        <div id="mol-viewer" class="mol-container"></div>
        <script>
            let viewer = $3Dmol.createViewer(document.getElementById("mol-viewer"), {{
                backgroundColor: 0x1e1e1e
            }});
            viewer.addModel(`{pdb}`, "pdb");
            viewer.setStyle({{}}, {{
                stick: {{}},
                sphere: {{radius: 0.5}}
            }});
            viewer.zoomTo();
            viewer.render();
        </script>
        """
        return viewer_html
    except Exception as e:
        st.error(f"3D visualization error: {str(e)}")
        return None

def create_2d_structure(smiles: str) -> Optional[str]:
    """
    Create 2D molecular structure
    """
    if not RDKIT_AVAILABLE or not smiles or smiles == "N/A":
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Generate SVG
        svg = Draw.MolToSVG(mol, width=400, height=300)
        
        # Style the SVG for dark theme - make it more visible
        svg = svg.replace('<svg', '<svg style="background-color: #1e1e1e; border-radius: 10px; display: block; margin: auto;"')
        
        # Add some basic styling to make sure it's visible
        svg = f"""
        <div style="text-align: center; padding: 20px; background-color: #1e1e1e; border-radius: 10px;">
            {svg}
        </div>
        """
        
        return svg
    except Exception as e:
        return None

def create_property_charts(mol_data: Dict) -> Dict:
    """
    Create comprehensive property visualization charts
    """
    charts = {}
    
    # Drug-likeness radar chart
    properties = ['Molecular Weight', 'LogP', 'H-Bond Donors', 'H-Bond Acceptors', 'Rotatable Bonds']
    values = [
        min(mol_data.get('Molecular_Weight', 0) / 500, 1) * 100,  # Normalize to 0-100
        min(max(mol_data.get('LogP', 0) + 1, 0) / 6, 1) * 100,    # Normalize -1 to 5
        min(mol_data.get('NumHDonors', 0) / 5, 1) * 100,          # Normalize to 0-5
        min(mol_data.get('NumHAcceptors', 0) / 10, 1) * 100,      # Normalize to 0-10
        min(mol_data.get('NumRotatableBonds', 0) / 10, 1) * 100   # Normalize to 0-10
    ]
    
    fig_radar = go.Figure()
    fig_radar.add_trace(go.Scatterpolar(
        r=values,
        theta=properties,
        fill='toself',
        name='Drug-likeness',
        line_color='#667eea'
    ))
    
    fig_radar.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 100]
            )),
        showlegend=False,
        title="Drug-likeness Profile",
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        font=dict(color='#fafafa')
    )
    
    charts['radar'] = fig_radar
    
    # Lipinski's Rule of Five compliance
    lipinski_rules = {
        'Molecular Weight ≤ 500': mol_data.get('Molecular_Weight', 0) <= 500,
        'LogP ≤ 5': mol_data.get('LogP', 0) <= 5,
        'H-Bond Donors ≤ 5': mol_data.get('NumHDonors', 0) <= 5,
        'H-Bond Acceptors ≤ 10': mol_data.get('NumHAcceptors', 0) <= 10
    }
    
    fig_lipinski = go.Figure(data=[
        go.Bar(
            x=list(lipinski_rules.keys()),
            y=[1 if rule else 0 for rule in lipinski_rules.values()],
            marker_color=['#28a745' if rule else '#dc3545' for rule in lipinski_rules.values()],
            text=[f"{'✓' if rule else '✗'}" for rule in lipinski_rules.values()],
            textposition='auto'
        )
    ])
    
    fig_lipinski.update_layout(
        title="Lipinski's Rule of Five Compliance",
        yaxis=dict(range=[0, 1], ticktext=['Fail', 'Pass'], tickvals=[0, 1]),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0.1)',
        font=dict(color='#fafafa')
    )
    
    charts['lipinski'] = fig_lipinski
    
    # Molecular property distribution (if we have multiple molecules)
    return charts

def create_advanced_analytics(molecules: list) -> Dict:
    """
    Create advanced analytics for multiple molecules
    """
    charts = {}
    
    if len(molecules) < 2:
        return charts
    
    df = pd.DataFrame(molecules)
    
    # Property correlation matrix
    numeric_cols = ['Molecular_Weight', 'LogP', 'TPSA', 'NumRotatableBonds', 'NumHDonors', 'NumHAcceptors']
    available_cols = [col for col in numeric_cols if col in df.columns and df[col].notna().any()]
    
    if len(available_cols) >= 2:
        corr_matrix = df[available_cols].corr()
        
        fig_corr = go.Figure(data=go.Heatmap(
            z=corr_matrix.values,
            x=corr_matrix.columns,
            y=corr_matrix.columns,
            colorscale='RdBu',
            zmid=0
        ))
        
        fig_corr.update_layout(
            title="Property Correlation Matrix",
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0.1)',
            font=dict(color='#fafafa')
        )
        
        charts['correlation'] = fig_corr
    
    # Property distributions
    if 'Molecular_Weight' in df.columns:
        fig_dist = make_subplots(
            rows=2, cols=2,
            subplot_titles=("Molecular Weight", "LogP", "TPSA", "Rotatable Bonds"),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        # Molecular Weight
        fig_dist.add_trace(
            go.Histogram(x=df['Molecular_Weight'].dropna(), name="MW", nbinsx=10),
            row=1, col=1
        )
        
        # LogP
        if 'LogP' in df.columns:
            fig_dist.add_trace(
                go.Histogram(x=df['LogP'].dropna(), name="LogP", nbinsx=10),
                row=1, col=2
            )
        
        # TPSA
        if 'TPSA' in df.columns:
            fig_dist.add_trace(
                go.Histogram(x=df['TPSA'].dropna(), name="TPSA", nbinsx=10),
                row=2, col=1
            )
        
        # Rotatable Bonds
        if 'NumRotatableBonds' in df.columns:
            fig_dist.add_trace(
                go.Histogram(x=df['NumRotatableBonds'].dropna(), name="RotBonds", nbinsx=10),
                row=2, col=2
            )
        
        fig_dist.update_layout(
            height=600,
            showlegend=False,
            title_text="Property Distributions",
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0.1)',
            font=dict(color='#fafafa')
        )
        
        charts['distributions'] = fig_dist
    
    # Scatter plots
    if all(col in df.columns for col in ['Molecular_Weight', 'LogP']):
        fig_scatter = px.scatter(
            df, x='Molecular_Weight', y='LogP',
            title="Molecular Weight vs LogP",
            hover_data=['Name'],
            color='Drug_Likeness_Score' if 'Drug_Likeness_Score' in df.columns else None
        )
        fig_scatter.update_layout(
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0.1)',
            font=dict(color='#fafafa')
        )
        charts['scatter'] = fig_scatter
    
    # Drug-likeness score distribution
    if 'Drug_Likeness_Score' in df.columns:
        fig_score = go.Figure(data=[
            go.Histogram(
                x=df['Drug_Likeness_Score'].dropna(),
                nbinsx=10,
                marker_color='#667eea'
            )
        ])
        fig_score.update_layout(
            title="Drug-likeness Score Distribution",
            xaxis_title="Score",
            yaxis_title="Count",
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0.1)',
            font=dict(color='#fafafa')
        )
        charts['score_dist'] = fig_score
    
    return charts

def calculate_statistical_summary(molecules: list) -> Dict:
    """
    Calculate comprehensive statistical summary
    """
    if not molecules:
        return {}
    
    df = pd.DataFrame(molecules)
    summary = {}
    
    # Basic statistics
    numeric_cols = ['Molecular_Weight', 'LogP', 'TPSA', 'NumRotatableBonds', 'NumHDonors', 'NumHAcceptors']
    available_cols = [col for col in numeric_cols if col in df.columns]
    
    if available_cols:
        stats_df = df[available_cols].describe()
        summary['basic_stats'] = stats_df
    
    # Drug-likeness analysis
    if 'Drug_Likeness_Score' in df.columns:
        scores = df['Drug_Likeness_Score'].dropna()
        if len(scores) > 0:
            summary['drug_likeness'] = {
                'mean': scores.mean(),
                'median': scores.median(),
                'std': scores.std(),
                'min': scores.min(),
                'max': scores.max(),
                'excellent': len(scores[scores >= 80]),
                'good': len(scores[(scores >= 60) & (scores < 80)]),
                'poor': len(scores[scores < 60])
            }
    
    # Lipinski compliance
    if all(col in df.columns for col in ['Molecular_Weight', 'LogP', 'NumHDonors', 'NumHAcceptors']):
        lipinski_compliance = {
            'mw_ok': len(df[df['Molecular_Weight'] <= 500]),
            'logp_ok': len(df[df['LogP'] <= 5]),
            'hbd_ok': len(df[df['NumHDonors'] <= 5]),
            'hba_ok': len(df[df['NumHAcceptors'] <= 10]),
            'all_rules': len(df[(df['Molecular_Weight'] <= 500) & (df['LogP'] <= 5) & 
                               (df['NumHDonors'] <= 5) & (df['NumHAcceptors'] <= 10)])
        }
        summary['lipinski'] = lipinski_compliance
    
    return summary

# ============================================================================
# MAIN APPLICATION
# ============================================================================

def main():
    # Initialize session state
    if 'molecules' not in st.session_state:
        st.session_state.molecules = []
    
    # Main header
    st.markdown("""
    <div style="text-align: center; padding: 2rem; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); border-radius: 15px; margin-bottom: 2rem;">
        <h1>🧪 Enhanced Molecule Analyzer</h1>
        <p style="font-size: 1.2rem; margin: 0;">Advanced Chemical Compound Analysis with 3D Visualization</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Create tabs
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "🔍 Search & Analysis",
        "🧬 Structure Viewer",
        "📊 Analytics",
        "📋 Batch Processing",
        "💾 Export"
    ])
    
    with tab1:
        st.header("🔍 Search & Analysis")
        
        # Search input
        query = st.text_input(
            "Enter molecule name, SMILES, or InChIKey:",
            placeholder="e.g., aspirin, caffeine, C9H8O4"
        )
        search_button = st.button("🔍 Search", type="primary", use_container_width=True)
        
        if search_button and query:
            with st.spinner("Searching for molecule..."):
                # Search for molecule
                result = search_molecule(query)
                
                if result and result.get('InChIKey'):
                    # Get detailed information
                    mol_data = get_molecule_details(
                        result['InChIKey'],
                        smiles=result.get('CanonicalSMILES'),
                        iupac_name=result.get('IUPACName')
                    )
                    
                    # Store in session state with proper name and IUPAC
                    mol_data['Title'] = result.get('Title', mol_data.get('Name', 'Unknown'))
                    # Ensure IUPAC name is stored correctly
                    if result.get('IUPACName') and result.get('IUPACName') != 'N/A':
                        mol_data['IUPAC_Name'] = result.get('IUPACName')
                    if mol_data not in st.session_state.molecules:
                        st.session_state.molecules.append(mol_data)
                    
                    # Display results
                    st.success(f"✅ Found: {result.get('Title', mol_data.get('Name', 'Unknown'))}")
                    
                    # Basic information card
                    st.markdown("""
                    <div class="molecule-card">
                        <h3>🧪 Basic Information</h3>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown(f"""
                        <div class="property-box">
                            <p><strong>Name:</strong> {result.get('Title', mol_data.get('Name', 'Unknown'))}</p>
                            <p><strong>Formula:</strong> {result.get('MolecularFormula', mol_data.get('Formula', 'N/A'))}</p>
                            <p><strong>IUPAC Name:</strong> {mol_data.get('IUPAC_Name', result.get('IUPACName', 'N/A'))}</p>
                            <p><strong>SMILES:</strong> {mol_data.get('SMILES', result.get('CanonicalSMILES', 'N/A'))}</p>
                        </div>
                        """, unsafe_allow_html=True)
                    
                    with col2:
                        st.markdown(f"""
                        <div class="property-box">
                            <p><strong>Molecular Weight:</strong> {result.get('ExactMass', mol_data.get('Molecular_Weight', 'N/A'))} g/mol</p>
                            <p><strong>LogP:</strong> {mol_data.get('LogP', 'N/A')}</p>
                            <p><strong>TPSA:</strong> {mol_data.get('TPSA', 'N/A')} Å²</p>
                            <p><strong>InChIKey:</strong> {result.get('InChIKey', mol_data.get('InChIKey', 'N/A'))}</p>
                        </div>
                        """, unsafe_allow_html=True)
                    
                    # Molecular properties
                    st.markdown("""
                    <div class="molecule-card">
                        <h3>📊 Molecular Properties</h3>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.markdown(f"""
                        <div class="property-box">
                            <p><strong>H-Bond Donors:</strong> {mol_data.get('NumHDonors', 'N/A')}</p>
                            <p><strong>H-Bond Acceptors:</strong> {mol_data.get('NumHAcceptors', 'N/A')}</p>
                        </div>
                        """, unsafe_allow_html=True)
                    
                    with col2:
                        st.markdown(f"""
                        <div class="property-box">
                            <p><strong>Rotatable Bonds:</strong> {mol_data.get('NumRotatableBonds', 'N/A')}</p>
                            <p><strong>Aromatic Rings:</strong> {mol_data.get('NumAromaticRings', 'N/A')}</p>
                        </div>
                        """, unsafe_allow_html=True)
                    
                    with col3:
                        st.markdown(f"""
                        <div class="property-box">
                            <p><strong>Ring Count:</strong> {mol_data.get('RingCount', 'N/A')}</p>
                            <p><strong>Heavy Atoms:</strong> {mol_data.get('HeavyAtomCount', 'N/A')}</p>
                        </div>
                        """, unsafe_allow_html=True)
                    
                    # Drug-likeness analysis
                    score = mol_data.get('Drug_Likeness_Score', 0)
                    if score > 0:
                        st.markdown("""
                        <div class="molecule-card">
                            <h3>💊 Drug-likeness Analysis</h3>
                        </div>
                        """, unsafe_allow_html=True)
                        
                        if score >= 80:
                            st.markdown(f"""
                            <div class="drug-like">
                                <h4>Excellent Drug-likeness Score: {score}/100</h4>
                                <p>This compound has excellent drug-like properties!</p>
                            </div>
                            """, unsafe_allow_html=True)
                        elif score >= 60:
                            st.markdown(f"""
                            <div class="drug-like">
                                <h4>Good Drug-likeness Score: {score}/100</h4>
                                <p>This compound has good drug-like properties.</p>
                            </div>
                            """, unsafe_allow_html=True)
                        else:
                            st.markdown(f"""
                            <div class="not-drug-like">
                                <h4>Drug-likeness Score: {score}/100</h4>
                                <p>This compound may have unfavorable drug-like properties.</p>
                            </div>
                            """, unsafe_allow_html=True)
                    
                    # Toxicity warnings
                    warnings = mol_data.get('ToxicityWarnings', [])
                    if warnings:
                        st.markdown("""
                        <div class="molecule-card">
                            <h3>⚠️ Toxicity Warnings</h3>
                        </div>
                        """, unsafe_allow_html=True)
                        
                        for warning in warnings:
                            st.warning(warning)
                
                else:
                    st.error("❌ Could not find the molecule. Please try a different search term.")
    
    with tab2:
        st.header("🧬 Structure Viewer")
        
        if not st.session_state.molecules:
            st.info("🔍 Search for a molecule first to view its structure.")
        else:
            # Select molecule
            mol_names = []
            for i, m in enumerate(st.session_state.molecules):
                name = m.get('Title', m.get('Name', f"Molecule {i+1}"))
                mol_names.append(name)
            
            selected_idx = st.selectbox("Select molecule:", range(len(mol_names)), format_func=lambda i: mol_names[i])
            selected_mol = st.session_state.molecules[selected_idx]
            
            # 2D Structure
            st.subheader("2D Structure")
            smiles = selected_mol.get('SMILES', 'N/A')
            if smiles != 'N/A':
                svg = create_2d_structure(smiles)
                if svg:
                    st.markdown(svg, unsafe_allow_html=True)
                else:
                    st.info("2D structure could not be generated.")
            else:
                st.info("No SMILES available for 2D structure.")
            
            # 3D Structure
            st.subheader("3D Structure")
            if smiles != 'N/A':
                viewer_html = create_3d_visualization(smiles)
                if viewer_html:
                    st.components.v1.html(viewer_html, height=600)
                else:
                    st.info("3D structure could not be generated.")
            else:
                st.info("No SMILES available for 3D structure.")
            
            # Show molecule info
            st.subheader("Molecule Information")
            col1, col2 = st.columns(2)
            with col1:
                st.write(f"**Name:** {selected_mol.get('Name', selected_mol.get('Title', 'Unknown'))}")
                st.write(f"**IUPAC Name:** {selected_mol.get('IUPAC_Name', selected_mol.get('IUPACName', 'N/A'))}")
                st.write(f"**Formula:** {selected_mol.get('Formula', selected_mol.get('MolecularFormula', 'N/A'))}")
                st.write(f"**SMILES:** {smiles}")
            with col2:
                st.write(f"**Molecular Weight:** {selected_mol.get('Molecular_Weight', selected_mol.get('ExactMass', 'N/A'))} g/mol")
                st.write(f"**LogP:** {selected_mol.get('LogP', 'N/A')}")
                st.write(f"**TPSA:** {selected_mol.get('TPSA', 'N/A')} Å²")
    
    with tab3:
        st.header("📊 Advanced Analytics")
        
        if not st.session_state.molecules:
            st.info("🔍 Search for molecules first to view analytics.")
        else:
            # Statistical Summary
            st.subheader("📈 Statistical Summary")
            summary = calculate_statistical_summary(st.session_state.molecules)
            
            if summary:
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    if 'basic_stats' in summary:
                        st.markdown("**📊 Basic Statistics**")
                        st.dataframe(summary['basic_stats'], use_container_width=True)
                
                with col2:
                    if 'drug_likeness' in summary:
                        st.markdown("**💊 Drug-likeness Analysis**")
                        dl = summary['drug_likeness']
                        st.metric("Mean Score", f"{dl['mean']:.1f}")
                        st.metric("Excellent (≥80)", dl['excellent'])
                        st.metric("Good (60-79)", dl['good'])
                        st.metric("Poor (<60)", dl['poor'])
                
                with col3:
                    if 'lipinski' in summary:
                        st.markdown("**🔬 Lipinski Compliance**")
                        lip = summary['lipinski']
                        st.metric("All Rules Pass", lip['all_rules'])
                        st.metric("MW ≤ 500", lip['mw_ok'])
                        st.metric("LogP ≤ 5", lip['logp_ok'])
                        st.metric("HBD ≤ 5", lip['hbd_ok'])
            
            # Single molecule analysis
            if len(st.session_state.molecules) == 1:
                st.subheader("🧪 Single Molecule Analysis")
                mol_data = st.session_state.molecules[0]
                charts = create_property_charts(mol_data)
                
                col1, col2 = st.columns(2)
                with col1:
                    if 'radar' in charts:
                        st.plotly_chart(charts['radar'], use_container_width=True)
                
                with col2:
                    if 'lipinski' in charts:
                        st.plotly_chart(charts['lipinski'], use_container_width=True)
            
            # Multiple molecules analysis
            if len(st.session_state.molecules) > 1:
                st.subheader("🔬 Multiple Molecules Analysis")
                
                # Advanced charts
                advanced_charts = create_advanced_analytics(st.session_state.molecules)
                
                # Correlation matrix
                if 'correlation' in advanced_charts:
                    st.plotly_chart(advanced_charts['correlation'], use_container_width=True)
                
                # Property distributions
                if 'distributions' in advanced_charts:
                    st.plotly_chart(advanced_charts['distributions'], use_container_width=True)
                
                # Scatter plots
                col1, col2 = st.columns(2)
                with col1:
                    if 'scatter' in advanced_charts:
                        st.plotly_chart(advanced_charts['scatter'], use_container_width=True)
                
                with col2:
                    if 'score_dist' in advanced_charts:
                        st.plotly_chart(advanced_charts['score_dist'], use_container_width=True)
                
                # Molecule comparison table
                st.subheader("📋 Molecule Comparison Table")
                df = pd.DataFrame(st.session_state.molecules)
                st.dataframe(df, use_container_width=True)
    
    with tab4:
        st.header("📋 Batch Processing")
        st.markdown("""
        ### Batch Analysis
        Upload a file containing multiple molecule identifiers for bulk analysis:
        - CSV, Excel, or TXT files supported
        - One molecule identifier per line
        - Supports names, InChIKeys, SMILES, and formulas
        """)
        
        uploaded_file = st.file_uploader("Upload file:", type=["csv", "xlsx", "txt"])
        pasted_text = st.text_area("Or paste content here:", height=150, placeholder="Enter one molecule identifier per line (e.g., aspirin, caffeine, paracetamol)")
        
        if uploaded_file is not None or pasted_text:
            if st.button("🚀 Process Batch", type="primary"):
                with st.spinner("Processing molecules..."):
                    molecules_to_process = []
                    
                    if uploaded_file is not None:
                        if uploaded_file.name.endswith('.csv'):
                            df = pd.read_csv(uploaded_file)
                            molecules_to_process = df.iloc[:, 0].tolist()
                        elif uploaded_file.name.endswith('.xlsx'):
                            df = pd.read_excel(uploaded_file)
                            molecules_to_process = df.iloc[:, 0].tolist()
                        else:
                            content = uploaded_file.read().decode()
                            molecules_to_process = [line.strip() for line in content.split('\n') if line.strip()]
                    
                    if pasted_text:
                        molecules_to_process.extend([line.strip() for line in pasted_text.split('\n') if line.strip()])
                    
                    # Remove duplicates
                    molecules_to_process = list(set(molecules_to_process))
                    
                    # Process each molecule
                    successful = 0
                    failed = 0
                    
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    for i, molecule in enumerate(molecules_to_process):
                        status_text.text(f"Processing {molecule}... ({i+1}/{len(molecules_to_process)})")
                        
                        result = search_molecule(molecule)
                        if result and result.get('InChIKey'):
                            mol_data = get_molecule_details(
                                result['InChIKey'],
                                smiles=result.get('CanonicalSMILES'),
                                iupac_name=result.get('IUPACName')
                            )
                            mol_data['Title'] = result.get('Title', mol_data.get('Name', 'Unknown'))
                            if mol_data not in st.session_state.molecules:
                                st.session_state.molecules.append(mol_data)
                            successful += 1
                        else:
                            failed += 1
                        
                        progress_bar.progress((i + 1) / len(molecules_to_process))
                    
                    status_text.text(f"✅ Completed! {successful} successful, {failed} failed")
                    st.success(f"Batch processing complete! Added {successful} molecules to analysis.")
        
        # Quick batch examples
        st.subheader("💡 Quick Examples")
        
        # Show current status
        if st.session_state.molecules:
            st.info(f"📊 Currently loaded: {len(st.session_state.molecules)} molecules")
        else:
            st.info("📊 No molecules loaded yet. Use quick examples or search to add molecules!")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("Add Common Drugs", type="primary"):
                with st.spinner("Adding common drugs..."):
                    common_drugs = ["aspirin", "caffeine", "paracetamol", "ibuprofen", "acetaminophen"]
                    added_count = 0
                    failed_count = 0
                    
                    for drug in common_drugs:
                        result = search_molecule(drug)
                        if result and result.get('InChIKey'):
                            mol_data = get_molecule_details(
                                result['InChIKey'],
                                smiles=result.get('CanonicalSMILES'),
                                iupac_name=result.get('IUPACName')
                            )
                            mol_data['Title'] = result.get('Title', mol_data.get('Name', 'Unknown'))
                            # Ensure IUPAC name is stored correctly
                            if result.get('IUPACName') and result.get('IUPACName') != 'N/A':
                                mol_data['IUPAC_Name'] = result.get('IUPACName')
                            if mol_data not in st.session_state.molecules:
                                st.session_state.molecules.append(mol_data)
                                added_count += 1
                        else:
                            failed_count += 1
                    
                    if added_count > 0:
                        st.success(f"✅ Added {added_count} common drugs! ({failed_count} failed)")
                    else:
                        st.error(f"❌ Failed to add any drugs. {failed_count} failed.")
        
        with col2:
            if st.button("Add Natural Products", type="primary"):
                with st.spinner("Adding natural products..."):
                    natural_products = ["caffeine", "morphine", "cocaine", "nicotine", "quinine"]
                    added_count = 0
                    failed_count = 0
                    
                    for product in natural_products:
                        result = search_molecule(product)
                        if result and result.get('InChIKey'):
                            mol_data = get_molecule_details(
                                result['InChIKey'],
                                smiles=result.get('CanonicalSMILES'),
                                iupac_name=result.get('IUPACName')
                            )
                            mol_data['Title'] = result.get('Title', mol_data.get('Name', 'Unknown'))
                            # Ensure IUPAC name is stored correctly
                            if result.get('IUPACName') and result.get('IUPACName') != 'N/A':
                                mol_data['IUPAC_Name'] = result.get('IUPACName')
                            if mol_data not in st.session_state.molecules:
                                st.session_state.molecules.append(mol_data)
                                added_count += 1
                        else:
                            failed_count += 1
                    
                    if added_count > 0:
                        st.success(f"✅ Added {added_count} natural products! ({failed_count} failed)")
                    else:
                        st.error(f"❌ Failed to add any products. {failed_count} failed.")
        
        with col3:
            if st.button("Add Simple Molecules", type="primary"):
                with st.spinner("Adding simple molecules..."):
                    simple_molecules = ["water", "ethanol", "methane", "glucose", "sucrose"]
                    added_count = 0
                    failed_count = 0
                    
                    for molecule in simple_molecules:
                        result = search_molecule(molecule)
                        if result and result.get('InChIKey'):
                            mol_data = get_molecule_details(
                                result['InChIKey'],
                                smiles=result.get('CanonicalSMILES'),
                                iupac_name=result.get('IUPACName')
                            )
                            mol_data['Title'] = result.get('Title', mol_data.get('Name', 'Unknown'))
                            # Ensure IUPAC name is stored correctly
                            if result.get('IUPACName') and result.get('IUPACName') != 'N/A':
                                mol_data['IUPAC_Name'] = result.get('IUPACName')
                            if mol_data not in st.session_state.molecules:
                                st.session_state.molecules.append(mol_data)
                                added_count += 1
                        else:
                            failed_count += 1
                    
                    if added_count > 0:
                        st.success(f"✅ Added {added_count} simple molecules! ({failed_count} failed)")
                    else:
                        st.error(f"❌ Failed to add any molecules. {failed_count} failed.")
        
        # Clear button in a new row
        col4, col5, col6 = st.columns(3)
        with col5:
            if st.button("🗑️ Clear All Molecules", type="secondary"):
                st.session_state.molecules = []
                st.success("✅ All molecules cleared!")
                st.rerun()
    
    with tab5:
        st.header("💾 Export")
        
        if not st.session_state.molecules:
            st.info("🔍 Search for molecules first to export data.")
        else:
            df = pd.DataFrame(st.session_state.molecules)
            
            # CSV Export
            csv = df.to_csv(index=False)
            st.download_button(
                label="📥 Download CSV",
                data=csv,
                file_name="molecule_analysis.csv",
                mime="text/csv"
            )
            
            # Excel Export
            excel_buffer = io.BytesIO()
            df.to_excel(excel_buffer, index=False, engine='openpyxl')
            excel_data = excel_buffer.getvalue()
            st.download_button(
                label="📥 Download Excel",
                data=excel_data,
                file_name="molecule_analysis.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            
            # Display data table
            st.subheader("Data Preview")
            st.dataframe(df, use_container_width=True)

if __name__ == "__main__":
    main()