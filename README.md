<!-- Enhanced Molecule Analyzer README -->

<div align="center">
  <img src="https://img.shields.io/badge/Streamlit-App-blueviolet?style=for-the-badge&logo=streamlit" />
  <img src="https://img.shields.io/badge/RDKit-Chemoinformatics-green?style=for-the-badge" />
  <img src="https://img.shields.io/badge/PubChem-API-success?style=for-the-badge&logo=pubchem" />
</div>

<h1 align="center">🧪 Enhanced Molecule Analyzer</h1>

<p align="center">Advanced chemical compound analysis with beautiful dark mode UI, 2D/3D structure visualization, drug-likeness scores, and batch processing.</p>

---

## 🌐 Live App (Optional)
<!-- Uncomment if deployed -->
<!-- 🔗 [Try It on Streamlit Cloud](https://share.streamlit.io/your-username/enhanced-molecule-analyzer/main/streamlit.py) -->

## 🚀 Features

- 🔍 **Search** any molecule by name, formula, InChIKey, or SMILES
- 🧬 **2D and 3D structure viewer** using RDKit + py3Dmol
- 📊 **Molecular property charts**, radar plots, and Lipinski compliance
- 💊 **Drug-likeness scoring** based on rules of five
- ⚠️ **Toxicity warnings** and structural flags
- 📦 **Batch processing** via CSV/XLSX or text input
- 📁 **Export** analyzed data to CSV or Excel
- 🎨 **Fully custom dark UI** with gradient cards and visual elements

---

## 🖼️ Screenshots

| Search & Analysis | 3D Visualization | Drug Score |
|------------------|------------------|-------------|
| ![search](https://github.com/shubhae/PokeSQL/assets/placeholder1.png) | ![3d](https://github.com/shubhae/PokeSQL/assets/placeholder2.png) | ![drug](https://github.com/shubhae/PokeSQL/assets/placeholder3.png) |

---

## 🛠️ Tech Stack

- `Python 3.9+`
- `Streamlit`
- `RDKit`
- `Plotly`
- `PubChem PUG REST API`
- `CACTUS Chemical Identifier Resolver`
- `py3Dmol` (for 3D structures)

---

## 📦 Installation

```bash
# Clone the repo
git clone https://github.com/your-username/enhanced-molecule-analyzer.git
cd enhanced-molecule-analyzer

# (Recommended) Create virtual env
python -m venv venv
source venv/bin/activate   # On Windows: venv\Scripts\activate

# Install requirements
pip install -r requirements.txt

# Run the app
streamlit run streamlit.py
```

---

## 🧠 How It Works

- Uses **PubChem** and **CACTUS** APIs to fetch molecule data  
- Uses **RDKit** for parsing SMILES and calculating chemical descriptors  
- Calculates **drug-likeness** using **Lipinski's Rule of Five**  
- Visualizes chemical structures and properties using **Plotly**, **py3Dmol**, and **Streamlit components**

---

## 📂 Project Structure

```plaintext
.
├── streamlit.py              # Main app
├── requirements.txt          # Dependencies
├── assets/                   # (Optional) Custom logos, screenshots
└── README.md                 # You're here!
```

---

## ✅ Example Inputs

```text
aspirin
C9H8O4
BSYNRYMUTXBXSQ-UHFFFAOYSA-N
O=C(C)Oc1ccccc1C(=O)O
```

---

## 📑 License

MIT License © [shubhae](https://github.com/shubhae)

---

## 🌟 Star if useful!

If this project helped you in chemoinformatics or bioinformatics workflows, consider giving it a ⭐️ on GitHub!
