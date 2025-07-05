🧪 Enhanced Molecule Analyzer ✨




 <!-- IMPORTANT: Replace with your actual Streamlit app URL -->

📝 Table of Contents
🚀 About The Project

✨ Features

🛠️ Built With

💻 Getting Started

Prerequisites

Installation

💡 Usage

🗺️ Roadmap

🤝 Contributing

⭐ Show Your Support

📄 License

📞 Contact

🚀 About The Project
The Enhanced Molecule Analyzer is a powerful Streamlit web application designed for advanced chemical compound analysis. It provides a comprehensive suite of tools for searching, visualizing, and analyzing molecular properties, including 3D structures, drug-likeness, and potential toxicity warnings. Whether you're a chemist, a student, or a researcher, this tool aims to streamline your molecular exploration by offering an intuitive interface and robust analytical capabilities.

✨ Features
🔍 Intelligent Molecule Search: Find molecules by name, SMILES, InChIKey, or formula using robust API integrations (PubChem, CACTUS).

🧬 Interactive 2D & 3D Structure Viewer: Visualize molecular structures in both 2D and interactive 3D, powered by RDKit and 3Dmol.js.

📊 Comprehensive Property Analysis:

Detailed display of molecular weight, LogP, TPSA, hydrogen bond counts, rotatable bonds, and more.

💊 Drug-likeness Scoring: Evaluate compounds based on Lipinski's Rule of Five and a custom drug-likeness score (0-100).

⚠️ Toxicity Warnings: Automated flagging of potential toxicity concerns based on structural alerts and property thresholds.

📈 Advanced Analytics & Visualization:

Radar Charts: Visualize drug-likeness profiles.

Compliance Bars: See Lipinski's Rule of Five compliance at a glance.

Correlation Matrices & Distributions: Analyze relationships and distributions of properties across multiple molecules.

Comparison Tables: Easily compare properties of loaded molecules.

📋 Batch Processing: Upload CSV, Excel, or TXT files for bulk analysis of multiple molecules, saving valuable time.

💾 Data Export: Download your analysis results as CSV or Excel files for further use.

🎨 Sleek Dark Theme: A modern, eye-pleasing dark interface for comfortable viewing.

🛠️ Built With
This project is built using the following technologies and libraries:

Streamlit - For building interactive web applications in Python.

Python

Pandas - For data manipulation and analysis.

Plotly - For interactive data visualizations.

RDKit - For cheminformatics functionalities (molecular properties, 2D/3D structures).

Requests - For making HTTP requests to external APIs.

NumPy - For numerical operations.

SciPy - For scientific computing (e.g., statistics).

Scikit-learn - For machine learning utilities (e.g., clustering, scaling).

PubChem PUG REST API - For chemical information.

CACTUS Chemical Structure Lookup Service - For chemical name to structure conversions.

3Dmol.js - For interactive 3D molecular visualization.

💻 Getting Started
To get a local copy up and running, follow these simple steps.

Prerequisites
Ensure you have Python 3.8+ installed on your system.

Installation
Clone the repository:

git clone https://github.com/shubhae/your-repo.git
cd your-repo

(Remember to replace your-repo with your actual GitHub repository name.)

Create a virtual environment (recommended):

python -m venv venv
source venv/bin/activate  # On Windows: `venv\Scripts\activate`

Install dependencies:

pip install -r requirements.txt

(Make sure you have a requirements.txt file in your project root with all the necessary dependencies listed.)

Run the Streamlit application:

streamlit run streamlit.py

The application will open in your default web browser at http://localhost:8501.

💡 Usage
Once the application is running, you can:

Search & Analyze: Use the "Search & Analysis" tab to input molecule names, SMILES strings, or InChIKeys. The app will fetch detailed information, calculate properties, and display drug-likeness and toxicity warnings.

View Structures: Navigate to the "Structure Viewer" tab to see interactive 2D and 3D representations of your selected molecule.

Explore Analytics: The "Analytics" tab provides statistical summaries, property distributions, and correlation charts for single or multiple loaded molecules.

Batch Process: Use the "Batch Processing" tab to upload a file (CSV, Excel, TXT) containing multiple molecule identifiers for bulk analysis.

Export Data: In the "Export" tab, you can download all the analyzed molecule data as CSV or Excel files.

🗺️ Roadmap
[ ] Integrate more external chemical databases (e.g., ChEMBL, ZINC).

[ ] Implement QSAR/QSPR model predictions (e.g., solubility, binding affinity).

[ ] Add substructure search and similarity search capabilities.

[ ] Enhance 3D visualization with more rendering options and molecular dynamics.

[ ] Develop a user authentication system for saving and loading sessions.

[ ] Improve error handling and user feedback for API calls.

🤝 Contributing
Contributions are what make the open-source community such an amazing place to learn, inspire, and create. Any contributions you make are greatly appreciated.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

Fork the Project

Create your Feature Branch (git checkout -b feature/AmazingFeature)

Commit your Changes (git commit -m 'Add some AmazingFeature')

Push to the Branch (git push origin feature/AmazingFeature)

Open a Pull Request

⭐ Show Your Support
If you find this project useful, please consider giving it a star! Your support helps us grow and improve.

📄 License
Distributed under the MIT License. See LICENSE for more information.

📞 Contact
Shubhae - shubhamytedt@gmail.com

Project Link: https://github.com/shubhae/Molecule_Analyser <!-- IMPORTANT: Replace with your actual GitHub repo URL -->
