🧪 Enhanced Molecule Analyzer ✨
🚀 Overview
The Enhanced Molecule Analyzer is a powerful Streamlit web application designed for advanced chemical compound analysis. It provides a comprehensive suite of tools for searching, visualizing, and analyzing molecular properties, including 3D structures, drug-likeness, and potential toxicity warnings. Whether you're a chemist, a student, or a researcher, this tool aims to streamline your molecular exploration.

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

💻 Installation & Setup
To run this application locally, follow these steps:

Clone the repository:

git clone https://github.com/shubhae/your-repo.git
cd your-repo


(Remember to replace your-repo with your actual GitHub repository details.)

Create a virtual environment (recommended):

python -m venv venv
source venv/bin/activate  # On Windows: `venv\Scripts\activate`


Install dependencies:

pip install -r requirements.txt


(You'll need to create a requirements.txt file. See the "Dependencies" section below.)

Run the Streamlit application:

streamlit run streamlit.py


The application will open in your default web browser.

⚙️ Dependencies
Create a requirements.txt file in your project root with the following content:

streamlit
requests
pandas
numpy
plotly
openpyxl # For Excel export
rdkit-pypi # For RDKit functionality
scipy # For statistical analysis (if used beyond basic pandas)
scikit-learn # For clustering/scaling (if used)


🤝 Contributing
Contributions are welcome! If you have suggestions for improvements, new features, or bug fixes, please feel free to:

Fork the repository.

Create a new branch (git checkout -b feature/YourFeature).

Make your changes.

Commit your changes (git commit -m 'Add some feature').

Push to the branch (git push origin feature/YourFeature).

Open a Pull Request.

⭐ Show Your Support
If you find this project useful, please consider giving it a star! Your support helps us grow and improve.

📄 License
This project is licensed under the MIT License - see the LICENSE file for details.

📞 Contact
For any questions or feedback, please open an issue on GitHub or reach out to shubhamytedt@gmail.com
