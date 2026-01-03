Link: https://loschmidt.chemi.muni.cz/soluprot/

This tool's scores are useful for predicting expression based on the wild type sequence of amino acids.

About
SoluProt is a web application and standalone software for prediction of soluble protein expression in Escherichia coli.

SoluProt is one of the latest additions to the family of solubility predictors based on machine learning [1]. The training set is based on the TargetTrack database [2], which was carefully filtered to keep only targets expressed in Escherichia coli. The negative and positive samples were balanced and equalized for the protein lengths. The independent validation set is derived from the NESG dataset [3].

The predictor is based on gradient boosting machines and employs 96 sequence-based features, e.g., amino acid content, sequence identity to PDB sequences and several aggregated physico-chemical properties. SoluProt achieves accuracy of 58.5% and AUC of 0.62, higher than other comparable tools.

References
Musil, M., Konegger, H., Hon, J., Bednar, D., Damborsky, J. (2019). Computational Design of Stable and Soluble Biocatalysts. ACS Catalysis 9: 1033−1054.
Berman, H. M., Gabanyi, M. J., Kouranov, A., Micallef, D. I., Westbrook, J., Protein Structure Initiative Network Of Investigators (2017). Protein Structure Initiative – Targettrack 2000–2017.
Price, W. N., Handelman, S. K., Everett, J. K., Tong, S. N., Bracic, A., Luff, J. D., Hunt, J. F. (2011). Large-Scale Experimental Studies Show Unexpected Amino Acid Effects on Protein Expression and Solubility in vivo in E. coli. Microbial Informatics and Experimentation 1: 6.
