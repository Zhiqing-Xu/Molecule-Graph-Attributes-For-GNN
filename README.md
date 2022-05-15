# Obtain Molecule Graph Attributes For GNN

## Inputs: List of SMILES String

![SMILES_list_phosphatase](https://user-images.githubusercontent.com/47986787/168453834-250e0121-f753-435c-bf25-f381d5696862.png)


## Outputs: 
## Annotaion of molecules
![Molecule_annotation_test](https://user-images.githubusercontent.com/47986787/168453803-13b4c809-3828-418b-9154-95e575d7075b.png)

## Get Morgan Substrctures
![MorganFP_test](https://user-images.githubusercontent.com/47986787/168453822-20f5a344-26b4-428c-89c3-cbd09f12d509.png)

    # 1. Atom_Attributes             size: (n_node, n_attr)               type: numpy.array
    # 2. Atom_RXN_Portion            size: (n_rxn_portion)                type: List
    # 3. Bond_Adjacency_Matrix       size: (n_node, n_node)               type: numpy.array
    # 4. Bond_Attributes             size: (n_node, n_node, dim_attr)     type: numpy.array
