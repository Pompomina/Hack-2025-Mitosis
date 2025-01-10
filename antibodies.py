from Bio import SeqIO
from rich.console import Console
from rich.text import Text


def search_and_store_indices(query, fasta_file_path):
    results = {}
    subsections = ["CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]

    # parse the FASTA file and look for the query
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        sequence_str = str(record.seq)
        matches = {}

        # find all occurrences of the query and store indices and characters
        start = 0
        while True:
            start = sequence_str.find(query, start)
            if start == -1:
                break
            for i in range(start, start + len(query)):
                matches[i] = sequence_str[i]  # Map index to character
            start += len(query)  # Move past this match

        # find location of query in antibody
        subsections_matches = {}
        description_parts = record.description.split(" |")  # split description by " |"
        for part in description_parts:
            for subsection in subsections:
                if part.startswith(subsection):  # check if the part is the expected subsection
                    subseq = part.split('=')[1] if '=' in part else ""
                    if query in subseq:
                        subsections_matches[subsection] = subseq
        # Add matches to results if any were found
        if matches:
            results[record.id] = {"description": record.description, "matches": matches, "subsections_matches": subsections_matches }

    return results


# Example usage
search_query = input("Enter the sequence query: ").upper()
fasta_file_path = "sequences.fasta"
results = search_and_store_indices(search_query, fasta_file_path)

# Print the results
for seq_id, data in results.items():
    print(f"Sequence ID: {seq_id}")
    print(f"Description: {data['description']}")
    print("Matches:")
    for index, char in data["matches"].items():
        print(f"  Index: {index}, Character: {char}")
    print()

''' reference code from the demo

def smiles_to_svg(smiles: str, width: int = 400, height: int = 400) -> bytes:
    """
    makes an SVG image of a molecule
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError("Invalid SMILES")

    Chem.rdCoordGen.AddCoords(mol)
    drawer = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    # set drawing options on drawer.getOptions()
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().encode()


def get_substructure_fingerprint(mol):
    """
    TODO: substructure fingerprints
    returns substructure fingerprint for mol
    """
    fp = ExplicitBitVect(SUBSTRUCTURE_FP_SIZE, True)
    return fp  # this is currently empty and useless


def get_highlighted_image(
    target_smiles: str, query_smiles: str, width: int = 400, height: int = 400
):
    """
    TODO: Creates image of highlighted molecule
    """
    target_mol = Chem.MolFromSmiles(target_smiles)
    query_mol = Chem.MolFromSmiles(query_smiles)
    match = target_mol.GetSubstructMatch(query_mol)

    Chem.rdCoordGen.AddCoords(target_mol)
    Chem.rdCoordGen.AddCoords(query_mol)

    drawer = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(target_mol, highlightAtoms=match)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().encode()


def search_compounds(
    query_smiles: str, compound_list: List[str] = EXAMPLE_COMPOUNDS
) -> List[List[int]]:
    """
    search the list of smiles and return substructure match indices

    is empty list if not match and that index

    it would be nice to fingerprint and store fingerprints in memory
    """
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        raise RuntimeError("Invalid query SMILES")

    compounds = [Chem.MolFromSmiles(s) for s in compound_list]
    matches = []
    for m in compounds:
        if m is None:
            matches.append([])
        else:
            matches.append(m.GetSubstructMatch(query_mol))
    return matches

'''
# print(search_compounds("C"))
# print(search_compounds("C1CCCCC1"))

# import rdkit
# print(rdkit.__version__)

# print(smiles_to_svg("CC"))
