from Bio import SeqIO
from rich.console import Console
from rich.text import Text
from Bio.Align import PairwiseAligner


def search_and_store_indices(query, fasta_file_path):
    results = {}

    # Parse the FASTA file and look for the query
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        sequence_str = str(record.seq)
        matches = {}

        # Find all occurrences of the query and store indices and characters
        start = 0
        while True:
            start = sequence_str.find(query, start)
            if start == -1:
                break
            for i in range(start, start + len(query)):
                matches[i] = sequence_str[i]  # Map index to character
            start += len(query)  # Move past this match

        # Add matches to results if any were found
        if matches:
            results[record.id] = {"description": record.description, "matches": matches}

    return results


# Example usage
search_query = input("Enter the sequence query: ").upper()
fasta_file_path = "sequences.fasta"
results = search_and_store_indices(search_query, fasta_file_path)

# Print the results
"""
for seq_id, data in results.items():
    print(f"Sequence ID: {seq_id}")
    print(f"Description: {data['description']}")
    print("Matches:")
    for index, char in data["matches"].items():
        print(f"  Index: {index}, Character: {char}")
    print()
"""


# (Naive) function for sequence comparison
from Bio.Align import PairwiseAligner


def compare_regions(region1, region2, name):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    alignments = aligner.align(region1, region2)

    print(f"\n{name} Comparison:")
    for alignment in alignments:
        print(alignment)
        print(f"Score: {alignment.score}")
        break  # Show only the best alignment


# Define regions
antibody1 = {
    "CDR1": "DYY",
    "FR2": "IHWVRQRPEQGLEWIGW",
    "CDR2": "LDPENGDTESAP",
    "FR3": "KFQGKATMTADTSSNTAYLQLSSLTSEASAVYYC",
    "CDR3": "NAISTTRDYYALDY",
    "FR4": "WGQGTSVTVSS",
}

antibody2 = {
    "CDR1": "NYW",
    "FR2": "IDWIKQRPGHGLEWIGE",
    "CDR2": "ILPGSGSTNYNEKF",
    "FR3": "RGKATFTADTSSNTAYMQLSSLTSEDSAVYY",
    "CDR3": "TRRGYWAYDFDY",
    "FR4": "WGQGTTLTVSS",
}

# Compare regions
for region in antibody1.keys():
    compare_regions(antibody1[region], antibody2[region], region)

# (BLAST) function for sequence comparison
from Bio.Blast.Applications import NcbiblastnCommandline

import subprocess


def create_blast_database(input_fasta, db_name, db_type="nucl"):
    """
    Creates a BLAST database from the input FASTA file.
    """
    command = [
        "makeblastdb",  # BLAST+ command
        "-in",
        input_fasta,  # Input FASTA file
        "-dbtype",
        db_type,  # Database type: 'nucl' for nucleotide, 'prot' for protein
        "-out",
        db_name,  # Name of the output database
    ]

    try:
        subprocess.run(command, check=True)
        print(f"BLAST database '{db_name}' created successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error while creating BLAST database: {e}")


# Define the BLAST command
blastn_cline = NcbiblastnCommandline(
    query="query.fasta",  # Input query file
    db="local_db",  # Local BLAST database
    evalue=0.001,  # E-value threshold
    outfmt=5,  # Output format (XML)
    out="local_blast_results.xml",  # Output file
)

# Run BLAST
stdout, stderr = blastn_cline()
print(stdout)
print(stderr)

# Parse the results
from Bio.Blast import NCBIXML

with open("local_blast_results.xml") as result_file:
    blast_records = NCBIXML.parse(result_file)
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                print(f"****Alignment****")
                print(f"Sequence: {alignment.title}")
                print(f"Length: {alignment.length}")
                print(f"E-value: {hsp.expect}")
                print(hsp.query)
                print(hsp.match)
                print(hsp.sbjct)


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
