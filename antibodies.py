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
        description_parts = record.description.split("|")  # split description by " |"
        subsections_parts = description_parts[4].split(", ")

#make a dictionary of the item and key
        start = 0
        #section_name_matches = []
        for part in subsections_parts:
                section_name, sequence = part.split("=")
                if query in sequence.upper():
                    subsections_matches[section_name] = sequence

        # Add matches to results if any were found
        if matches:
            results[record.id] = {
                "description": record.description,
                "matches": matches,
                "subsections_matches": subsections_matches,
            }

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
    

    for subsection, subseq in data["subsections_matches"].items():
        print(f"  {subsection}: {subseq}")
    print()
