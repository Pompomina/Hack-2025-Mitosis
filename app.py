from flask import Flask, json, render_template, request, send_from_directory, jsonify
import os
from antibodies import search_and_store_indices

app = Flask(__name__)


@app.route("/", methods=["GET"])
def get_antibodies():

    antibodies_query = request.args.get("query", "").lower()
    all_antibodies = []

    try:
        with open("sequences.json", "r") as file:
            all_antibodies = json.load(file)
    except FileNotFoundError:
        app.logger.error("sequences.json file not found.")
        return "Error: sequences.json file not found.", 500

    app.logger.info(f"Loaded {len(all_antibodies)} antibodies.")

    highlights = []
    if len(antibodies_query) > 0:
        matches = search_and_store_indices(
            antibodies_query, "sequences.fasta"
        )  # fasta file path hardcoded for demo
        matches_list = []
        for match in matches:
            id, chain, type, genus, _ = match["description"].split("|")
            sequence = match["sequence"]
            matches_list.append(
                {
                    "id": id,
                    "chain": chain,
                    "type": type,
                    "genus": genus,
                    "sequence": sequence,
                }
            )
        for match in matches:
            highlights.append(match["matches"])
        antibodies_matching_query = matches_list
    else:
        antibodies_matching_query = all_antibodies

    return render_template(
        "antibodies.html",
        antibodies_list=antibodies_matching_query,
        antibodies_query=antibodies_query,
        highlights=highlights,
    )


@app.route("/compound-image", methods=["GET"], endpoint="get_compound_image")
def get_compound_image():

    smiles = request.args.get("smiles", "")

    if not smiles:
        return jsonify({"error": "Missing smiles parameter"}), 400

    try:
        id_value = smiles.split("=")[1]
    except IndexError:
        app.logger.error(f"Invalid smiles parameter: {smiles}")
        return jsonify({"error": "Invalid smiles parameter"}), 400

    image_dir = "static"
    image_file = f"{id_value}.jpeg"

    try:
        return send_from_directory(image_dir, image_file, mimetype="image/jpeg")
    except FileNotFoundError:
        app.logger.error(f"Image not found: {image_file}")
        return jsonify({"error": f"Image not found for ID {id_value}"}), 404


if __name__ == "__main__":
    app.run(debug=True)
