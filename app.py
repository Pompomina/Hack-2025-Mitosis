from flask import Flask, json, render_template, request, send_from_directory, jsonify
import os

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

    if antibodies_query:
        antibodies_matching_query = [
            antibody
            for antibody in all_antibodies
            if antibodies_query in antibody["type"].lower()
            or antibodies_query in antibody["id"].lower()
        ]
    else:
        antibodies_matching_query = all_antibodies

    return render_template(
        "antibodies.html",
        antibodies_list=antibodies_matching_query,
        antibodies_query=antibodies_query,
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
