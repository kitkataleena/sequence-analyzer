from flask import Flask, render_template, request
from Bio.Seq import Seq
from Bio import SeqIO
import os
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = "uploads"
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)


def analyze_sequence(seq):
    user_seq = seq.upper().strip()
    if not user_seq:
        raise ValueError("The sequence is empty.")
    if any(base not in "ATGC" for base in user_seq):
        raise ValueError(
            "The sequence contains invalid characters. Only A, T, C, and G are allowed.")
    if len(user_seq) < 3:
        raise ValueError("The sequence is not long enough to translate")
    return {
        "nucleotide_count": {base: user_seq.count(base) for base in "ATGC"},
        "gc_content": round((user_seq.count("G") + user_seq.count("C")) / len(user_seq)*100, 2),
        "reverse_complement": str(Seq(user_seq).reverse_complement()),
        "transcription": user_seq.replace("T", "U"),
        "translation": str(Seq(user_seq).translate())
    }


@app.route("/", methods=["GET", "POST"])
def index():
    user_seq = ""
    data = None
    error = None
    filename = None

    if request.method == "POST":
        try:
            if request.form.get("sequence"):
                user_seq = request.form["sequence"].upper()
                data = analyze_sequence(user_seq)

            elif "file" in request.files:
                file = request.files["file"]
                if file.filename.endswith(".fasta"):
                    filename = secure_filename(file.filename)
                    filepath = os.path.join(
                        app.config['UPLOAD_FOLDER'], filename)
                    file.save(filepath)

                    with open(filepath, "r") as f:
                        records = list(SeqIO.parse(f, "fasta"))
                        if not records:
                            raise ValueError(
                                "No sequences found in the FASTA file.")
                        user_seq = str(records[0].seq)
                        data = analyze_sequence(user_seq)
                else:
                    raise ValueError("Please upload a valid .fasta file.")
        except ValueError as e:
            error = str(e)

    return render_template("templates.html", sequence=user_seq, data=data, error=error, filename=filename)

if __name__ == "__main__":
    app.run()
