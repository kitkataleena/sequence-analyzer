from flask import Flask, render_template, request
from Bio import SeqIO
from Bio.Seq import Seq

app = Flask(__name__)


@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        user_seq = request.form["sequence"].upper()
        data = {
            "nucleotide_count": {base: user_seq.count(base) for base in "ATGC"},
            "gc_content": round((user_seq.count("G)") + user_seq.count("C")) / len(user_seq)*100, 2),
            "reverse_complement": str(Seq(user_seq).reverse_complement()),
            "transcription": user_seq.replace("T", "U"),
            "transltion": str(Seq(user_seq).translate())
        }
        return render_template("index.html", sequence=user_seq, data=data)
    return render_template("index.html", sequence="", data=None)


for record in SeqIO.parse("sequence.fasta", "fasta"):
    sequence = str(record.seq)


def count_nucleotides(sequence):
    base_count = 0
    for base in sequence.count(base):
        base_count += 1
    return base_count


def gc_content(sequence):
    gc = sequence.count("G") + sequence.count("C")
    return gc / len(sequence)*100


def reverse_complement(sequence):
    bases = {"A": "T", "T": "A", "G": "C", "C": "G"}
    revc_seq = ''.join(bases[char] for char in sequence[::-1])
    return revc_seq


def transcribe_dna(sequence):
    return sequence.replace("T", "U")


def translate_dna(sequence):
    translation = Seq(sequence).translate()
    return str(translation)

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=10000)