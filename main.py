from flask import Flask, render_template, request, redirect, url_for, flash, send_file
import os, io, sys
import zipfile
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline

app = Flask(__name__)
app.secret_key = 'O&*^Y*Liuhr3f78ihfed'  # Required for flashing messages

# Add the path to the BLAST executables to sys.path
sys.path.append(os.path.join(os.getcwd(),"ncbi-blast-2.16.0+/bin"))


# Path to the blastn executable
BLASTN_EXECUTABLE = os.path.join(os.getcwd(),"ncbi-blast-2.16.0+/bin/blastn")

FASTA_OUTPUT_DIR = os.path.join("fasta_output")
BLAST_OUTPUT_DIR = os.path.join("blast_output")

QUALITY_THRESHOLD = int(30)


UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

ALLOWED_EXTENSIONS = {'zip'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def upload_form():
    return render_template('upload.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        flash('No file part')
        return redirect(url_for('upload_form'))

    file = request.files['file']
    threshold = int(request.form['qscore'])

    if file.filename == '':
        flash('No selected file')
        return redirect(url_for('upload_form'))

    if file and allowed_file(file.filename):
        #DO SOMETHING WITH THE FILE

        zip_fie = process_sequence_file(threshold, file)
        #Default behavior was to save the file
        flash('File uploaded successfully')
        return send_file(zip_fie, 
                     download_name='files.zip', 
                     as_attachment=True,
                     mimetype='application/zip')
        
    
    flash('Invalid file type. Only ZIP files are allowed.')
    return redirect(request.url)

def process_sequence_file(threshold, file):
    if os.path.exists('fasta_output'):
        shutil.rmtree('fasta_output')
    os.mkdir('fasta_output')
    if os.path.exists('blast_output'):
        shutil.rmtree('blast_output')
    os.mkdir('blast_output')

    zip_data = io.BytesIO(file.read())
    unzipped_files = {}
    with zipfile.ZipFile(zip_data) as zip_file:
        for file_info in zip_file.infolist():
            with zip_file.open(file_info) as extracted_file:
                unzipped_files[file_info.filename] = extracted_file.read()  # Read file content
            print(file_info)

    fasta_files = []
    # Process each .ab1 file
    for file in unzipped_files:
        if file.endswith(".ab1"):
            f = io.BytesIO(unzipped_files[file])
            ab1_path = file
            print(f"Processing file: {ab1_path}")
            
            # Read the .ab1 file
            try:
                record = SeqIO.read(f, "abi")
            except Exception as e:
                print(f"Error reading {ab1_path}: {e}")
                continue
            # Trim low-quality ends
            qual_scores = record.letter_annotations["phred_quality"]
            
            # Determine trimming positions
            start, end = 0, len(qual_scores)
            while start < end and qual_scores[start] < threshold:
                start += 1
            while end > start and qual_scores[end - 1] < threshold:
                end -= 1
            
            # Trim sequence and quality scores
            trimmed_seq = record.seq[start:end]
            
            # Create a FASTA record
            fasta_record = SeqRecord(trimmed_seq, id=record.id, description="Trimmed AB1 sequence")

            # Save as FASTA in the single output folder
            fasta_filename = os.path.splitext(os.path.basename(file))[0] + ".fasta"
            fasta_path = os.path.join(FASTA_OUTPUT_DIR, fasta_filename)
            
            fasta_filename = os.path.splitext(os.path.basename(file))[0] + ".fasta"
            fasta_path = os.path.join(FASTA_OUTPUT_DIR, fasta_filename)
            SeqIO.write(fasta_record, fasta_path, "fasta")
            fasta_files.append(fasta_filename)

            print(f"Processed and saved: {fasta_filename}")

    print(f"Fasta files converted.")

    # Convert reference.txt to FASTA format and save in fasta_output
    reference_txt_path = "./reference.txt"
    reference_fasta_path = "./reference.fasta"

    # Check if reference.txt exists
    if not os.path.exists(reference_txt_path):
        print(f"Error: {reference_txt_path} does not exist. Please ensure the file is present.")
        exit()

    with open(reference_txt_path, 'r') as ref_txt, open(reference_fasta_path, 'w') as ref_fasta:
        ref_fasta.write(">reference_sequence\n")
        for line in ref_txt:
            ref_fasta.write(line.strip() + "\n")

    print(f"Converted reference.txt to reference.fasta and saved in: {FASTA_OUTPUT_DIR}")

    # Perform BLAST alignment for each FASTA file against reference.fasta
    for fasta_file in os.listdir(FASTA_OUTPUT_DIR):
        if fasta_file.endswith(".fasta") and fasta_file != "reference.fasta":
            fasta_path = os.path.join(FASTA_OUTPUT_DIR, fasta_file)
            th = os.path.join(FASTA_OUTPUT_DIR, fasta_file)
            blast_output = os.path.join(BLAST_OUTPUT_DIR, os.path.splitext(fasta_file)[0] + "_blast.txt")
            
            # Check if the query sequence is not empty
            query_seq = SeqIO.read(fasta_path, "fasta").seq
            if len(query_seq) == 0:
                print(f"Warning: {fasta_file} contains no data. Skipping BLAST alignment.")
                continue

            blastn_cline = NcbiblastnCommandline(cmd=BLASTN_EXECUTABLE, query=fasta_path, subject=reference_fasta_path, outfmt=0, out=blast_output)
            print(f"Running BLAST for: {fasta_file}")
            stdout, stderr = blastn_cline()
            print(f"BLAST stdout: {stdout}")
            print(f"BLAST stderr: {stderr}")
            print(f"BLAST alignment complete for: {fasta_file}")

    print(f"BLAST processing complete! All BLAST results are in: {BLAST_OUTPUT_DIR}")

    files = []
    zip_buffer = io.BytesIO()

    # Create a ZipFile object in write mode using the BytesIO buffer
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
        # Iterate through all files in the specified directory
        for foldername, subfolders, filenames in os.walk(BLAST_OUTPUT_DIR):
            for filename in filenames:
                # Create the full path to the file
                file_path = os.path.join(foldername, filename)
                
                # Add the file to the zip file, preserving its relative path
                zipf.write(file_path, os.path.relpath(file_path, BLAST_OUTPUT_DIR))
    
    # Move the buffer's cursor to the beginning so it can be read later
    zip_buffer.seek(0)
    return zip_buffer



if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port='8013')