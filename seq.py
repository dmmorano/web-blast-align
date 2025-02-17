import os
import zipfile
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline

# Add the path to the BLAST executables to sys.path
sys.path.append(os.path.join(os.getcwd,"ncbi-blast-2.16.0+/bin"))

# Set the directory to /Users/justinmorano/Desktop/Seq
SCRIPT_DIR = "/Users/justinmorano/Desktop/Seq"

# Identify the first .zip file in the directory
zip_files = [f for f in os.listdir(SCRIPT_DIR) if f.endswith(".zip")]

if not zip_files:
    print("No .zip file found in the script directory. Please add a .zip file and try again.")
    exit()

ZIP_FILE = os.path.join(SCRIPT_DIR, zip_files[0])  # Use the first found .zip file

# Define extraction and output directories
EXTRACT_DIR = os.path.join(SCRIPT_DIR, "extracted_ab1")
FASTA_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "fasta_output")
BLAST_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "blast_output")

# Get quality score from user
QUALITY_THRESHOLD = int(input("Enter the quality score threshold: "))

# Ensure directories exist
os.makedirs(EXTRACT_DIR, exist_ok=True)
os.makedirs(FASTA_OUTPUT_DIR, exist_ok=True)
os.makedirs(BLAST_OUTPUT_DIR, exist_ok=True)

# Extract .ab1 files from the found ZIP file
print(f"Extracting files from: {ZIP_FILE}")
with zipfile.ZipFile(ZIP_FILE, 'r') as zip_ref:
    zip_ref.extractall(EXTRACT_DIR)

# List all files in the extracted directory
extracted_files = []
for root, dirs, files in os.walk(EXTRACT_DIR):
    for file in files:
        extracted_files.append(os.path.join(root, file))

print(f"Extracted files: {extracted_files}")

# Process each .ab1 file
for file in extracted_files:
    if file.endswith(".ab1"):
        ab1_path = file
        print(f"Processing file: {ab1_path}")
        
        # Read the .ab1 file
        try:
            record = SeqIO.read(ab1_path, "abi")
        except Exception as e:
            print(f"Error reading {ab1_path}: {e}")
            continue
        
        # Trim low-quality ends
        qual_scores = record.letter_annotations["phred_quality"]
        
        # Determine trimming positions
        start, end = 0, len(qual_scores)
        while start < end and qual_scores[start] < QUALITY_THRESHOLD:
            start += 1
        while end > start and qual_scores[end - 1] < QUALITY_THRESHOLD:
            end -= 1
        
        # Trim sequence and quality scores
        trimmed_seq = record.seq[start:end]
        
        # Create a FASTA record
        fasta_record = SeqRecord(trimmed_seq, id=record.id, description="Trimmed AB1 sequence")

        # Save as FASTA in the single output folder
        fasta_filename = os.path.splitext(os.path.basename(file))[0] + ".fasta"
        fasta_path = os.path.join(FASTA_OUTPUT_DIR, fasta_filename)
        SeqIO.write(fasta_record, fasta_path, "fasta")

        print(f"Processed and saved: {fasta_filename}")

print(f"Processing complete! All FASTA files are in: {FASTA_OUTPUT_DIR}")

# Convert reference.txt to FASTA format and save in fasta_output
reference_txt_path = os.path.join(SCRIPT_DIR, "reference.txt")
reference_fasta_path = os.path.join(FASTA_OUTPUT_DIR, "reference.fasta")

# Check if reference.txt exists
if not os.path.exists(reference_txt_path):
    print(f"Error: {reference_txt_path} does not exist. Please ensure the file is present.")
    exit()

with open(reference_txt_path, 'r') as ref_txt, open(reference_fasta_path, 'w') as ref_fasta:
    ref_fasta.write(">reference_sequence\n")
    for line in ref_txt:
        ref_fasta.write(line.strip() + "\n")

print(f"Converted reference.txt to reference.fasta and saved in: {FASTA_OUTPUT_DIR}")

# Path to the blastn executable
BLASTN_EXECUTABLE = "/Users/justinmorano/Downloads/ncbi-blast-2.16.0+/bin/blastn"

# Perform BLAST alignment for each FASTA file against reference.fasta
for fasta_file in os.listdir(FASTA_OUTPUT_DIR):
    if fasta_file.endswith(".fasta") and fasta_file != "reference.fasta":
        fasta_path = os.path.join(FASTA_OUTPUT_DIR, fasta_file)
        blast_output = os.path.join(BLAST_OUTPUT_DIR, os.path.splitext(fasta_file)[0] + "_blast.txt")
        
        blastn_cline = NcbiblastnCommandline(cmd=BLASTN_EXECUTABLE, query=fasta_path, subject=reference_fasta_path, outfmt=0, out=blast_output)
        print(f"Running BLAST for: {fasta_file}")
        stdout, stderr = blastn_cline()
        print(f"BLAST stdout: {stdout}")
        print(f"BLAST stderr: {stderr}")
        print(f"BLAST alignment complete for: {fasta_file}")

print(f"BLAST processing complete! All BLAST results are in: {BLAST_OUTPUT_DIR}")
