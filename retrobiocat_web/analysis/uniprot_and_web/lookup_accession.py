import requests
from Bio import Entrez, SeqIO

Entrez.email = 'william.finnigan@manchester.ac.uk'

codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

def translate(DNA):
    protein = ''
    for i in range(0,len(DNA), 3):
        protein += codon_table[DNA[i:i+3]]
    return protein


def lookup_uniprot(accession):
    url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
    req = requests.get(url)

    if req.status_code in [200]:
        fasta = req.text
        seq = fasta[fasta.find('\n'):]
        seq = seq.replace('\n', '')
    else:
        seq = ''

    return seq

def lookup_ncbi(accession):
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, 'fasta')
        seq = str(record.seq)
    except:
        seq = ''

    return seq

def lookup_genbank_nucleotide(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, 'fasta')
        dna_seq = str(record.seq)
        seq = translate(dna_seq)

    except:
        seq = ''

    return seq


def get_sequence(accession):
    source = ''
    seq = lookup_uniprot(accession)

    # try uniprot first
    if seq != '':
        source = 'Uniprot'
    else:
        seq = lookup_ncbi(accession)

        # if no sequence, then try NCBI
        if seq != '':
            source = 'NCBI'
        else:
            seq = lookup_genbank_nucleotide(accession)

            # if still no sequence, try genbank
            if seq != '':
                source = 'GenBank'

    return seq, source



if __name__ == '__main__':
    test_accession = "MF540819"
    seq = lookup_genbank_nucleotide(test_accession)
    print(seq)



