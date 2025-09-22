# CHEETAH Designer

A Python tool for designing **euRptr** and **peuRptr** constructs based on sequence motifs.  
This project was developed for RNA recording using genome editing tools (Cas9, ABE, PE).

---

## Features
- **euRptr candidate design**
  - Two modes:
    - **Cas9/PE**: standard design
    - **ABE**: requires user-defined editing window (default: 1–10 bp); only designs if an **A** is present in that window
  - Detects motifs of type **VHHNNN**
  - Extracts upstream 20 nt for target sequence
  - Generates modified antirepeat according to defined rules:
    1. Reverse complement of downstream 20 nt
    2. Replace the **last (20th) base** with `T`
    3. Remove the **7th–8th bases from the end** and insert `AAGT`

- **peuRptr design**
  - Assemble sequences in the order:  
    `antirepeat + scaffold + RTT + PBS + (optional linker) + tevopreQ1`
  - Outputs results as a tab-separated table
  - First row contains the **full joined peuRptr sequence**

---

## Usage
Run the program:
```bash
python CHEETAH_euRptr_designer.py
```

You will be prompted to select:
1. **Design mode**  
   - `1`: euRptr candidate design (Cas9/PE/ABE)  
   - `2`: peuRptr design  

2. **Inputs**  
   - For euRptr design: DNA sequence, optional gene name, optional output path  
   - For ABE mode: also enter the editing window (e.g., `1-10`)  
   - For peuRptr design: antirepeat, RTT, PBS (required), linker (optional)

3. **Output**  
   - Results are saved as a `.txt` file (tab-delimited) in the chosen folder  
   - If no path is given, the current working directory is used

---

## Output Example

### euRptr candidate design (Cas9/PE)
```
euRptr Name   Position   euRptr antirepeat              Target sequence
example_euRptr1  21         TCCCAAAAGTTGCCACGCTAGC        ATCGTTGCGGGCTAGCGTGG
```

### peuRptr design
```
peuRptr sequence:   CCTTGAGTAACTaaggctagtccgtt...cgcggttctatctagttacgcgttaaaccaactagaa
antirepeat   scaffold   RTT   PBS   linker   tevopreQ1
CCTTGAGTAACT aaggct...  AGCCTTGGATT   TTGAGC   GGCC   cgcggt...
```
