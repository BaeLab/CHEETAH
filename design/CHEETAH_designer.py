import os
import re

def reverse_complement(seq):
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq[::-1].translate(table)

def is_nxxnnn(seq):
    if len(seq) != 6:
        return False
    seq = seq.upper()
    # Exclude only motifs with GG in positions 2–3
    return not (seq[1] == "G" and seq[2] == "G")

def generate_antirepeat_corrected(target_seq):
    if len(target_seq) != 20:
        return None
    rc = reverse_complement(target_seq)
    rc = rc[:-1] + 'T'  # last sequence to T
    if len(rc) < 8:
        return None
    cut_pos = -8
    return rc[:cut_pos] + 'AAGT' + rc[cut_pos+2:]

def design_targets(seq, gene, mode="cas9", edit_window=(1, 10)):
    results = [["euRptr Name", "Position", "euRptr antirepeat", "Target sequence"]]
    count = 1
    for i in range(len(seq) - 5):
        motif = seq[i:i+6]
        if not is_nxxnnn(motif):
            continue
        if i < 20 or i + 20 > len(seq):
            continue

        upstream = seq[i - 20:i]
        target_seq = seq[i:i + 20]

        if mode == "abe":
            window_start = edit_window[0] - 1
            window_end = edit_window[1]
            if window_end > len(upstream):
                continue
            window_seq = upstream[window_start:window_end]
            if "A" not in window_seq.upper():
                continue

        antirepeat = generate_antirepeat_corrected(target_seq)
        if not antirepeat:
            continue

        name = f"{gene}_euRptr{count}" if gene else f"euRptr{count}"
        results.append([name, str(i + 1), antirepeat, upstream])
        count += 1
    return results

def save_to_txt(rows, filepath):
    with open(filepath, "w") as f:
        for row in rows:
            f.write("\t".join(row) + "\n")

def parse_editing_window(raw):
    try:
        start, end = map(int, raw.split("-"))
        if start >= 1 and end >= start:
            return (start, end)
    except:
        pass
    return (1, 10)

def run_eurptr_design():
    print("\n[euRptr antirepeat design mode]")
    print("1. Cas9 / PE")
    print("2. ABE")
    mode_input = input("Enter sub-mode number [1 or 2]: ").strip()
    mode = "cas9" if mode_input != "2" else "abe"

    edit_window = (1, 10)
    if mode == "abe":
        ew_input = input("Enter editing window (e.g. 1-10, optional) [default: 1-10]: ").strip()
        if ew_input:
            edit_window = parse_editing_window(ew_input)

    seq = input("Enter the gene sequence: ").strip()
    gene = input("Enter the gene name (optional): ").strip()
    path = input("Enter output path (optional) [default: .]: ").strip()

    filename = f"cheetah_output_{gene}.txt" if gene else "cheetah_output.txt"
    if not path:
        path = os.path.join(os.getcwd(), filename)
    else:
        path = os.path.expanduser(path)
        if os.path.isdir(path):
            path = os.path.join(path, filename)

    rows = design_targets(seq, gene, mode, edit_window)
    if len(rows) == 1:
        print("No valid euRptr targets found.")
    else:
        save_to_txt(rows, path)
        print(f"Saved {len(rows) - 1} targets to '{path}'.")

def run_peurptr_design():
    print("\n[peuRptr design mode]")

    antirepeat = input("Enter euRptr antirepeat sequence: ").strip()
    rtt = input("Enter RTT sequence: ").strip()
    pbs = input("Enter PBS sequence: ").strip()
    linker = input("Enter Linker sequence (only if using epegRNA form): ").strip()

    if not antirepeat or not rtt or not pbs:
        print("Error: antirepeat, RTT, and PBS are required.")
        return

    scaffold = "aaggctagtccgttatcaacttGGACTTCGGTCCaagtggcaccgagtcggtgc"
    tevopreq1 = "cgcggttctatctagttacgcgttaaaccaactagaa"

    full_seq = antirepeat + scaffold + rtt + pbs + linker + tevopreq1
    columns = ["antirepeat", "scaffold", "RTT", "PBS", "linker", "tevopreQ1"]
    values = [antirepeat, scaffold, rtt, pbs, linker, tevopreq1]

    print("\npeuRptr sequence:\t" + full_seq)
    print("\t".join(columns))
    print("\t".join(values))


def run_cheetah():
    print("CHEETAH euRptr designer (Baelab)")
    print("Select mode:")
    print("1. euRptr antirepeat design")
    print("2. peuRptr design")
    mode = input("Enter mode number [1 or 2]: ").strip()

    if mode == "2":
        run_peurptr_design()
    else:
        run_eurptr_design()

if __name__ == "__main__":
    run_cheetah()
