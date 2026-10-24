"""In-silico PCR ringan untuk dashboard (anchor 3' 12-mer + cek mismatch)."""
RC = str.maketrans("ACGT", "TGCA")

def revcomp(s):
    return s.translate(RC)[::-1]

def mismatches(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)

def count_amplicons(seqs, fwd, rev, max_mm=2, min_amp=50, max_amp=1000):
    """seqs: iterable urutan (string). Hitung amplikon valid."""
    rrev = revcomp(rev)
    fa, ra = fwd[-12:], rrev[-12:]
    count = 0
    for genome_seq in seqs:
        i = genome_seq.find(fa)
        while i != -1:
            s = i + 12 - len(fwd)
            if s >= 0:
                sub = genome_seq[s:s + len(fwd)]
                if len(sub) == len(fwd) and mismatches(sub, fwd) <= max_mm and mismatches(sub[-3:], fwd[-3:]) == 0:
                    j = genome_seq.find(ra, s)
                    while j != -1 and j - s <= max_amp + len(rev):
                        rs = j + 12 - len(rev)
                        if rs >= 0:
                            rsub = genome_seq[rs:rs + len(rev)]
                            if len(rsub) == len(rev) and mismatches(rsub, rrev) <= max_mm and mismatches(rsub[-3:], rrev[-3:]) == 0:
                                amp = j + 12 - s
                                if min_amp <= amp <= max_amp:
                                    count += 1
                        j = genome_seq.find(ra, j + 1)
            i = genome_seq.find(fa, i + 1)
    return count
