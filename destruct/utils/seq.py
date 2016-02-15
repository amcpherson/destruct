

def read_sequences(fasta):
    id = None
    sequences = []
    for line in fasta:
        line = line.rstrip()
        if len(line) == 0:
            continue
        if line[0] == '>':
            if id is not None:
                yield (id, ''.join(sequences))
            id = line[1:].split()[0]
            sequences = []
        else:
            sequences.append(line)
    if id is not None:
        yield (id, ''.join(sequences))

