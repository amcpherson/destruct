import sys

max_read_length = 0

for line in sys.stdin:
    fields = line.split()
    if len(fields) < 10:
        continue
    max_read_length = max(max_read_length, len(fields[9]))

print max_read_length
