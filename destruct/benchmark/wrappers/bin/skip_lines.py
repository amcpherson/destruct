import sys

num_lines = int(sys.argv[1])

for idx, line in enumerate(sys.stdin):
    if idx < num_lines:
        continue
    try:
        sys.stdout.write(line)
    except IOError as e:
        if e.errno == 32:
            break
