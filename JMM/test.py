import sys
import time

from pathlib import Path

print('one line')
print('two lines')

time.sleep(1)
print('slept 1 second')

with open(Path(sys.argv[1])) as f:
    i = f.readline()
    print(f'read integer from Python: {int(i)}')
