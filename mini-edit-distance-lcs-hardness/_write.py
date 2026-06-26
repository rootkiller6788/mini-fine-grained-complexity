import sys, os
BASE = os.path.dirname(os.path.abspath(__file__))
if len(sys.argv) < 3:
    print("Usage: python _write.py <relpath> <content>")
    sys.exit(1)
path = os.path.join(BASE, sys.argv[1])
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, 'a', encoding='utf-8') as f:
    f.write(sys.argv[2])
    f.write('\n')
print(f"appended to {sys.argv[1]}")
