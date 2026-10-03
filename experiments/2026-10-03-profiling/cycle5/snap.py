# snap.py TAG -> prints per-file line counts as json
import glob,sys,json
d={f:sum(1 for _ in open(f)) for f in sorted(glob.glob(f'/workspace/ponder-this-sept/experiments/remote/{sys.argv[1]}_s*.jsonl'))}
print(json.dumps(d))
