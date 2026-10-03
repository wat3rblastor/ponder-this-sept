import json, sys
K, M = sys.argv[1], sys.argv[2]; T = f"K{K}_m{M}"
eh = set(); nun = 0
for l in open(f"en_{T}.jsonl"):
    j = json.loads(l)
    if j.get("hit"):
        assert j.get("mode", 0) == int(M); eh.add((j["n"], j["a"], j["d"]))
    else:
        assert j["mode"] == int(M); nun += 1
oh = set(tuple(map(int, l.split())) for l in open(f"or_{T}.hits"))
ew = set(int(l.split()[2]) for l in open(f"en_{T}.win"))
ow = set(int(l) for l in open(f"or_{T}.win"))
print(f"{T}: engine units={nun} hits={len(eh)} oracle hits={len(oh)} common={len(eh&oh)} engine-only={len(eh-oh)} oracle-only={len(oh-eh)}  -> HITS {'MATCH' if eh==oh else 'DIFFER'}")
print(f"{T}: stage-3 windows engine={len(ew)} oracle sieve-survivors={len(ow)} engine-only={len(ew-ow)} oracle-only={len(ow-ew)} -> WINDOWS {'MATCH' if ew==ow else 'DIFFER'}")
