"""Run the SA subset search on a list of universes; save results/<name>.json"""
import sys, os, json, time
import numpy as np
import universe as uv
from universes_list import UNIVERSES
from search import run_universe, formula

os.makedirs('results', exist_ok=True)

def main(names, ns, runs, iters, seed=0):
    for name in names:
        t = time.time()
        U = UNIVERSES[name]()
        print(U.summary(), f"built in {time.time()-t:.1f}s", flush=True)
        t = time.time()
        res = run_universe(U, ns, runs=runs, iters=iters, seed=seed)
        out = {'universe': name, 'summary': U.summary(), 'N': U.N, 'runs': runs, 'iters': iters, 'seconds': round(time.time() - t, 1),
               'results': {str(n): [{'value': v, 'subset': list(k), 'O': o, 'formula': formula(n)} for v, k, o in items[:5]] for n, items in res.items()}}
        fn = f'results/{name}_s{seed}.json'
        json.dump(out, open(fn, 'w'), indent=1)
        print(f"{name}: {time.time()-t:.1f}s -> {fn}", flush=True)

if __name__ == "__main__":
    names = sys.argv[1].split(',') if len(sys.argv) > 1 and sys.argv[1] != 'all' else list(UNIVERSES)
    ns = range(6, 21)
    runs = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    iters = int(sys.argv[3]) if len(sys.argv) > 3 else 300_000
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    main(names, ns, runs, iters, seed)
