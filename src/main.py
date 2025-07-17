import json
import sys
from pathlib import Path

from ipinstance import IPInstance
from model_timer import Timer


def main(filepath: str):

    filename = Path(filepath).name
    watch = Timer()
    watch.start()
    solver = IPInstance(filepath)
    solution = solver.solve()
    inst_str = solver.toString()
    watch.stop()
    print(inst_str)

    sol_dict = {
        "Instance": filename,
        "Time": str(round(watch.getElapsed(), 2)),
        "Result": str(solution),
        "Solution": "OPT"
    }
    print(json.dumps(sol_dict))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python main.py <input_file>")
    main(sys.argv[1])
