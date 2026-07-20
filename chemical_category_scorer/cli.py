from __future__ import annotations

import argparse

from . import available_models, details_smiles


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Chemical Category Scorer command line interface.")
    parser.add_argument("--list-models", action="store_true", help="List built-in public model ids.")
    parser.add_argument("--score", help="Score a SMILES string.")
    parser.add_argument("--model-id", default="final_pesticides", help="Model id used with --score.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.list_models:
        for model_id in available_models():
            print(model_id)
        return
    if args.score:
        print(details_smiles(args.score, model_id=args.model_id))
        return
    raise SystemExit("Use --list-models or --score.")


if __name__ == "__main__":
    main()
