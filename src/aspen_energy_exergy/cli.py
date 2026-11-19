import argparse
from pathlib import Path

from .api import analyze, plan_reference_flashes, write_native_exergy_template
from .flash_exchange import write_flash_plan
from .report import summary, write_results


def main(argv=None):
    parser = argparse.ArgumentParser(description="Analyze a saved Aspen Plus case or normalized JSON export.")
    parser.add_argument("source", type=Path, help="Solved .apwz, .bkp, or normalized .json")
    parser.add_argument("--settings", help="Settings JSON path, biogas-ft, or as-saved (raw transfers)")
    parser.add_argument("--references", action="append", type=Path, metavar="FILE",
                        help="Matching ambient-flash JSON or flash_returns.csv (repeat for partial batches)")
    parser.add_argument("--native-exergy", type=Path, metavar="XLSX_OR_CSV",
                        help="Full Aspen material workbook or completed EXERGYFL return sheet")
    parser.add_argument("--report-text", type=Path, metavar="TXT", help="Aspen text report to cross-check with the workbook")
    parser.add_argument("--native-splits-only", action="store_true", help="Disable library-estimated splits")
    parser.add_argument("--balance-engine", choices=("builtin","exerpy-tespy"), default="builtin",
                        help="Evaluate imported Aspen balances using the pinned ExerPy and TESPy engines")
    outputs = parser.add_mutually_exclusive_group()
    outputs.add_argument("-o", "--output", type=Path, help="Results ZIP (default: aspen-analysis.zip)")
    outputs.add_argument("--write-exergy-template", type=Path, metavar="CSV",
                         help="Create a blank EXERGYFL return sheet for this case and exit")
    outputs.add_argument("--write-flash-plan", type=Path, metavar="ZIP_OR_JSON",
                         help="Group outstanding ambient flashes into a request packet or fillable JSON")
    args = parser.parse_args(argv)
    destination = args.write_flash_plan or args.write_exergy_template or args.output or Path("aspen-analysis.zip")
    if args.write_exergy_template and (args.references or args.native_exergy or args.report_text or args.settings or args.native_splits_only):
        parser.error("--write-exergy-template uses only the source case and destination CSV")
    if not (args.write_exergy_template or args.write_flash_plan) and destination.suffix.lower() != ".zip":
        parser.error("Results destination must end in .zip")
    inputs = [args.source, args.native_exergy, args.report_text, *(args.references or [])]
    if args.settings and args.settings not in {"biogas-ft", "as-saved"}:
        inputs.append(Path(args.settings))
    if any(p and p.resolve() == destination.resolve() for p in inputs):
        parser.error("Results destination must differ from every input")
    try:
        if args.write_flash_plan:
            plan = plan_reference_flashes(args.source, settings=args.settings, references=args.references,
                                          native_exergy=args.native_exergy, report_text=args.report_text)
            write_flash_plan(plan, destination)
            coverage = plan["summary"]
            print(f"{coverage['flash_groups']} candidate ambient flashes for {coverage['plannable_streams']} streams.")
            print(f"Flash plan saved to {destination.resolve()}")
            return 0
        if args.write_exergy_template:
            path = write_native_exergy_template(args.source, destination)
            print(f"EXERGYFL return sheet saved to {path.resolve()}")
            return 0
        result = analyze(args.source, settings=args.settings, references=args.references,
                         estimate_splits=not args.native_splits_only, native_exergy=args.native_exergy, report_text=args.report_text,
                         balance_engine=args.balance_engine)
        write_results(result, destination)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, f"Analysis failed: {exc}\n")
    print(summary(result), end="")
    print(f"Results saved to {destination.resolve()}")
    return 0
