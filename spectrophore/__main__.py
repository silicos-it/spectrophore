"""
Convert an sdf file to spectrophore fingerprints
The input sdf file should contain 3D poses with explicit hydrogens
"""

import argparse
import os
import time

from concurrent import futures
from rdkit import Chem
from spectrophore import spectrophore
from tqdm import tqdm


def __check_resolution(value):
    """
    Check if the resolution value is valid
    """
    v = float(value)
    if v <= 0:
        raise argparse.ArgumentTypeError("should be larger than 0")
    return v


def __processCommandline():
    """
    Process the command-line arguments
    """
    parser = argparse.ArgumentParser(
        prog="spectrophore",
        description="Convert a SDF-file containing 3D structures to a csv with the corresponding spectrophores fingerprints",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-n",
        "--norm",
        dest="norm",
        help="Normalization setting",
        default="all",
        choices=["none", "mean", "all", "std"],
    )
    parser.add_argument(
        "-m",
        "--mode",
        dest="mode",
        help="Which type of spectrophore to calculate",
        default="classic",
        choices=["classic", "full"],
    )
    parser.add_argument(
        "-s",
        "--stereo",
        dest="stereo",
        help="Stereo setting",
        default="none",
        choices=["none", "unique", "mirror", "all"],
    )
    parser.add_argument(
        "-a",
        "--accuracy",
        dest="accuracy",
        help="Accuracy setting",
        type=int,
        default=20,
        choices=[1, 2, 3, 4, 5, 6, 9, 10, 12, 15, 18, 20, 30, 36, 45, 60, 90, 180],
    )
    parser.add_argument(
        "-r",
        "--resolution",
        dest="resolution",
        help="Resolution setting (>0)",
        default=3,
        type=__check_resolution,
    )
    parser.add_argument(
        "-p",
        "--np",
        dest="max_workers",
        help="Number of processors to use; -1 is all processors",
        default=-1,
        type=int,
    )
    parser.add_argument(
        "--silent",
        dest="silent",
        help="Don't show a progressbar",
        default=False,
        action="store_true",
    )

    requiredNamed = parser.add_argument_group("required arguments")
    requiredNamed.add_argument(
        "-i", "--in", dest="infile", help="Input sdf file", required=True
    )
    requiredNamed.add_argument(
        "-o",
        "--out",
        dest="outfile",
        help="Output spectrophore file",
        required=True,
    )
    return parser.parse_args()


def main():
    args = __processCommandline()

    args.outfile = os.path.abspath(args.outfile)

    calculator = spectrophore.SpectrophoreCalculator(
        normalization=args.norm,
        mode=args.mode,
        stereo=args.stereo,
        accuracy=args.accuracy,
        resolution=args.resolution,
    )
    print("Normalization: ", args.norm)
    print("Mode:          ", args.mode)
    print("Stereo:        ", args.stereo)
    print("Accuracy:      ", args.accuracy)
    print("Resolution:    ", args.resolution)
    print("Input file:    ", args.infile)
    print("Output file:   ", args.outfile)
    print("Processors:    ", args.max_workers)

    startTime = time.time()

    successCounter = 0

    supplier = Chem.SDMolSupplier(args.infile, removeHs=False)
    of = open(args.outfile, "w")
    mw = None

    if args.max_workers > 0:
        mw = args.max_workers

    with futures.ProcessPoolExecutor(max_workers=mw) as executor:
        jobs = []
        for mol in supplier:
            if mol is None:
                continue

            job = executor.submit(calculator.calculate_string, mol)
            jobs.append(job)

        for job in tqdm(
            futures.as_completed(jobs), total=len(jobs), disable=args.silent
        ):
            spec = job.result()
            if spec is not None:
                of.write(f"{spec}\n")
                successCounter += 1

    of.close()

    print("\nSpectrophore calculations complete!")
    print(
        f"* Successfully processed {successCounter} out of {len(supplier)} compounds ({(successCounter/len(supplier))*100:.2f}%)"
    )
    print(f"* Time elapsed: {time.time()-startTime:.2f} seconds")
    print(f"* Output file written to: {args.outfile}")


if __name__ == "__main__":
    main()
