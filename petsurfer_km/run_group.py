"""Main entry point for petsurfer-km-group CLI."""

from __future__ import annotations

# Use pysqlite3 as a drop-in replacement for sqlite3 if available.
# This is needed for FreeSurfer's fspython, which lacks the _sqlite3 C extension.
# pysqlite3-binary only has wheels for Linux x86-64
#
# Be sure to run `fspython -m pip install pysqlite3-binary` if running from
# inside fspython
try:
    __import__('pysqlite3')
    import sys
    sys.modules['sqlite3'] = sys.modules.pop('pysqlite3')
except ImportError:
    pass

import argparse
import logging
import os
import shutil
import sys
import subprocess
from pathlib import Path

import petsurfer_km
from petsurfer_km import __version__
from petsurfer_km.cli.parser import build_parser
from petsurfer_km.inputs import InputGroup, discover_inputs
from petsurfer_km.steps import (
    run_bidsify,
    run_kinetic_modeling,
    run_preprocessing,
    run_report,
    run_surface,
    run_volumetric,
)
from bids import BIDSLayout

logger = logging.getLogger("petsurfer_km")

def setup_logging(level: str) -> None:
    """
    Configure logging for petsurfer-km-group

    Args:
        level: Log level string ('error', 'warn', 'info', or 'debug')
    """
    level_map = {
        "error": logging.ERROR,
        "warn": logging.WARNING,
        "info": logging.INFO,
        "debug": logging.DEBUG,
    }
    log_level = level_map.get(level, logging.WARNING)

    # Configure format based on level
    if log_level == logging.DEBUG:
        fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    else:
        fmt = "%(levelname)s: %(message)s"

    logging.basicConfig(
        level=log_level,
        format=fmt,
        handlers=[logging.StreamHandler(sys.stderr)],
    )


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """
    Validate parsed arguments for consistency.

    Args:
        args: Parsed arguments namespace.
        parser: ArgumentParser instance for error reporting.
    """
    # Check that tstar is provided for Logan methods
    logan_methods = {"logan", "logan-ma1"}
    selected_logan = logan_methods.intersection(args.km_method)
    if selected_logan and args.tstar is None:
        parser.error(
            f"--tstar is required when using {', '.join(sorted(selected_logan))} method(s)"
        )

    # Cannot disable both volumetric and surface analysis
    if args.no_vol and args.no_surf:
        parser.error(
            "Cannot disable both volumetric (--no-vol) and surface (--no-surf) analysis"
        )

    # Check bloodstream-dir exists for Logan methods (uses arterial input function)
    if selected_logan:
        bloodstream_dir = args.bloodstream_dir
        if bloodstream_dir is None:
            bloodstream_dir = args.bids_dir / "derivatives" / "bloodstream"
        if not bloodstream_dir.exists():
            parser.error(
                f"bloodstream-dir does not exist: {bloodstream_dir}\n"
                f"Logan methods require arterial input function data from bloodstream."
            )

    # Check petprep-dir exists
    petprep_dir = args.petprep_dir
    if petprep_dir is None:
        petprep_dir = args.bids_dir / "derivatives" / "petprep"
    if not petprep_dir.exists():
        parser.error(
            f"petprep-dir does not exist: {petprep_dir}\n"
            f"Please run petprep first or specify --petprep-dir."
        )


def set_defaults(args: argparse.Namespace) -> argparse.Namespace:
    """Set default values that depend on other arguments."""
    # Set default petprep-dir
    if args.petprep_dir is None:
        args.petprep_dir = args.bids_dir / "derivatives" / "petprep"

    # Set default bloodstream-dir
    if args.bloodstream_dir is None:
        args.bloodstream_dir = args.bids_dir / "derivatives" / "bloodstream"

    # Set default work-dir
    if args.work_dir is None:
        args.work_dir = Path(f"/tmp/petsurfer-km-{os.getpid()}")

    # Handle hemisphere selection
    if not args.lh and not args.rh:
        # Default: process both hemispheres
        args.hemispheres = ["lh", "rh"]
    else:
        args.hemispheres = []
        if args.lh:
            args.hemispheres.append("lh")
        if args.rh:
            args.hemispheres.append("rh")

    return args


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse and validate command-line arguments.

    Args:
        argv: Command-line arguments (defaults to sys.argv[1:])

    Returns:
        Validated namespace of arguments.

    Raises:
        SystemExit: On parsing or validation errors.
    """

    parser = argparse.ArgumentParser(
        prog="petsurfer-km-group",
        description=("BIDS App to perform group analysis on PET kinetic model output of PetSurfer")
    )
    parser.add_argument("bids_dir",type=Path,help="Root directory of the BIDS dataset.");
    parser.add_argument("output_dir",type=Path,help="Output directory for results.");
    parser.add_argument("--km",
        choices=["MRTM1", "MRTM2", "Logan", "MA1"],
        help=("Kinetic modeling method to get data from. Only one method can be specified. "
            "MRTM1, MRTM2, Logan, MA1. "))
    parser.add_argument('--ses',default="baseline",help="Session label.")
    parser.add_argument('--tracer',default="11CPS13",help="Tracer to analyze (eg,11CPS13).")
    parser.add_argument("--space",choices=["fsaverage", "MNI152NLin2009cAsym","ROI"],help=("Group space(s)."))
    parser.add_argument("--hemi", choices=["L", "R"],help="Single hemisphere when using fsaverage.")
    parser.add_argument('--fwhm',help="FWHM when running the kinetic modeling.")

    args = parser.parse_args(argv)

    if(args.space == "fsaverage" and args.hemi is None):
        print("ERROR: must spec hemi with fsaverage");
        sys.exit(1);

    print(f"bids_dir {args.bids_dir}")
    print(f"output_dir {args.output_dir}")
    print(f"km {args.km}")
    print(f"ses {args.ses}")
    print(f"tracer {args.tracer}")
    print(f"space {args.space}")
    print(f"hemi {args.hemi}")
    print(f"fwhm {args.fwhm}")

    #args = set_defaults(args)
    #validate_args(args, parser)
    return args

def ensure_fsaverage() -> None:
    """
    Set SUBJECTS_DIR to $FREESURFER/subjects to ensure fsaverage exists

    Raises:
        RuntimeError: If FREESURFER_HOME is not set
    """

    freesurfer_home = os.environ.get("FREESURFER_HOME")
    if not freesurfer_home:
        raise RuntimeError(f"FREESURFER_HOME is not set. Cannot locate fsaverage.")

    subjects_dir = Path(freesurfer_home) / "subjects"
    os.environ.get("SUBJECTS_DIR",subjects_dir)

def run(args: argparse.Namespace) -> int:
    """
    Execute petsurfer-km processing.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Exit code (0 for success, non-zero for errors).
    """
    logger.info(f"petsurfer_km_group version: {__version__}")
    logger.debug(f"BIDS directory: {args.bids_dir}")
    logger.debug(f"Output directory: {args.output_dir}")
    logger.debug(f"Kinetic methods: {args.km_method}")
    logger.debug(f"PetPrep directory: {args.petprep_dir}")
    logger.debug(f"Work directory: {args.work_dir}")
    if args.session_label:
        logger.debug(f"Sessions: {args.session_label}")
    logger.debug(f"Volumetric FWHM: {args.vol_fwhm} mm")
    logger.debug(f"Surface FWHM: {args.surf_fwhm} mm")
    logger.debug(f"Hemispheres: {args.hemispheres}")
    logger.debug(f"Skip volumetric: {args.no_vol}")
    logger.debug(f"Skip surface: {args.no_surf}")

    # Create working directory
    args.work_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created working directory: {args.work_dir}")

    # Ensure fsaverage is available under SUBJECTS_DIR
    ensure_fsaverage()

    # Track processing results
    failed_subjects: list[str] = []
    successful_subjects: list[str] = []

    # Process each subject/session
    for input_group in input_groups:
        try:
            process_subject(input_group, args)
            successful_subjects.append(input_group.label)
        except Exception as e:
            failed_subjects.append(input_group.label)
            logger.error(f"Failed to process {input_group.label}: {e}")
            if args.abort_on_error:
                logger.error("Aborting due to --abort-on-error")
                break

    # Log summary
    if successful_subjects:
        logger.info(f"Successfully processed: {', '.join(successful_subjects)}")
    if failed_subjects:
        logger.error(f"Failed to process: {', '.join(failed_subjects)}")

    # Cleanup
    if args.nocleanup:
        logger.info(f"Skipping cleanup (--nocleanup). Working directory: {args.work_dir}")
    else:
        logger.info(f"Cleaning up working directory: {args.work_dir}")
        try:
            shutil.rmtree(args.work_dir)
        except Exception as e:
            logger.warning(f"Failed to clean up working directory: {e}")

    # Return exit code
    if failed_subjects:
        return 1
    return 0


def main(argv: list[str] | None = None) -> None:
    """
    Main entry point for petsurfer-km-group CLI.

    Args:
        argv: Command-line arguments (defaults to sys.argv[1:])
    """
    args = parse_args(argv)

    if(args.km in ["MRTM1", "MRTM2"]): meas = "BPND";
    if(args.km in ["Logan", "MA1"]):   meas = "VT";
    print(meas)

    desc = "sm"+args.fwhm;
    print(desc)

    #pskmconfig = "/autofs/space/iddhi_005/users/greve/petsurfer-bids/petsurfer_km/petsurfer-km-bids-config.json";
    pkgdir = os.path.dirname(petsurfer_km.__file__)
    pskmconfig = os.path.join(pkgdir, "petsurfer-km-bids-config.json")
    layout = BIDSLayout(args.bids_dir,validate=False,config=["bids", "derivatives", pskmconfig]);

    fshemi = None;
    if(args.hemi == "L"): fshemi = "lh";
    if(args.hemi == "R"): fshemi = "rh";

    gfiles = layout.get(
        session=args.ses,
        datatype="pet",
        tracer=args.tracer,
        hemi=args.hemi,
        space=args.space,
        desc=desc,
        model=args.km,
        meas=meas,
        suffix="mimap",
        extension=".nii.gz",
        return_type="filename")
    nfiles = len(gfiles);
    print(f"space={args.space} hemi={args.hemi} n={nfiles} ============")

    if(fshemi == None):
        stack = Path(args.output_dir,args.space+".nii.gz")
        glmdir = Path(args.output_dir,"glm.mni152")
    else:
        stack = Path(args.output_dir,"fsaveage."+fshemi+".nii.gz")
        glmdir = Path(args.output_dir,"glm.fsaverage."+fshemi);
    cmd = f"mri_concat --o {stack} {' '.join(map(str, gfiles))}"
    print(cmd)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    try:
        result = subprocess.run(cmd.split(), check=True, capture_output=True, text=True)
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
    except subprocess.CalledProcessError as e:
        print("Command failed with exit code", e.returncode)
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        sys.exit(1);

    cmd = f"mri_glmfit --o {glmdir} --y {stack} --osgm"
    if(args.space == "fsaverage"): cmd = f"{cmd} --surf fsaverage {fshemi}"
    print(cmd)
    try:
        result = subprocess.run(cmd.split(), check=True, capture_output=True, text=True)
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
    except subprocess.CalledProcessError as e:
        print("Command failed with exit code", e.returncode)
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        sys.exit(1);

    sys.exit(0);


if __name__ == "__main__":
    main()
