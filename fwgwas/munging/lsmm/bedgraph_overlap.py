#!/usr/bin/env python3

"""Retrieves values for annotations overlapping SNP coordinates.

Provides a '0' or '1' value for every SNP record in a single chromosome SNP
coordinates file. A '1' value indicates a given SNP overlaps with a genomic
region record that has a data value of '1' in a BedGraph file. A '0'
indicates no overlap. Records with data values of '0' in the BedGraph file are
ignored. All coordinates are assumed to be for the same chromosome. Genomic
regions in the BedGraph file are assumed to be non-overlapping and sorted.
"""

import argparse
import os.path

def get_args():
  """Parses command line arguments.

  Returns:
    A Namespace object with argument variables as attributes.
  """
  parser = argparse.ArgumentParser(
    description=(
      "Provides a '0' or '1' value for every SNP record in a single "
      "chromosome SNP coordinates file. A '1' value indicates a given "
      "SNP overlaps with a genomic region record that has a data value "
      "of '1' in a BedGraph file. A '0' indicates no overlap. Records "
      "with data values of '0' in the BedGraph file are ignored. All "
      "coordinates across files are assumed to be for the same "
      "chromosome."
    )
  )
  parser.add_argument(
    "snp_file",
    type=str,
    help=(
      "SNP file with chromosome positions in the 3rd column. "
      "Header line assumed."
    )
  )
  parser.add_argument(
    "anno_file",
    type=str,
    help="BedGraph format file. Header line assumed."
  )
  parser.add_argument(
    "anno_name",
    type=str,
    help="Header name for the output file."
  )
  parser.add_argument(
    "output_file",
    type=str,
    help="Name of the output file."
  )
  args = parser.parse_args()
  if not os.path.isfile(args.snp_file):
    raise SystemError(f"{args.snp_file} not found.")
  if not os.path.isfile(args.anno_file):
    raise SystemError(f"{args.anno_file} not found.")
  return args

def import_bedgraph(filename):
  """Imports BedGraph file data.

  Reads in BedGraph records with data values equal to 1. Requires the
  records to be non-overlapping and sorted by increasing position.

  Args:
    filename: str; BedGraph file with a header line.

  Returns:
    A tuple of ranges as tuples (i.e., tuple of tuples).
  """
  intervals = []
  with open(filename, "r") as fin:
    fin.readline() #assumes first line is a header
    for line in fin:
      items = line.strip().split("\t")
      if items[3] == "1":
        start = int(items[1])
        end = int(items[2])
        intervals.append((start, end))
  for i in range(1, len(intervals)):
    assert intervals[i][0] >= intervals[i-1][1], \
      (f"Intervals {intervals[i-1]} and {intervals[i]} overlap or are "
       "not in increasing order!")
  return tuple(intervals)

def has_overlap(query, db):
  """Checks if a value lies in a set of intervals.

  Identifies whether an integer value resides in any of the provided
  intervals. Assumes that an interval is in half-open format, so the end
  position value is not considered as part of the interval. Additionally,
  it is assumed that the intervals are non-overlapping and sorted by
  increasing position.

  Args:
    query: int; A value to check against the interval set.
    db: list-like; The interval set with start and end positions as first
      and second elements respectively in list-like data structures.

  Returns:
    A boolean. 'True' if the value resides in any interval,
    otherwise 'False'.
  """
  # Binary search works well since intervals are sorted and non-overlapping
  low = 0
  high = len(db) - 1
  mid = 0
  while low <= high:
    mid = (high + low) // 2
    if db[mid][1] <= query:
      low = mid + 1
    elif db[mid][0] > query:
      high = mid - 1
    else:
      return True
  return False

def main():
  """Main function for runtime execution."""
  args = get_args()
  intervals = import_bedgraph(args.anno_file)
  anno_values = []
  with open(args.snp_file, "r") as fin:
    fin.readline() #assumes first line is a header
    for line in fin:
      items = line.strip().split("\t")
      indicator = has_overlap(int(items[2]), intervals)
      anno_values.append(int(indicator))
  with open(args.output_file, "w") as fout:
    fout.write(args.anno_name + "\n")
    for value in anno_values:
      fout.write(f"{value}\n")

if __name__ == "__main__":
  main()
