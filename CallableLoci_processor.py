#!/usr/bin/env python3

import sys
if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3 to run.")

import os
import argparse
import pathlib
import errno

parser = argparse.ArgumentParser(description = "This script processes CallableLoci files, to generate merged sorted intersections of sample and control files. It also calculates the surveyed region of the genome. Requires bedtools to be in your path.")
parser.add_argument("dir_in", help = "The input directory. This directory should contain a seperate folder for each sample, whith the CallableLoci outputs")
parser.add_argument("dir_out", default = os.getcwd(), help = "The output directory")
parser.add_argument("name", help = "The name of the experiment")
parser.add_argument("control", help = "The name of the control sample")
parser.add_argument("--samples", nargs = "*", help = "The names of the samples. The script will normally try to get these from a vcf, but they can also be supplied manually.")
parser.add_argument("--vcf", help = "The name of the vcf from which to get the samples. If not provided the script will try to guess the name of the vcf. Alternatively the sample names can be supplied manually")
parser.add_argument("--version", action='version', version="%(prog)s 1.0.0")
parser.add_argument("-v", "--verbose", action = "store_true")
parser.add_argument("--overwrite", action = "store_true")
args = parser.parse_args()

dir_in = args.dir_in
exp = args.name
bulk = args.control
out = args.dir_out

#Get the sample names
if args.samples:
    samples = args.samples
else:
    if args.vcf:
        vcf = args.vcf
    else:
        vcf = dir_in + exp + ".filtered_variants.vcf"
    with open(vcf) as vcf_file:
        for line in vcf_file:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                line = line.strip()
                header = line.split("\t")
                samples_bulk = header[9:]
                samples = samples_bulk.copy()
                samples.remove(bulk)
                break

#Create output directory
out_dir = out + exp + "_CallableLoci/"
pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

#Function to determine whether or not an output should be written
def write_or_not(fname):
    write = not os.path.exists(fname) or args.overwrite
    return write

#Get callable loci for the bulk
bulk_bed = dir_in + bulk + "/" + bulk + "_CallableLoci.bed"
bulk_callable_fname = out_dir + bulk + "CallableLoci_CALLABLE.bed"
if not os.path.isfile(bulk_bed):
    bulk_bed = dir_in + bulk + "_dedup.realigned/" + bulk + "_dedup.realigned_CallableLoci.bed"
if os.path.isfile(bulk_bed) and write_or_not(bulk_callable_fname):
    command = "grep 'CALLABLE' {0} > {1}".format(bulk_bed, bulk_callable_fname)
    os.system(command)
    if args.verbose:
        print("Filtered CALLABLE file for exp: {0}".format(exp))

bulk_autosomal_fname = out_dir + bulk + "CallableLoci_autosomal.bed"
if os.path.isfile(bulk_callable_fname) and write_or_not(bulk_autosomal_fname):
    command = "sed '/[XYMT]/d' {0} > {1}".format(bulk_callable_fname, bulk_autosomal_fname)
    os.system(command)
    if args.verbose:
        print("Filtered for autosomal regions on the CALLABLE bulk file for exp: {0}".format(exp))


#Get callable regions of each sample and merge them with the callable regions of the bulk
for sample in samples:

    #GEt callable region
    bed_fname = dir_in + sample + "/" + sample + "_CallableLoci.bed"
    bed_callable_fname = out_dir + sample + "CallableLoci_CALLABLE.bed"
    if not os.path.isfile(bed_fname):
        bed_fname = dir_in + sample + "_dedup.realigned/" + sample + "_dedup.realigned_CallableLoci.bed"
    if os.path.isfile(bed_fname) and write_or_not(bed_callable_fname):
        command = "grep 'CALLABLE' {0} > {1}".format(bed_fname, bed_callable_fname)
        os.system(command)
        if args.verbose:
            print("Created CALLABLE file for sample {0} in exp {1}".format(sample, exp))

    #Intersect bed sample with bed bulk
    bed_intersected = out_dir + sample + "_" + bulk + "_CallableLoci.bed"
    if os.path.isfile(bed_callable_fname) and os.path.isfile(bulk_callable_fname) and write_or_not(bed_intersected):
        command = "intersectBed -a {0} -b {1} > {2}".format(bed_callable_fname, bulk_callable_fname, bed_intersected)
        os.system(command)
        if args.verbose:
            print("Intersected Callable bed file for sample {0} in exp {1} with the bulk {2}".format(sample, exp, bulk))

    #Sort the bedfile
    bed_sorted = out_dir + sample + "_" + bulk + "_CallableLoci_sorted.bed"
    if os.path.isfile(bed_intersected) and write_or_not(bed_sorted):
        command = "sortBed -i {0} > {1}".format(bed_intersected, bed_sorted)
        os.system(command)
        if args.verbose:
            print("Sorted intersected bed for sample {0} in exp {1}".format(sample, exp))

    #Merge overlapping regions of the bedfile
    bed_merged = out_dir + sample + "_" + bulk + "_CallableLoci_merged.bed"
    if os.path.isfile(bed_sorted) and write_or_not(bed_merged):
        command = "mergeBed -i {0} > {1}".format(bed_sorted, bed_merged)
        os.system(command)
        if args.verbose:
            print("Merged overlapping regions in bed for sample {0} in exp {1}".format(sample, exp))

    #Filter for autosomal regions
    bed_merged_autosomal = out_dir + sample + "_" + bulk + "_CallableLoci_merged_autosomal.bed"
    if os.path.isfile(bed_merged) and write_or_not(bed_merged_autosomal):
        command = "sed '/[XYMT]/d' {0} > {1}".format(bed_merged, bed_merged_autosomal)
        os.system(command)
        if args.verbose:
            print("Filtered for autosomal sites for sample {0} in exp {1}".format(sample, exp))

    #Count surveyed region of the genome
    surveyed_region = out_dir + sample + "_" + bulk + "_surveyed.txt"
    if os.path.isfile(bed_merged_autosomal) and write_or_not(surveyed_region):
        command = "cat " + bed_merged_autosomal + " | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' > " + surveyed_region
        os.system(command)
        if args.verbose:
            print("Counted surveyed region for sample {0} in exp {1}".format(sample, exp))
