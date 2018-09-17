#!/usr/bin/env python
# justin le tourneau Sep 16 2018
# transforms Illumina Genotyping Report into a PED file with the help of SNP_Map file
# version 1.0

import argparse, sys, os, subprocess, re, datetime

def print_error(message):
    print "ERROR: %s " % message
    sys.exit()

def check_file(value, name):
    try:os.path.exists(value)
    except:print_error("%s was not found" % name)

def create_sorted_SNP_Map(file):
	print "Sorting SNP_Map"
	# remove the header line
	# sort and get unique npm_map file name's
	os.system('sed 1d %s | sort -k2 | uniq -f1 -i -c > /tmp/tmp_snp_map.txt' % file)

def get_SNP_Map_marker_ids():
	print "Loading SNP_Map marker ids"
	# markers to return
	markers = []
	dupes = []
	# get the duplicates
	# found_duplicates = subprocess.check_output(['grep "^\s*\([2-9]\|\d\d\d*\)\s" /tmp/tmp_snp_map.txt | awk \'{$1=""}1\''], shell=True)
	marker_ids = subprocess.check_output(['awk \'{print $1 " " $3}\' /tmp/tmp_snp_map.txt'], shell=True)
	marker_ids = marker_ids.rstrip().split('\n')
	# go through each marker
	for marker in marker_ids:
		mark = marker.split(' ')
		# if it is only in the SNP file once, add it
		markers.append(mark[1]) if mark[0] == '1' else dupes.append(mark[1])
	#warn the user of dupes
	if len(dupes) > 0:
		print "\nDuplicates in SNP Map found. \nThe following will be omitted from the PED file: \n"
		for d in dupes:
			print d
		print "\n"
	# return list(map(lambda x: (x.split())[1], dupes))
	return markers

def createPED(snp_map, report, output_name):
	check_file(snp_map, "SNP Map")
	check_file(report, "Illumina Genotyping Report")
	create_sorted_SNP_Map(snp_map)
	markers = get_SNP_Map_marker_ids()
	print "Generating PED file..."
	# walk through report file
	with open(report, 'r') as report_file, open(output_name, 'w') as output_file:
		# find the 'data' start
		data_start = subprocess.check_output(['grep -n "\[Data\]" ' + report], shell=True)
		# get the [Data] heading
		start_pos = int((data_start.split(':'))[0])
		current_sample = ''
		sample = {}
		#write header
		output_file.write("#Generated at %s\n#family sample parent1 parent2 sex phenotype" % datetime.datetime.now())
		for marker in markers:
			output_file.write(" %s" % marker)
		output_file.write("\n")
		# loop over each line in report
		for line_num, report_line in enumerate(report_file, 1):
			# get the row after the table heading
			if line_num > start_pos + 1:
				line = report_line.split()
				if line[1] != current_sample:
					if sample:
						output_file.write('test %s 0 0 0 100 ' % current_sample)
						for marker in markers:
							if sample.has_key(marker):
								output_file.write('%s ' % sample[marker])
							else:
								print_error ("key %s not found for %s" % (marker, current_sample))

						output_file.write("\n")

					sample = {}
					current_sample = line[1]

				# set 00's if we find --
				sample[line[0]] = (line[6] if line[6] != '-' else '0') + (line[7] if line[7] != '-' else '0')

	print "Cleaning up tmp files..."
	os.remove('/tmp/tmp_snp_map.txt')
	print "\nDone!"
	print "Your output is located in this directory in the file '%s'" % output_name

# setup arguments / help
parser = argparse.ArgumentParser(description='Program to combine a SNP_Map and Illumina Genotyping Report to PED file')
parser.add_argument('SNP Map', help='SNP_Map file to parse')
parser.add_argument('Genotype Report', help='Illumina Genotyping Report')
parser.add_argument('-o', '--output', default='output.ped', type=str, help='Name of the output PED file')
args = parser.parse_args()

#check we have the appropriate arguments
if len(sys.argv) >= 2:
	var_arguments = vars(args)
	# start loading variants
	createPED(var_arguments['SNP Map'],
			  var_arguments['Genotype Report'],
			  var_arguments['output'])
else:
	print_error("please specify a SNP_Map and Illumina Genotype Report")
