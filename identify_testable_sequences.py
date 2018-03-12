#!/usr/bin/python
p6_list = []
p5_list = []
p4_list = []
p3_list = []
p2_list = []

with open('pep_seq_lists/all_HCV_uncleaved.txt', 'r') as r:
	list_uncleaved = [i.strip() for i in r.readlines()]

with open('pep_seq_lists/all_HCV_cleaved.txt', 'r') as s:
	list_cleaved = [i.strip() for i in s.readlines()]

uncleaved = []
untested = []

for p6 in p6_list:
	for p5 in p5_list:
		for p4 in p4_list:
			for p3 in p3_list:
				for p2 in p2_list:
					seq = p6 + p5 + p4 + p3 + p2 + 'C'
					if seq in list_uncleaved:
						uncleaved.append(seq)
					elif seq not in list_cleaved:
						untested.append(seq)

print "uncleaved:"
for seq in uncleaved:
	print '\t', seq

print "untested:"
for seq in untested:
	print '\t', seq

